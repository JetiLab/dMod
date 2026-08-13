## compile.R -- C/C++ compilation and DLL (de)registration helpers

## Windows-only: temp Makevars (existing user Makevars + `lines`) for R_MAKEVARS_USER.
.compileMakevarsUser <- function(lines) {
  f <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
  if (is.na(f) || !file.exists(f)) {
    cand <- path.expand(c("~/.R/Makevars.ucrt", "~/.R/Makevars.win64",
                          "~/.R/Makevars.win", "~/.R/Makevars"))
    cand <- cand[file.exists(cand)]
    f <- if (length(cand)) cand[1] else NA
  }
  prev <- if (!is.na(f) && file.exists(f)) readLines(f, warn = FALSE) else character()
  mv <- tempfile(fileext = ".mk")
  writeLines(c(prev, lines), mv)
  mv
}



## Largest command string we hand to system()/system2(). Both route the whole
## command through a single `sh -c <cmd>` argument (CreateProcess on Windows),
## so the ceiling is the per-argument limit -- 128 KiB on Linux
## (MAX_ARG_STRLEN), ~32 KiB for a Windows command line -- not ARG_MAX. A model
## with a few thousand conditions blows past it, so we stay well below either.
## `dMod.compile.cmdlimit` lowers the threshold so the archive route can be
## exercised without generating a thousand source files.
.compileCmdLimit <- function() {
  lim <- suppressWarnings(as.integer(getOption("dMod.compile.cmdlimit")))
  if (length(lim) == 1L && !is.na(lim)) return(lim)
  if (.Platform$OS.type == "windows") 24000L else 96000L
}

## Split `x` into chunks whose quoted length stays under `maxChars` and whose
## element count stays at or below `maxN`, so each chunk can be passed to one
## compiler/archiver invocation.
.compileChunks <- function(x, maxChars = .compileCmdLimit(), maxN = 100L) {
  if (!length(x)) return(list())
  n <- nchar(x) + 3L                    # shell quotes plus separator
  grp <- integer(length(x)); g <- 1L; acc <- 0L; cnt <- 0L
  for (i in seq_along(x)) {
    if (cnt >= maxN || (cnt > 0L && acc + n[i] > maxChars)) {
      g <- g + 1L; acc <- 0L; cnt <- 0L
    }
    grp[i] <- g; acc <- acc + n[i]; cnt <- cnt + 1L
  }
  unname(split(x, factor(grp, levels = seq_len(g))))
}

## Linker incantation that pulls *every* member of a static archive into the
## shared object. Without it the linker keeps only members that resolve an
## undefined symbol, and R looks up all model entry points by name at run time
## -- so nothing at all would be kept.
.compileWholeArchive <- function(lib) {
  if (Sys.info()[["sysname"]] == "Darwin")
    paste0("-Wl,-force_load,", shQuote(lib))
  else
    paste("-Wl,--whole-archive", shQuote(lib), "-Wl,--no-whole-archive")
}

## Budget of additional shared objects R can still dyn.load(); R caps the
## number of simultaneously loaded DLLs at R_MAX_NUM_DLLS (614 by default).
.compileDLLBudget <- function() {
  lim <- suppressWarnings(as.integer(Sys.getenv("R_MAX_NUM_DLLS", "614")))
  if (is.na(lim)) lim <- 614L
  max(0L, lim - length(getLoadedDLLs()))
}


#' Compile model-related C/C++ code
#'
#' @description
#' Compiles model objects ([parfn], [obsfn], [prdfn]) related C/C++ files into shared libraries via `R CMD SHLIB`.
#'
#' @details
#' Per-file compile and link flags are taken from the `"compileInfo"`
#' attribute that [odemodel()], [Xs()], [Xf()], [Y()] and [P()] attach to
#' their return values. Each entry carries the source file together with the
#' `compileArgs` and `linkArgs` reported by the backend that produced it
#' (`cOde::funC`, `cppDE::cppODE`, `cppDE::cvode`, ...), so solver-specific
#' libraries reach only the files that need them. Objects without
#' `compileInfo` fall back to modelname-based file discovery in the current
#' working directory.
#'
#' @param ... One or more model objects.
#' @param output Optional name for a combined shared library. When set, all
#'   files are linked into one object and the union of their `linkArgs` is
#'   applied.
#' @param args Additional compiler/linker flags applied to every file.
#' @param cores Parallel compilation jobs (Unix only, requires `cores > 1`).
#' @param chunkSize Maximum number of files handed to a single archiver call
#'   when the combined link goes through a static archive (see Details).
#' @param verbose If `TRUE`, print compiler commands.
#'
#' @section Many conditions:
#' Objects built per condition -- e.g. a [parfn] summed over a few thousand
#' experiments -- carry one source file per condition. Sources are always
#' compiled one file per command, so that side scales, but naming every object
#' on the final `R CMD SHLIB` command line does not: the command travels as a
#' single argument to the shell and overruns the OS limit (128 KiB on Linux,
#' 32 KiB on Windows) somewhere around a thousand files. Past that limit
#' `compile()` bundles the objects into a static archive instead, appended in
#' chunks of `chunkSize`, and links that archive whole into the shared object.
#' The result is identical; only the command lines get shorter. This requires
#' `output` to be set -- without it every source becomes its own shared object
#' and R's `R_MAX_NUM_DLLS` cap (614 by default) is hit first.
#'
#' @return Invisibly `TRUE` on success.
#' @export
compile <- function(..., output = NULL, args = NULL, cores = 1, chunkSize = 100,
                    verbose = FALSE) {

  ## save & restore env
  old <- Sys.getenv(c("PKG_CFLAGS","PKG_CXXFLAGS","PKG_CPPFLAGS","PKG_LIBS"), unset = NA)
  on.exit({
    for (n in names(old))
      if (is.na(old[n])) Sys.unsetenv(n) else Sys.setenv(structure(old[n], names = n))
  }, add = TRUE)

  objs <- list(...)
  if (!length(objs)) stop("No objects")
  obj.names <- as.character(substitute(list(...)))[-1]
  Rbin  <- shQuote(file.path(R.home("bin"), "R"))
  so    <- .Platform$dynlib.ext
  cfg   <- function(x) trimws(system(paste(Rbin, "CMD config", x), intern = TRUE))
  strip <- function(x) trimws(gsub("(^| )-std=[^ ]+", "", x))

  ## classify objects
  is_dmod <- vapply(objs, inherits, logical(1), c("obsfn","parfn","prdfn"))
  is_cpp  <- vapply(objs, function(o) !is.null(attr(o, "srcfile")), logical(1))

  ## Collect per-file build info.
  ## Primary source is `attr(o, "compileInfo")` carrying
  ## (srcfile, compileArgs, linkArgs) as reported by cOde/cppDE/CVODE.
  ## Falls back to modelname-based file discovery for objects that lack
  ## compileInfo, and to the bare `srcfile` attribute for raw cppDE objects.
  info_from_compileInfo <- unlist(
    lapply(objs, function(o) attr(o, "compileInfo")),
    recursive = FALSE
  )

  info_fallback <- list()
  for (i in seq_along(objs)) {
    o <- objs[[i]]
    if (!is.null(attr(o, "compileInfo"))) next
    if (is_dmod[i]) {
      b <- outer(modelname(o), c("","_deriv","_s","_s2","_sdcv","_dfdx","_dfdp"), paste0)
      cand <- c(paste0(b, ".c"), paste0(b, ".cpp"))
      src <- cand[file.exists(cand)]
      for (s in src)
        info_fallback[[length(info_fallback) + 1]] <-
          list(srcfile = normalizePath(s, winslash = "/", mustWork = TRUE),
               compileArgs = "", linkArgs = "")
    } else if (is_cpp[i]) {
      s <- attr(o, "srcfile")
      if (length(s) && nzchar(s) && file.exists(s))
        info_fallback[[length(info_fallback) + 1]] <- list(
          srcfile     = normalizePath(s, winslash = "/", mustWork = TRUE),
          compileArgs = attr(o, "compileArgs") %||% "",
          linkArgs    = attr(o, "linkArgs")    %||% "",
          sparse      = isTRUE(attr(o, "sparse")))
    }
  }

  info <- c(info_from_compileInfo, info_fallback)

  ## Expand entries with multiple srcfiles (e.g. cOde spills _deriv.c
  ## alongside the main .c) into one entry per file.
  info <- unlist(lapply(info, function(e) {
    if (!length(e$srcfile)) return(list())
    if (length(e$srcfile) == 1L) return(list(e))
    lapply(e$srcfile, function(s) list(srcfile = s, compileArgs = e$compileArgs,
                                       linkArgs = e$linkArgs, sparse = e$sparse))
  }), recursive = FALSE)

  info <- Filter(function(e) length(e$srcfile) == 1L && nzchar(e$srcfile) && file.exists(e$srcfile), info)
  if (!length(info)) stop("No source files found")

  ## Deduplicate by srcfile, keeping the first (non-empty) flags we saw.
  ord <- order(vapply(info, function(e) e$srcfile, character(1)))
  info <- info[ord]
  keep <- !duplicated(vapply(info, function(e) e$srcfile, character(1)))
  info <- info[keep]

  files      <- vapply(info, function(e) e$srcfile, character(1))
  roots      <- sub("\\.[^.]+$", "", basename(files))
  roots_full <- sub("\\.[^.]+$", "", files)

  ## compiler flags
  if (.Platform$OS.type == "windows") cores <- 1
  pic  <- if (.Platform$OS.type == "windows") "" else "-fPIC"
  base <- paste("-O2 -DNDEBUG -w", pic)

  ## KLU flags for sparse models, mirroring cppDE::compile(). The flag lives on
  ## the cppDE object inside an odemodel, so a bare `attr(o, "sparse")` on the
  ## dMod fn objects handed to compile() never sees it -- go through the
  ## per-file info. The `-DKLU*` fallback covers objects whose compileInfo was
  ## built before `sparse` was recorded there.
  uses_klu <- any(vapply(objs, function(o) isTRUE(attr(o, "sparse")), logical(1))) ||
    any(vapply(info, function(e) isTRUE(e$sparse) ||
                 grepl("(^|\\s)-DKLU", e$compileArgs %||% ""), logical(1)))
  klu_flag <- ""; klu_lib <- ""
  if (uses_klu) {
    cfgCppDE <- .cppDE_config()
    if (!isTRUE(cfgCppDE$klu_available))
      stop("A sparse Jacobian was requested, but cppDE was installed without the ",
           "KLU linear solver.\n  Install SuiteSparse/KLU and re-install cppDE -- ",
           "see cppDE::install_libs(\"suitesparse\").", call. = FALSE)
    klu_flag <- trimws(paste("-DKLU", cfgCppDE$klu_cflags))
    klu_lib  <- cfgCppDE$klu_libs
  }

  ## shared pieces (compiler/linker) that apply to every file
  cxx_base <- paste(base, klu_flag)
  extra_args <- paste(c(args), collapse = " ")
  if (nzchar(extra_args)) {
    base     <- paste(base,     extra_args)
    cxx_base <- paste(cxx_base, extra_args)
  }
  ## BLAS/LAPACK: on Windows `R CMD config BLAS_LIBS` returns a value with
  ## unexpanded `$(R_HOME)`/`$(R_ARCH)` references. Those go into PKG_LIBS as
  ## an env var, and make should re-expand them, but in practice the
  ## expansion is unreliable inside SHLIB-generated link commands -- the
  ## final g++ invocation comes out without any BLAS libs. We sidestep that
  ## by building an absolute -L path here and skipping `R CMD config`.
  if (.Platform$OS.type == "windows") {
    ## R.home("bin") already resolves to the arch-specific bin dir (…/bin/x64)
    ## on R >= 4.2, which is where Rblas.dll / Rlapack.dll live. Appending
    ## .Platform$r_arch again produced a bogus …/bin/x64/x64 -L path (it only
    ## ever linked by accident, via the default -L…/bin/x64 that -lR adds).
    ## dMod requires R >= 4.5.0, so R.home("bin") is always the correct dir.
    r_bin   <- R.home("bin")
    blaslapack <- paste0("-L", shQuote(r_bin), " -lRlapack -lRblas")
  } else {
    blaslapack <- paste(cfg("LAPACK_LIBS"), cfg("BLAS_LIBS"))
  }
  base_libs <- paste(klu_lib, blaslapack)
  cppflags  <- paste0("-I", system.file("include", package = "cppDE"))

  ## Compiler invocation bits cached up front so parallel forks don't each
  ## re-spawn R-CMD-config. Used by compile_one_obj() for the direct
  ## $CC/$CXX -c path.
  cc_bin      <- cfg("CC")
  cxx_bin     <- cfg("CXX")
  cflags_R    <- cfg("CFLAGS")
  cxxflags_R  <- cfg("CXXFLAGS")
  cpicflags   <- cfg("CPICFLAGS")
  cxxpicflags <- cfg("CXXPICFLAGS")
  r_inc       <- paste0("-I", shQuote(R.home("include")))

  ## toolchain report (use a representative entry for display)
  if (any(grepl("\\.c$",   files))) cat(sprintf("using C compiler:   %s [%s]\n", strip(cc_bin),  trimws(base)))
  if (any(grepl("\\.cpp$", files))) cat(sprintf("using C++ compiler: %s [%s]\n", strip(cxx_bin), trimws(cxx_base)))

  ## unload stale DLLs
  loaded <- getLoadedDLLs()
  for (i in seq_along(roots))
    if (roots[i] %in% names(loaded)) try(dyn.unload(loaded[[roots[i]]][["path"]]), silent = TRUE)
  if (!is.null(output)) try(dyn.unload(paste0(output, so)), silent = TRUE)

  ## Compile one file with its own compile/link flags applied via PKG_*.
  ## Each invocation sets the env just before shelling out to R CMD SHLIB,
  ## so per-file linkArgs (e.g. Sundials libs for CVODE) reach only the
  ## files that need them. Works inside mclapply because each fork has its
  ## own env.
  compile_one <- function(entry) {
    extra_c <- entry$compileArgs %||% ""
    pkg_c  <- trimws(paste(base,     extra_c))
    pkg_cx <- trimws(paste(cxx_base, extra_c))
    pkg_l  <- trimws(paste(base_libs, entry$linkArgs %||% ""))
    Sys.setenv(
      PKG_CFLAGS   = pkg_c,
      PKG_CXXFLAGS = pkg_cx,
      PKG_CPPFLAGS = cppflags,
      PKG_LIBS     = pkg_l
    )
    if (.Platform$OS.type == "windows") {
      mv <- .compileMakevarsUser(c(
        paste("PKG_CFLAGS =",   pkg_c),
        paste("PKG_CXXFLAGS =", pkg_cx),
        paste("PKG_CPPFLAGS =", cppflags),
        paste("PKG_LIBS =",     pkg_l)
      ))
      old_mu <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
      Sys.setenv(R_MAKEVARS_USER = mv)
      on.exit({
        if (is.na(old_mu)) Sys.unsetenv("R_MAKEVARS_USER")
        else Sys.setenv(R_MAKEVARS_USER = old_mu)
        unlink(mv)
      }, add = TRUE)
    }
    cmd <- paste(Rbin, "CMD SHLIB", shQuote(entry$srcfile))
    if (verbose) cat(cmd, "\n")
    if (system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose) != 0)
      stop("Compilation failed: ", entry$srcfile)
  }

  ## Compile a single source to a .o object only, via a direct $CC/$CXX -c
  ## invocation. Used by the combined-output path so the per-file compile
  ## phase can run in parallel; the subsequent R CMD SHLIB link then sees
  ## the .o files are up-to-date and skips recompilation.
  compile_one_obj <- function(entry) {
    src     <- entry$srcfile
    extra_c <- entry$compileArgs %||% ""
    is_cpp  <- grepl("\\.cpp$", src, ignore.case = TRUE)
    obj     <- sub("\\.[^.]+$", ".o", src)

    if (is_cpp) {
      cmd <- paste(
        cxx_bin, r_inc, cppflags,
        trimws(paste(cxx_base, extra_c)),
        cxxpicflags, cxxflags_R,
        "-c", shQuote(src),
        "-o", shQuote(obj)
      )
    } else {
      cmd <- paste(
        cc_bin, r_inc, cppflags,
        trimws(paste(base, extra_c)),
        cpicflags, cflags_R,
        "-c", shQuote(src),
        "-o", shQuote(obj)
      )
    }

    if (verbose) cat(cmd, "\n")
    if (system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose) != 0)
      stop("Compilation failed: ", src)
    obj
  }

  if (is.null(output)) {
    ## One shared object per source, all of them dyn.load()ed. R caps how many
    ## DLLs it can hold, so refuse up front with an actionable message instead
    ## of failing halfway through a long build.
    budget <- .compileDLLBudget()
    if (length(info) > budget)
      stop(length(info), " source files would need as many shared objects, but R can ",
           "load at most ", budget, " more (R_MAX_NUM_DLLS).\n  Pass output = <name> ",
           "to link them into a single shared object.", call. = FALSE)
    if (.Platform$OS.type == "unix" && cores > 1)
      parallel::mclapply(info, compile_one, mc.cores = cores)
    else for (e in info) compile_one(e)
    for (r in roots_full) dyn.load(paste0(r, so))
  } else {
    ## Combined output: per-file compile to .o (parallel on Unix when cores>1,
    ## serial otherwise -- including on Windows), then a single R CMD SHLIB
    ## link over the original sources. Because every .o is freshly written
    ## above, make sees them as up-to-date and only runs the link recipe;
    ## passing the source list lets SHLIB pick the C++ linker when any
    ## source is .cpp, which a .o-only invocation would miss. The pre-compile
    ## also has to run on Windows: the single-call SHLIB (compile + link in
    ## one go) was occasionally producing .dll files that LoadLibrary
    ## couldn't resolve when the source pulled in BLAS via the symbolic-
    ## mode chain wrapper -- splitting compile and link sidesteps that.
    if (.Platform$OS.type == "unix" && cores > 1)
      parallel::mclapply(info, compile_one_obj, mc.cores = cores)
    else
      for (e in info) compile_one_obj(e)

    ## Link step: union of every entry's linkArgs (dedup) so Sundials-dependent
    ## files still pull their libs.
    all_link <- unique(unlist(lapply(info, function(e) strsplit(trimws(e$linkArgs %||% ""), "\\s+")[[1]])))
    all_link <- all_link[nzchar(all_link)]
    all_compile <- unique(unlist(lapply(info, function(e) strsplit(trimws(e$compileArgs %||% ""), "\\s+")[[1]])))
    all_compile <- all_compile[nzchar(all_compile)]

    output <- sub(paste0("\\", so, "$"), "", output)

    ## Link inputs. Naming every source on the SHLIB command line is the direct
    ## route, but that command is passed to the shell as one argument, so a few
    ## thousand conditions overrun the per-argument limit. Past that point the
    ## objects go into a static archive -- appended in chunks so no single `ar`
    ## call overruns it either -- and SHLIB gets one anchor source plus the
    ## archive. Anchor on a C++ source whenever there is one: SHLIB picks the
    ## C++ linker (and libstdc++) from the sources it is handed, never from the
    ## archive contents.
    link_files <- files
    if (nchar(paste(shQuote(files), collapse = " ")) > .compileCmdLimit()) {
      is_cxx  <- grepl("\\.cpp$", files, ignore.case = TRUE)
      anchor  <- if (any(is_cxx)) which(is_cxx)[1] else 1L
      objects <- sub("\\.[^.]+$", ".o", files)
      lib     <- file.path(dirname(files[1]), paste0(output, "_objects.a"))
      unlink(lib)
      on.exit(try(unlink(lib), silent = TRUE), add = TRUE)
      ## Route the archiver through `R CMD` rather than calling it directly:
      ## on Windows that is what puts the Rtools toolchain on PATH, exactly as
      ## it does for the R CMD SHLIB call below.
      ar_bin     <- cfg("AR");     if (!nzchar(ar_bin))     ar_bin     <- "ar"
      ranlib_bin <- cfg("RANLIB"); if (!nzchar(ranlib_bin)) ranlib_bin <- "ranlib"
      chunks <- .compileChunks(objects[-anchor], maxN = as.integer(chunkSize))
      for (i in seq_along(chunks)) {
        ## `q` appends without the index; the index is written once by ranlib.
        cmd <- paste(Rbin, "CMD", ar_bin, if (i == 1L) "qc" else "q", shQuote(lib),
                     paste(shQuote(chunks[[i]]), collapse = " "))
        if (verbose) cat(cmd, "\n")
        if (system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose) != 0)
          stop("Archiving failed at chunk ", i, " of ", length(chunks), call. = FALSE)
      }
      if (system(paste(Rbin, "CMD", ranlib_bin, shQuote(lib)),
                 ignore.stdout = !verbose, ignore.stderr = !verbose) != 0)
        stop("Building the archive index failed: ", lib, call. = FALSE)
      cat(sprintf("archived %d objects into %s (%d chunks)\n",
                  length(objects) - 1L, basename(lib), length(chunks)))
      link_files <- files[anchor]
      base_libs  <- paste(.compileWholeArchive(lib), base_libs)
    }

    pkg_cflags   <- trimws(paste(base,     paste(all_compile, collapse = " ")))
    pkg_cxxflags <- trimws(paste(cxx_base, paste(all_compile, collapse = " ")))
    pkg_libs     <- trimws(paste(base_libs, paste(all_link, collapse = " ")))
    Sys.setenv(
      PKG_CFLAGS   = pkg_cflags,
      PKG_CXXFLAGS = pkg_cxxflags,
      PKG_CPPFLAGS = cppflags,
      PKG_LIBS     = pkg_libs
    )

    ## Belt-and-suspenders for Windows: env-imported PKG_LIBS has been
    ## observed to vanish from SHLIB's generated link command on some R/rtools
    ## combinations, leaving the .dll unlinked against BLAS/LAPACK. Drop a
    ## per-link Makevars(.win) alongside the source files so make picks it up
    ## even if the environment doesn't make it through. We clean it up after
    ## the link so the directory state stays hermetic.
    mv_dir  <- dirname(files[1])
    mv_name <- if (.Platform$OS.type == "windows") "Makevars.win" else "Makevars"
    mv_path <- file.path(mv_dir, mv_name)
    mv_pre  <- if (file.exists(mv_path)) readLines(mv_path, warn = FALSE) else NULL
    writeLines(c(
      paste("PKG_CFLAGS =",   pkg_cflags),
      paste("PKG_CXXFLAGS =", pkg_cxxflags),
      paste("PKG_CPPFLAGS =", cppflags),
      paste("PKG_LIBS =",     pkg_libs)
    ), mv_path)
    on.exit({
      if (is.null(mv_pre)) try(unlink(mv_path), silent = TRUE)
      else                 try(writeLines(mv_pre, mv_path), silent = TRUE)
    }, add = TRUE)

    ## Windows fallback for BLAS/LAPACK: inject PKG_* via R_MAKEVARS_USER.
    if (.Platform$OS.type == "windows") {
      mv <- .compileMakevarsUser(c(
        paste("PKG_CFLAGS =",   pkg_cflags),
        paste("PKG_CXXFLAGS =", pkg_cxxflags),
        paste("PKG_CPPFLAGS =", cppflags),
        paste("PKG_LIBS =",     pkg_libs)
      ))
      old_mu <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
      Sys.setenv(R_MAKEVARS_USER = mv)
      on.exit({
        if (is.na(old_mu)) Sys.unsetenv("R_MAKEVARS_USER")
        else Sys.setenv(R_MAKEVARS_USER = old_mu)
        unlink(mv)
      }, add = TRUE)
    }

    out <- file.path(dirname(files[1]), paste0(output, so))
    try(dyn.unload(out), silent = TRUE)
    if (file.exists(out)) unlink(out)
    ## Link, capturing stdout+stderr via system2() pipes. Strip the compiler
    ## banner afterwards: the .o files are already fresh from compile_one_obj,
    ## so make only runs the link recipe and "using C/C++ compiler:" would be
    ## misleading. We pass the streams through system2(stdout/stderr = TRUE)
    ## rather than appending a `2>&1` token to a command string: on Windows
    ## that trailing token is not consumed by a shell but swallowed by
    ## R CMD SHLIB as the make override `PKG_LIBS=2>&1`, which beats every
    ## Makevars/R_MAKEVARS_USER assignment and strips the BLAS/LAPACK libs
    ## (breaking the symbolic-mode chain_jac link with "undefined reference to
    ## dgemm_"). system2() keeps the argument vector clean and is identical on
    ## Linux/macOS, where the shell would have consumed the redirection anyway.
    Rexe <- file.path(R.home("bin"), "R")
    shlib_args <- c("CMD", "SHLIB", shQuote(link_files), "-o", shQuote(out))
    if (verbose) cat(shQuote(Rexe), paste(shlib_args, collapse = " "), "\n")
    out_lines <- suppressWarnings(
      system2(Rexe, shlib_args, stdout = TRUE, stderr = TRUE)
    )
    status <- attr(out_lines, "status")
    if (verbose) {
      out_lines <- out_lines[!grepl("^using (C|C\\+\\+) compiler:", out_lines)]
      writeLines(out_lines)
    }
    if (!is.null(status) && status != 0L)
      stop("Compilation failed:\n", paste(out_lines, collapse = "\n"))
    if (!file.exists(out))
      stop("R CMD SHLIB returned exit 0 but did not produce ", out, ":\n",
           paste(out_lines, collapse = "\n"))
    dyn.load(out)
    for (i in which(is_dmod))
      eval.parent(parse(text = paste0("modelname(", obj.names[i], ") <- '", output, "'")))
  }

  invisible(TRUE)
}




#' Determine loaded DLLs available in working directory
#' 
#' @return Character vector with the names of the loaded DLLs available in the working directory
#' @export
getLocalDLLs <- function() {
  
  all.dlls <- getLoadedDLLs()
  is.local <- sapply(all.dlls, function(x) grepl(getwd(), unclass(x)$path, fixed = TRUE))
  names(is.local)[is.local]
  
}




#' Load shared object for a dMod object
#' 
#' Usually when restarting the R session, although all objects are saved in
#' the workspace, the dynamic libraries are not linked any more. `loadDLL`
#' is a wrapper for `dyn.load` that uses the "modelname" attribute of
#' dMod objects like prediction functions, observation functions, etc. to
#' load the corresponding shared object.
#' 
#' @param ... objects of class prdfn, obsfn, parfn, objfn, ...
#' 
#' @export
loadDLL <- function(...) {
  
  .so <- .Platform$dynlib.ext
  models <- modelname(...)
  files <- paste0(outer(models, c("", "_s", "_s2", "_sdcv", "_deriv", "_dfdx", "_dfdp"), paste0), .so)
  files <- files[file.exists(files)]
  
  for (f in files) {
    try(dyn.unload(f), silent = TRUE)
    dyn.load(f)
  }
  message("The following local files were dynamically loaded: ", paste(files, collapse = ", "))
}


## compileInfo plumbing (moved from classes.R) ----------------------------------------

## ODE model class -------------------------------------------------------------------

## Collect per-sub-object build info from ODE model pieces.
## Each backend (cOde::funC, cppDE::cppODE, cppDE::cvode) may attach
## `srcfile`, `compileArgs`, and `linkArgs` to its func/extended result. For
## backends that don't (cOde), we fall back to modelname-based file discovery
## in the current working directory. The resulting list is the single
## authoritative source consulted by `compile()` when given dMod fn objects.

## Per-entry dedup key, cached on the list itself. Summing a parfn over a few
## thousand conditions merges pairwise, so recomputing every key on every merge
## would make building the model quadratic in the number of conditions.
.compileInfoKeys <- function(x) {
  k <- attr(x, "srckeys")
  if (!is.null(k) && length(k) == length(x)) return(k)
  vapply(x, function(e) if (length(e$srcfile)) paste(e$srcfile, collapse = "\x1f")
                        else NA_character_, character(1))
}

## Merge two compileInfo lists, deduplicating by srcfile (per file, the first
## occurrence wins -- that keeps the originating compile/link flags). Returns
## NULL when both inputs are empty so the attribute stays absent on objects
## that never had native code to begin with.
.mergeCompileInfo <- function(a, b) {
  if (!length(a) && !length(b)) return(NULL)
  ka <- .compileInfoKeys(a)
  kb <- .compileInfoKeys(b)
  keep_a <- !is.na(ka) & !duplicated(ka)
  keep_b <- !is.na(kb) & !duplicated(kb) & !(kb %in% ka[keep_a])
  out  <- c(a[keep_a], b[keep_b])
  if (!length(out)) return(NULL)
  attr(out, "srckeys") <- c(ka[keep_a], kb[keep_b])
  out
}

.collectCompileInfo <- function(...) {
  objs <- list(...)
  objs <- objs[!vapply(objs, is.null, logical(1))]
  out <- list()
  for (o in objs) {
    src <- attr(o, "srcfile")
    if (is.null(src) || !length(src) || !nzchar(src)) {
      mname <- attr(o, "modelname")
      if (is.null(mname) && is.character(o)) mname <- unname(o[1])
      if (is.null(mname) || !nzchar(mname)) next
      b <- outer(mname, c("", "_deriv", "_s", "_s2", "_sdcv", "_dfdx", "_dfdp"), paste0)
      cand <- c(paste0(b, ".c"), paste0(b, ".cpp"))
      src <- cand[file.exists(cand)]
      if (!length(src)) next
      src <- normalizePath(src, winslash = "/", mustWork = FALSE)
    } else {
      src <- normalizePath(src, winslash = "/", mustWork = FALSE)
    }
    out[[length(out) + 1]] <- list(
      srcfile     = src,
      compileArgs = attr(o, "compileArgs") %||% "",
      linkArgs    = attr(o, "linkArgs")    %||% "",
      ## cppDE flags a sparse-Jacobian model here; the flag has to travel with
      ## the file because `compile()` only ever sees the wrapping dMod fn.
      sparse      = isTRUE(attr(o, "sparse"))
    )
  }
  out
}


## odemodel() constructor lives in R/odeClass.R.


