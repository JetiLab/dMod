## Past the shell's argument limit compile() links through a chunked static
## archive instead of naming every object on the SHLIB command line.

test_that(".compileChunks caps both count and command length", {

  files <- sprintf("/a/b/file%04d.o", 1:250)

  chunks <- .compileChunks(files, maxN = 100L)
  expect_equal(lengths(chunks), c(100L, 100L, 50L))
  expect_identical(unlist(chunks, use.names = FALSE), files)

  ## An over-long name forms a chunk of its own rather than an empty one.
  long <- c(strrep("x", 500), files[1:3])
  expect_equal(lengths(.compileChunks(long, maxChars = 100L, maxN = 100L)),
               c(1L, 3L))

  expect_identical(.compileChunks(character(0)), list())

})


test_that(".mergeCompileInfo dedups by srcfile, first flags win", {

  a <- list(list(srcfile = "x.cpp", compileArgs = "A"),
            list(srcfile = "y.cpp", compileArgs = "B"))
  b <- list(list(srcfile = "y.cpp", compileArgs = "OVERRIDDEN"),
            list(srcfile = "z.cpp", compileArgs = "C"),
            list(srcfile = "z.cpp", compileArgs = "D"))

  m <- .mergeCompileInfo(a, b)
  expect_equal(vapply(m, function(e) e$srcfile,     ""), c("x.cpp", "y.cpp", "z.cpp"))
  expect_equal(vapply(m, function(e) e$compileArgs, ""), c("A", "B", "C"))

  ## The cached keys must stay in step with the list they annotate.
  expect_equal(attr(m, "srckeys"), c("x.cpp", "y.cpp", "z.cpp"))
  expect_equal(.compileInfoKeys(m), attr(m, "srckeys"))
  expect_equal(.mergeCompileInfo(m, m), m)

  expect_null(.mergeCompileInfo(NULL, NULL))
  expect_null(.mergeCompileInfo(list(list(srcfile = character(0))), NULL))

})


test_that("compile() links many conditions through a chunked archive", {

  withr::local_dir(tempdir())
  ## Force the archive route without generating a thousand sources.
  withr::local_options(dMod.compile.cmdlimit = 1L)

  conditions <- sprintf("chk%02d", 1:12)
  p <- Reduce("+", lapply(conditions, function(cn)
    P(eqnvec(k1 = "exp(logk1)", A0 = "exp(logA0)*scale"),
      condition = cn, compile = FALSE, modelname = paste0("chunked_p_", cn))))

  expect_length(attr(p, "compileInfo"), length(conditions))
  expect_output(compile(p, output = "chunked_all", cores = 1), "archived 11 objects")

  expect_true(file.exists(paste0("chunked_all", .Platform$dynlib.ext)))
  expect_equal(list.files(pattern = "^chunked_all.*\\.a$"), character(0))
  expect_equal(unique(modelname(p)), "chunked_all")

  ## Every condition's entry point resolves out of the single shared object.
  out <- p(c(logk1 = log(2), logA0 = log(3), scale = 5), deriv = TRUE)
  expect_equal(vapply(conditions, function(cn) unname(out[[cn]]["A0"]), 0),
               setNames(rep(15, length(conditions)), conditions))
  J <- attr(out[[conditions[length(conditions)]]], "deriv")
  expect_equal(J["A0", "scale"], 3)
  expect_equal(J["k1", "logk1"], 2)

})


test_that("compile() refuses more shared objects than R can load", {

  withr::local_dir(tempdir())
  withr::local_envvar(R_MAX_NUM_DLLS = "1")

  p <- P(eqnvec(k1 = "exp(logk1)"), condition = "dll",
         compile = FALSE, modelname = "chunked_budget")

  expect_error(compile(p), "R_MAX_NUM_DLLS")
  expect_error(compile(p), "output = <name>", fixed = TRUE)

})


## The remote build script hits the same limit. These run the generated bash
## locally, no ssh involved.

test_that("the remote build script builds a shared object through the archive", {

  skip_on_os("windows")
  skip_if_not_installed("cppDE")
  withr::local_dir(withr::local_tempdir())

  ## <vector> ahead of the R headers: libc++ has a member `length`, which the
  ## Rf_ remapping macro would swallow.
  for (i in 1:5)
    writeLines(c("#include <vector>", "#include <R.h>", "#include <Rinternals.h>",
                 sprintf('extern "C" SEXP rmt_%d(SEXP x) { std::vector<double> v(2, %d.0); return Rf_ScalarReal(v[0]); }', i, i)),
               sprintf("rmt_%d.cpp", i))
  writeLines(c("#include <R.h>", "#include <Rinternals.h>",
               "SEXP rmt_c(SEXP x) { return ScalarReal(42.0); }"), "rmt_c.c")

  files <- c(list.files(pattern = "\\.c$"), list.files(pattern = "\\.cpp$"))
  writeLines(files, "job_files.txt")

  ## Below the limit the script keeps the plain single-SHLIB route.
  withr::with_options(list(dMod.compile.cmdlimit = 96000L), {
    plain <- .remoteBuildScript(files, output = "plain.so", cxx = TRUE,
                                filelist = "job_files.txt")
    expect_false(grepl("whole-archive", plain))
    expect_match(plain, "R CMD SHLIB rmt_c\\.c", fixed = FALSE)
  })

  withr::local_options(dMod.compile.cmdlimit = 1L)
  script <- .remoteBuildScript(files, output = "remote_all.so", cxx = TRUE,
                               cores = 2, filelist = "job_files.txt")
  ## A C++ anchor, otherwise SHLIB links with the C driver.
  expect_match(script, 'ANCHOR=.?rmt_1\\.cpp')
  expect_match(script, "whole-archive")

  writeLines(script, "build.sh")
  expect_equal(system("bash build.sh", ignore.stdout = TRUE, ignore.stderr = TRUE), 0L)
  expect_true(file.exists("remote_all.so"))
  ## The archive is an intermediate.
  expect_equal(list.files(pattern = "\\.a$"), character(0))

  dyn.load("remote_all.so")
  on.exit(try(dyn.unload("remote_all.so"), silent = TRUE), add = TRUE)
  expect_equal(.Call("rmt_5", 1), 5)
  expect_equal(.Call("rmt_c", 1), 42)

})


test_that("the remote build script links pre-compiled objects the same way", {

  skip_on_os("windows")
  skip_if_not_installed("cppDE")
  withr::local_dir(withr::local_tempdir())

  for (i in 1:4)
    writeLines(c("#include <R.h>", "#include <Rinternals.h>",
                 sprintf('extern "C" SEXP rmo_%d(SEXP x) { return Rf_ScalarReal(%d.0); }', i, i)),
               sprintf("rmo_%d.cpp", i))
  objs <- vapply(list.files(pattern = "\\.cpp$"), function(f) {
    o <- sub("\\.cpp$", ".o", f)
    cc <- trimws(system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX"), stdout = TRUE))
    fl <- trimws(system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXXFLAGS"), stdout = TRUE))
    pic <- trimws(system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXXPICFLAGS"), stdout = TRUE))
    expect_equal(system(paste(cc, shQuote(paste0("-I", R.home("include"))), fl, pic,
                              "-c", shQuote(f), "-o", shQuote(o))), 0L)
    o
  }, character(1), USE.NAMES = FALSE)

  writeLines(objs, "job_files.txt")
  withr::local_options(dMod.compile.cmdlimit = 1L)
  script <- .remoteBuildScript(objs, output = "remote_link.so", link = TRUE,
                               cxx = TRUE, filelist = "job_files.txt")
  writeLines(script, "build.sh")

  expect_equal(system("bash build.sh", ignore.stdout = TRUE, ignore.stderr = TRUE), 0L)
  expect_true(file.exists("remote_link.so"))
  dyn.load("remote_link.so")
  on.exit(try(dyn.unload("remote_link.so"), silent = TRUE), add = TRUE)
  expect_equal(.Call("rmo_4", 1), 4)

})


test_that("compile() reuses objects whose source and command are unchanged", {

  withr::local_dir(withr::local_tempdir())

  trafo <- eqnvec(k1 = "exp(logk1)", A0 = "exp(logA0)*scale")
  conditions <- sprintf("reuse%02d", 1:9)
  mk <- function() Reduce("+", lapply(conditions, function(cn)
    P(trafo, condition = cn, compile = FALSE, modelname = paste0("reuse_p_", cn))))

  p <- mk()
  expect_no_match(capture.output(compile(p, output = "reuse_all", cores = 1)), "reusing")

  ## Codegen rewrites every source, so only the content decides.
  p <- mk()
  expect_output(compile(p, output = "reuse_all", cores = 1),
                paste("reusing", length(conditions), "unchanged"))
  out <- p(c(logk1 = log(2), logA0 = log(3), scale = 5), deriv = TRUE)
  expect_equal(as.numeric(out[[conditions[1]]]["A0"]), 15)
  expect_equal(attr(out[[conditions[9]]], "deriv")["k1", "logk1"], 2)

  ## A changed compile command has to invalidate the entry.
  expect_no_match(capture.output(compile(p, output = "reuse_all", cores = 1,
                                         args = "-DDMOD_TEST_FLAG")),
                  "reusing")

})


test_that(".compilePCHIncludes only accepts prologues it can safely prepend", {

  d <- withr::local_tempdir()
  clean <- file.path(d, "clean.cpp")
  writeLines(c("/** generated **/", "", "#include <cmath>", "#include <vector>",
               "double f() { return 0; }"), clean)
  expect_equal(.compilePCHIncludes(clean), c("#include <cmath>", "#include <vector>"))

  ## A macro ahead of the includes can change how they expand, so no header.
  guarded <- file.path(d, "guarded.cpp")
  writeLines(c("#define _USE_MATH_DEFINES", "#include <cmath>", "double f() { return 0; }"),
             guarded)
  expect_null(.compilePCHIncludes(c(clean, guarded)))
  expect_null(.compilePCHIncludes(file.path(d, "none.cpp") |>
                                    (\(f) { writeLines("int f();", f); f })()))

})
