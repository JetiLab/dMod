#' Detect number of free cores
#'
#' @description Estimates free cores from the 1-min load average.
#' Supports Linux, macOS, and remote machines via SSH.
#' On Windows, returns 1 with a warning (no load average available).
#' Result is floored at 1 so it can be fed straight into
#' `mclapply(mc.cores = ...)` without crashing under heavy load.
#' @param machine character vector of SSH hosts, e.g. "user@@localhost".
#' NULL (default) for the local machine.
#' @return numeric vector of free cores (>= 1) with attributes "ncores" and "used".
#' @export
detectFreeCores <- function(machine = NULL) {
  
  .getLoadAndCores <- function(prefix = NULL) {
    cmd <- function(x) {
      if (!is.null(prefix)) x <- paste("ssh", prefix, x)
      system(x, intern = TRUE)
    }
    
    # Detect OS: use uname for remote, Sys.info() for local
    os <- if (!is.null(prefix)) cmd("uname") else Sys.info()[["sysname"]]
    
    if (grepl("Windows", os, ignore.case = TRUE)) {
      # No load average on Windows -- return 1 free core as safe default
      warning("detectFreeCores: load average not available on Windows, returning 1")
      nCores <- parallel::detectCores()
      return(list(free = 1, nCores = nCores, occupied = NA_real_))
    }
    
    if (grepl("Darwin", os, ignore.case = TRUE)) {
      # macOS: sysctl provides load avg as "{ x.xx x.xx x.xx }"
      occupied <- as.numeric(strsplit(cmd("sysctl -n vm.loadavg"), " ")[[1]][2])
      nCores <- as.numeric(cmd("sysctl -n hw.ncpu"))
    } else {
      # Linux: read 1-min load average from /proc/loadavg
      occupied <- as.numeric(strsplit(cmd("cat /proc/loadavg"), " ", fixed = TRUE)[[1]][1])
      nCores <- as.numeric(cmd("nproc --all"))
    }
    
    # Floor at 1: callers feed `free` straight into mclapply(mc.cores = ...)
    # which rejects 0. Under heavy load the 1-min average can exceed nCores;
    # reporting 1 means "serialise, don't die" instead of crashing.
    list(free = max(1L, round(nCores - occupied)), nCores = nCores, occupied = occupied)
  }
  
  if (!is.null(machine)) {
    output <- lapply(machine, .getLoadAndCores)
    freeCores <- vapply(output, `[[`, numeric(1), "free")
    attr(freeCores, "ncores") <- vapply(output, `[[`, numeric(1), "nCores")
    attr(freeCores, "used") <- vapply(output, `[[`, numeric(1), "occupied")
  } else {
    res <- .getLoadAndCores()
    freeCores <- res$free
    attr(freeCores, "ncores") <- res$nCores
    attr(freeCores, "used") <- res$occupied
  }

  # CRAN policy: R CMD check sets _R_CHECK_LIMIT_CORES_ and parallel forbids
  # mc.cores > 2 in that mode. Cap before returning so all mclapply callsites
  # (P, normL2, nlme, compile) stay legal under check without each having to
  # know about the env var.
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(chk) && chk != "false") {
    if (length(freeCores) > 0L) freeCores[] <- pmin(freeCores, 2L)
  }

  freeCores
}


#' Run an R expression in the background (only on UNIX)
#' 
#' @description Generate an R code of the expression that is copied via `scp`
#' to any machine (ssh-key needed). Then collect the results.
#' @details `runbg()` generates a workspace from the `input` argument
#' and copies the workspace to the remote machines via `scp`. This will only
#' work if *an ssh-key had been generated and added to the authorized keys
#' on the remote machine*. The code snippet, i.e. the `...` argument, can
#' include several intermediate results but only the last call which is not
#' redirected into a variable is returned via the variable `.runbgOutput`,
#' see example below.
#'
#' Depending on the `compile` and `link` arguments, build-related files are
#' handled as follows:
#' \itemize{
#'   \item `compile = TRUE`: C/C++ source files are transferred and compiled
#'         remotely via `R CMD SHLIB`.
#'   \item `link = TRUE`: Pre-compiled object files (`.o`) are transferred and
#'         linked remotely.
#'   \item Both `FALSE` (default): Existing shared objects (`.so`) are copied
#'         directly.
#' }
#' When `compile` or `link` is `TRUE`, objects of class `obsfn`, `parfn`, or
#' `prdfn` in the workspace automatically get their `modelname` updated to
#' point to the newly built shared object.
#' @param ... Some R code.
#' @param machine Character vector, e.g. `"localhost"` or `"knecht1.fdm.uni-freiburg.de"`
#' or `c("localhost", "localhost")`.
#' @param filename Character, defining the filename of the temporary file. A random
#' file name is chosen if `NULL`.
#' @param input Character vector, the objects in the workspace that are stored
#' into an R data file and copied to the remote machine.
#' @param compile Logical. If `TRUE`, C/C++ source files (`.c`, `.cpp`) are
#' transferred to the remote machine and fully recompiled into a shared object
#' (`.so`). If set to `TRUE`, this overrides `link = TRUE`.
#' @param link Logical. If `TRUE`, only existing object files (`.o`) are
#' transferred to the remote machine and linked into a shared object (`.so`),
#' skipping compilation. If no `.o` files are found, an error is raised.
#' This option is ignored if `compile = TRUE`.
#' @param wait Logical. Wait until executed. If `TRUE`, the code checks if the
#' result file is already present in which case it is loaded. If not present,
#' `runbg()` starts, produces the result and loads it as `.runbgOutput` directly
#' into the workspace. If `wait = FALSE`, `runbg()` starts in the background
#' and the result is only loaded into the workspace when the `get()` function
#' is called, see Value section.
#' @param recover Logical. This option is useful to recover the functions
#' `check()`, `get()`, `purge()` and `terminate()`, e.g. when a session has
#' crashed. Then, the functions are recreated without restarting the job. They
#' can then be used to get the results of a job without having to do it manually.
#' Requires the correct filename, so if the previous `runbg()` was run with
#' `filename = NULL`, you have to specify the filename manually.
#' @param walltime Optional character. Maximum runtime in the format `"HH:MM:SS"`.
#' If exceeded, the job will be terminated.
#' @return List of functions `check()`, `get()`, `purge()` and `terminate()`. 
#' `check()` checks if the result is ready.
#' `get()` copies the result file
#' to the working directory and loads it into the workspace as an object called `.runbgOutput`. 
#' This object is a list named according to the machines that contains the results returned by each
#' machine.
#' `purge()` deletes the temporary folder
#' from the working directory and the remote machines.
#' `terminate()` kills all running processes associated with this job on the remote machines.
#' @export
#' @examples
#' \dontrun{
#' out_job1 <- runbg({
#'          M <- matrix(rnorm(1e2), 10, 10)
#'          solve(M)
#'          }, machine = c("localhost", "localhost"), filename = "job1")
#' out_job1$check()          
#' out_job1$get()
#' result <- .runbgOutput
#' print(result)
#' out_job1$purge()
#' }
#' \dontrun{
#' # Recover a runbg job with the option "recover"
#' out_job1 <- runbg({
#'          M <- matrix(rnorm(1e2), 10, 10)
#'          solve(M)
#'          }, machine = c("localhost", "localhost"), filename = "job1")
#' Sys.sleep(1)
#' remove(out_job1)
#' try(out_job1$check())
#' out_job1 <- runbg({
#'   "This code is not run"
#' }, machine = c("localhost", "localhost"), filename = "job1", recover = TRUE)
#' out_job1$get()
#' result <- .runbgOutput
#' print(result)
#' out_job1$purge()
#' }
runbg <- function(..., machine = "localhost", filename = NULL, input = ls(.GlobalEnv), compile = FALSE, link = FALSE, wait = FALSE, recover = FALSE, walltime = NULL) {
  
  expr <- as.expression(substitute(...))
  nmachines <- length(machine)
  
  # compile takes precedence over link
  if (compile) link <- FALSE
  
  # Generate a random filename if none is provided
  if (is.null(filename))
    filename <- paste0("tmp_", paste(sample(c(0:9, letters), 5, replace = TRUE), collapse = ""))
  
  filename0 <- filename
  filename <- paste(filename, 1:nmachines, sep = "_")
  
  # Initialize output list
  out <- structure(vector("list", 4), names = c("check", "get", "purge", "terminate"))
  
  # Check whether results are ready on all machines
  out[[1]] <- function() {
    
    check.out <- sapply(1:nmachines, function(m) length(suppressWarnings(
      system(paste0("ssh ", machine[m], " ls ", filename[m], "_folder/ | grep -x ", filename[m], "_result.RData"), 
             intern = TRUE))))
    
    if (all(check.out > 0)) {
      cat("Result is ready!\n")
      return(TRUE)
    } else if (any(check.out > 0)) {
      cat("Result from machines", paste(which(check.out > 0), collapse = ", "), "are ready.")
      return(FALSE)
    } else {
      cat("Not ready!\n") 
      return(FALSE)
    }
    
  }
  
  # Fetch result files from remote machines and load into workspace
  out[[2]] <- function() {
    
    result <- structure(vector(mode = "list", length = nmachines), names = machine)
    for (m in 1:nmachines) {
      .runbgOutput <- NULL
      system(paste0("scp ", machine[m], ":", filename[m], "_folder/", filename[m], "_result.RData ./"), ignore.stdout = TRUE, ignore.stderr = TRUE)
      check <- try(load(file = paste0(filename[m], "_result.RData")), silent = TRUE) 
      if (!inherits(check, "try-error")) result[[m]] <- .runbgOutput
    }
    
    .GlobalEnv$.runbgOutput <- result
    
  }
  
  # Remove temporary folders and files on remote machines and locally
  out[[3]] <- function() {
    
    for (m in 1:nmachines) {
      folder_exists <- suppressWarnings(
        system(paste0("ssh ", machine[m], " '[ -d ", filename[m], "_folder ] && echo 1 || echo 0'"), 
               intern = TRUE)
      )
      if (folder_exists == "1") {
        system(paste0("ssh ", machine[m], " rm -r ", filename[m], "_folder"))
      }
      
      rout_exists <- suppressWarnings(
        system(paste0("ssh ", machine[m], " '[ -f ", filename[m], ".Rout ] && echo 1 || echo 0'"), 
               intern = TRUE)
      )
      if (rout_exists == "1") {
        system(paste0("ssh ", machine[m], " rm ", filename[m], ".Rout"))
      }
    }
    
    local_files <- list.files(pattern = paste0(filename0, ".*"))
    if (length(local_files) > 0) {
      system(paste0("rm ", filename0, "*"))
    }
  }
  
  # Kill all processes associated with this job on the remote machines
  out[[4]] <- function() {
    for (m in 1:nmachines) {
      pids <- suppressWarnings(
        system(paste0("ssh ", machine[m], 
                      " 'ps aux | grep \"", filename[m], "\" | grep -v grep'"), 
               intern = TRUE)
      )
      
      if (length(pids) > 0) {
        running_pids <- sapply(strsplit(pids, "\\s+"), function(x) {
          if (grepl("R", x[8])) x[2] else NULL
        })
        running_pids <- running_pids[!sapply(running_pids, is.null)]
        
        if (length(running_pids) > 0) {
          system(paste0("ssh ", machine[m], " 'kill ", paste(running_pids, collapse = " "), "'"))
          cat("Terminated", length(running_pids), "running processes on", machine[m], "\n")
        } else {
          cat("No running processes found on", machine[m], "\n")
        }
      }
    }
  }
  
  # Recover control functions without re-submitting the job
  if (recover) return(out)
  
  # If result files already exist locally and wait = TRUE, load them directly
  resultfile <- paste(filename, "result.RData", sep = "_")
  if (all(file.exists(resultfile)) & wait) {
    
    result <- structure(vector(mode = "list", length = nmachines), names = machine)
    for (m in 1:nmachines) {
      load(file = resultfile[m])
      result[[m]] <- .runbgOutput
    }
    .GlobalEnv$.runbgOutput <- result
    return(out)
  }
  
  # Save current workspace to be transferred to remote machines
  save(list = input, file = paste0(filename0, ".RData"), envir = .GlobalEnv)
  
  # Collect currently loaded packages to replicate the library state remotely
  pack <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  pack <- pack[!is.na(pack)]
  pack <- paste(paste0("try(library(", pack, "))"), collapse = "\n")
  
  output <- ".runbgOutput"
  
  # Compiler flags mirroring compile() in tools.R
  # Written into a shell script to avoid quoting issues with nested SSH commands
  if (compile || link) {
    
    if (compile) {
      sourcefiles <- paste(
        c(list.files(pattern = glob2rx("*.c")), list.files(pattern = glob2rx("*.cpp"))),
        collapse = " "
      )
    } else {
      object_files <- Sys.glob("*.o")
      if (length(object_files) == 0)
        stop("No .o files found for linking! You must compile first.")
      sourcefiles <- paste(object_files, collapse = " ")
    }
    
  }
  
  if (compile || link) {
    # R code to load the newly built shared object and update modelnames of known function objects
    objfns <- 'obj.fns <- ls()[sapply(ls(), function(nm) inherits(get(nm, envir=.GlobalEnv), c("obsfn", "parfn", "prdfn")))]'
    setmn <- sprintf('for (o in obj.fns) eval(parse(text=paste0("modelname(", o, ") <- \'%s\'")))', paste0(filename0, "_shared_object"))
    load_so <- paste0("dyn.load('", filename0, "_shared_object.so')")
  } else {
    objfns <- NULL
    setmn <- NULL
    load_so <- NULL
  }
  
  # Assemble the R script to be executed on each remote machine
  program <- lapply(1:nmachines, function(m) paste(
    c(
      pack,
      paste0("setwd('~/", filename[m], "_folder')"),
      "rm(list = ls())",
      if (!is.null(walltime)) paste0("Sys.setenv(R_TIMEOUT='", walltime, "')"),
      paste0("load('", filename0, ".RData')"),
      objfns,
      setmn,
      if (!is.null(load_so)) {
        load_so
      } else {
        c("files <- list.files(pattern = '\\\\.so$')", "for (f in files) dyn.load(f)")
      },
      paste0(".node <- ", m),
      if (!is.null(walltime)) {
        paste0(".runbgOutput <- try(tools::pskill(Sys.getpid(), tools::SIGALRM, ", walltime, "); ", as.character(expr), ")")
      } else {
        paste0(".runbgOutput <- try(", as.character(expr), ")")
      },
      paste0("save(", output ,", file = '", filename[m], "_result.RData')")
    ),
    collapse = "\n"
  ))
  
  for (m in 1:nmachines) {
    
    # Write the R script for this machine
    cat(program[[m]], file = paste0(filename[m], ".R"))
    
    # Create remote working directory
    system(paste0("ssh ", machine[m], " mkdir -p ", filename[m], "_folder/"), 
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Clear any leftover files from previous runs
    system(paste0("ssh ", machine[m], " rm -r ", filename[m], "_folder/*"), 
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Transfer workspace
    system(paste0("scp ", getwd(), "/", filename0, ".RData* ", machine[m], ":", filename[m], "_folder/"))
    
    # Transfer R script
    system(paste0("scp ", getwd(), "/", filename[m], ".R* ", machine[m], ":", filename[m], "_folder/"))
    
    # Transfer C/C++ source and object files (always; harmless if none exist)
    system(paste0("scp ", getwd(), "/*.c ", machine[m], ":", filename[m], "_folder/"),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(paste0("scp ", getwd(), "/*.cpp ", machine[m], ":", filename[m], "_folder/"),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(paste0("scp ", getwd(), "/*.o ", machine[m], ":", filename[m], "_folder/"),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Transfer shared objects only when no remote build is requested
    if (!compile && !link) {
      system(paste0("scp ", getwd(), "/*.so ", machine[m], ":", filename[m], "_folder/"),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
    
    # When recompiling, remove stale .o files so R CMD SHLIB starts clean
    if (compile) {
      system(paste0("ssh ", machine[m], " 'rm -f ", filename[m], "_folder/*.o'"),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
    
    # Write and transfer compile script per machine, then build the remote command
    compile_cmd <- ""
    if (compile || link) {
      compile_script_content <- paste(
        "#!/bin/bash",
        "export PKG_CFLAGS='-O2 -DNDEBUG -w -fPIC'",
        "export PKG_CXXFLAGS='-O2 -DNDEBUG -w -fPIC'",
        paste0("export PKG_CPPFLAGS=-I$(Rscript -e ", "'cat(system.file(\"include\", package=\"CppODE\"))'", ")"),
        paste0("cd ", filename[m], "_folder"),
        paste0("R CMD SHLIB ", sourcefiles, " -o ", filename0, "_shared_object.so"),
        sep = "\n"
      )
      compile_script_file <- paste0(filename[m], "_compile.sh")
      cat(compile_script_content, file = compile_script_file)
      system(paste0("scp ", getwd(), "/", compile_script_file, " ", machine[m], ":", filename[m], "_folder/"))
      compile_cmd <- paste0("bash ", filename[m], "_folder/", compile_script_file)
    }
    
    # Execute: compile/link if requested, then run the R script
    # OMP_NUM_THREADS=1 and MKL_NUM_THREADS=1 ensure single-threaded execution per job
    system(paste0(
      "ssh ", machine[m], 
      " 'export OMP_NUM_THREADS=1 && export MKL_NUM_THREADS=1",
      if (nzchar(compile_cmd)) paste0(" && ", compile_cmd),
      " && R CMD BATCH --vanilla ", filename[m], "_folder/", filename[m], ".R'"
    ), intern = FALSE, wait = wait)
  }
  
  if (wait) {
    out$get()
    out$purge()
  } else {
    return(out)
  }
  
}



## ---- HPC/SLURM distributed computing (moved from toolsSeverin.R) ----------

#' Run any R function on a remote HPC system with SLURM
#'
#' @description
#' Generates R and bash scripts, transfers them to a remote HPC system via SSH,
#' and executes the given R code in parallel using the SLURM batch manager.
#' The function handles workspace export, job submission, and result retrieval.
#'
#' @details
#' `distributedComputing()` generates R and bash scripts designed to run
#' on an HPC system managed by SLURM. The current R workspace together with the
#' scripts are exported and transferred to the remote system via SSH.  
#' If ssh-key authentication is not possible, the SSH password can be provided and
#' is used by `sshpass` (which must be installed locally).
#'
#' The code to be executed remotely is passed to the `...` argument; its final
#' output is stored in `cluster_result`, which can be loaded in the local
#' workspace via the `get()` function.
#'
#' It is possible to run multiple repetitions of the same program (via `no_rep`)
#' or to pass a list of parameter arrays through `var_values`. Parameters that
#' vary between runs must be named `var_i`, where *i* matches the index
#' of the corresponding array in `var_values`.
#'
#' @param ... R code to be remotely executed. Parameters to be changed for each run
#' must be named `var_i` (see "Details").
#' @param jobname Unique name (character) for the run. Existing runs with the same
#' name will be overwritten. Must not contain the string "Minus".
#' @param partition SLURM partition name to use. Default is `"single"`.
#' @param cores Number of cores per node. Values above 16 may limit available nodes.
#' @param nodes Number of nodes per task. Default is 1; typically should not be changed.
#' @param mem_per_core Memory per CPU core in GB. Default is 2 GB.
#' @param walltime Maximum runtime in format `"hh:mm:ss"`. Default is 1 hour.
#' @param ssh_passwd Password string for SSH authentication via `sshpass`.
#' Optional, and only used if key-based authentication is unavailable.
#' @param machine SSH address of the remote HPC system, e.g. `"user@@cluster"`.
#' @param var_values List of parameter arrays. Each array corresponds to one variable
#' `var_i`. The length of each array determines the number of SLURM array jobs.
#' Mutually exclusive with `no_rep`.
#' @param no_rep Integer number of repetitions (mutually exclusive with `var_values`).
#' @param recover Logical; if `TRUE`, no computation is performed. Instead,
#' the returned list of functions `check()`, `get()`, and `purge()`
#' can be used to interact with previously submitted jobs.
#' @param purge_local Logical; if `TRUE`, the `purge()` function also
#' deletes local result files.
#' @param compile Logical; if `TRUE`, all C/C++ source files (`*.c`, `*.cpp`)
#' are transferred to the cluster and fully recompiled into shared objects (`.so`).
#' If set to `TRUE`, this overrides `link = TRUE`.
#' @param link Logical; if `TRUE`, only existing object files (`*.o`) are
#' transferred to the cluster and linked into shared objects (`.so`),
#' skipping compilation. If no `.o` files are found, an error is raised.
#' This option is ignored if `compile = TRUE`.
#' @param custom_folders Named vector with exactly three elements: `"compiled"`,
#' `"output"`, and `"tmp"`. Each value is a relative path specifying where
#' compiled files, temporary data, and output results should be stored.
#' If `NULL`, all operations occur in the current working directory.
#' @param resetSeeds Logical; if `TRUE` (default), removes `.Random.seed`
#' from the transferred workspace to ensure each node has independent random seeds.
#' @param returnAll Logical; if `TRUE` (default), retrieves all remote files.
#' If `FALSE`, only result files (`*result.RData`) are fetched.
#'
#' @return
#' A list containing three functions:
#' \itemize{
#'   \item `check()` – Checks whether all remote results are complete.
#'   \item `get()` – Downloads results and loads them into
#'         `cluster_result` in the local workspace.
#'   \item `purge()` – Deletes temporary remote files; optionally removes local ones.
#' }
#'
#' @examples
#' \dontrun{
#' outDistributedComputing <- distributedComputing(
#' {
#'   mstrust(
#'     objfun=objective_function,
#'     center=outer_pars,
#'     studyname = "study",
#'     rinit = 1,
#'     rmax = 10,
#'     fits = 48,
#'     cores = 16,
#'     iterlim = 700,
#'     sd = 4
#'   )
#' },
#' jobname = "my_name",
#' partition = "single",
#' cores = 16,
#' nodes = 1,
#' mem_per_core = 2,
#' walltime = "02:00:00",
#' ssh_passwd = "password",
#' machine = "cluster",
#' var_values = NULL,
#' no_rep = 20,
#' recover = F,
#' compile = F,
#' link = F
#' )
#' outDistributedComputing$check()
#' outDistributedComputing$get()
#' outDistributedComputing$purge()
#' result <- cluster_result
#' print(result)
#' 
#' 
#' # calculate profiles
#' var_list <- profileParsPerNode(best_fit, 4)
#' profile_jobname <- paste0(fit_filename,"_profiles_opt")
#' method <- "optimize"
#' profilesDistributedComputing <- distributedComputing(
#'   {
#'     profile(
#'       obj = obj,
#'       pars =  best_fit,
#'       whichPar = (as.numeric(var_1):as.numeric(var_2)),
#'       limits = c(-5, 5),
#'       cores = 16,
#'       method = method,
#'       stepControl = list(
#'         stepsize = 1e-6,
#'         min = 1e-4, 
#'         max = Inf, 
#'         atol = 1e-2,
#'         rtol = 1e-2, 
#'         limit = 100
#'       ),
#'       optControl = list(iterlim = 20)
#'     )
#'   },
#'   jobname = profile_jobname,
#'   partition = "single",
#'   cores = 16,
#'   nodes = 1,
#'   walltime = "02:00:00",
#'   ssh_passwd = "password",
#'   machine = "cluster",
#'   var_values = var_list,
#'   no_rep = NULL,
#'   recover = F,
#'   compile = F,
#'   link = F
#' )
#' profilesDistributedComputing$check()
#' profilesDistributedComputing$get()
#' profilesDistributedComputing$purge()
#' profiles  <- NULL
#' for (i in cluster_result) {
#'   profiles <- rbind(profiles, i)
#' }
#' }
#' 
#' @export
distributedComputing <- function(
    ...,
    jobname,
    partition = "single",
    cores = 16,
    nodes = 1,
    mem_per_core = 2,
    walltime = "01:00:00",
    ssh_passwd = NULL,
    machine = "cluster",
    var_values = NULL,
    no_rep = NULL,
    recover = TRUE,
    purge_local = FALSE,
    compile = FALSE,
    link = FALSE,
    custom_folders = NULL,
    resetSeeds = TRUE,
    returnAll = TRUE
){
  original_wd <- getwd()
  if (is.null(custom_folders)) {
    output_folder_abs <- "./"
  } else if(!is.null(custom_folders) & !all(length(custom_folders) == 3 & sort(names(custom_folders)) == c("compiled", "output", "tmp"))) {
    warning("'custom_folders' must be named vector with exact three elements:\n
            'compiled', 'output', 'tmp', containing relative paths to the resp folders\n
            input is wrong, ignored.\n")
  } else {
    compiled_folder <- custom_folders["compiled"]
    output_folder <- custom_folders["output"]
    tmp_folder <- custom_folders["tmp"]
    
    system(paste0("cp ", compiled_folder, "* ", tmp_folder))
    
    setwd(output_folder)
    output_folder_abs <- getwd()
    
    setwd(original_wd)
    setwd(tmp_folder)
  }
  
  # Safety rule: never compile and link at the same time
  if (compile) link <- FALSE
  
  
  on.exit(setwd(original_wd))
  
  # - definitions - #
  
  # relative path to the working directory, will now allways be used
  
  wd_path <- paste0("./",jobname, "_folder/")
  data_path <- paste0(getwd(),"/",jobname, "_folder/")
  
  # number of repetitions
  if(!is.null(no_rep) & is.null(var_values)) {
    num_nodes <- no_rep - 1
  } else if(is.null(no_rep) & !is.null(var_values)) {
    num_nodes <- length(var_values[[1]]) - 1
  } else {
    stop("I dont know what you want how often done. Please set either 'no_rep' or pass 'var_values' (_not_ both!)")
  }
  
  # define the ssh command depending on 'sshpass' being used
  if(is.null(ssh_passwd)){
    ssh_command <- "ssh "
    scp_command <- "scp "
  } else {
    ssh_command <- paste0("sshpass -p ", ssh_passwd, " ssh ")
    scp_command <- paste0("sshpass -p ", ssh_passwd, " scp ")
  }
  
  # - output functions - #
  # Structure of the output 
  out <- structure(vector("list", 3), names = c("check", "get", "purge"))
  
  # check function
  out[[1]] <- function() {
    
    result_length <- length(
      suppressWarnings(
        system(
          paste0(ssh_command, machine, " 'ls ", jobname, "_folder/ | egrep *result.RData'"),
          intern = TRUE)
      )
    )
    
    if (result_length == num_nodes +1) {
      cat("Result is ready!\n")
      return(TRUE)
    }
    else if (result_length  < num_nodes +1) {
      cat("Result from", result_length, "out of", (num_nodes +1), "nodes are ready.")
      return(FALSE)
    }
    setwd(original_wd)
    
  }
  
  # get function
  out[[2]] <- function () {
    # copy all files back
    if (returnAll == T) {
      system(
        paste0(
          "mkdir -p ", output_folder_abs, "/", jobname,"_folder/results/; ",
          ssh_command, "-n ", machine, # go to remote
          " 'ZSTD_NBTHREADS=0 tar -C ", jobname, "_folder", " -I zstd -cf - ./'", # compress all files on remote)
          " | ", # pipe to local
          "",
          "tar -C ", output_folder_abs, "/", jobname,"_folder/results/ -I zstd -xf -"
        )
      )
    } else {
      # copy only result files back
      system(
        paste0(
          "mkdir -p ", output_folder_abs, "/", jobname,"_folder/results/; ",
          ssh_command, "-n ", machine, " '",
          "ZSTD_NBTHREADS=0 find ", jobname, "_folder -type f -name \"*result.RData\" -exec tar -I zstd -cf - {} +'",
          " | tar --strip-components=1 -I zstd -x -C ", shQuote(paste0(output_folder_abs, "/", jobname,"_folder/results/"))
        )
      )
    }
    
    
    # get list of all currently available output files
    # setwd(paste0(jobname,"_folder/results"))
    result_list <- structure(vector(mode = "list", length = num_nodes+1))
    result_files <- list.files(path=paste0(output_folder_abs,"/",jobname,"_folder/results/"),pattern = glob2rx("*result.RData"))
    # setwd("../../")
    
    # result_files <- Sys.glob(file.path(paste0(wd_path, "/results/*RData")))
    
    for (i in seq(1, length(result_files))) {
      cluster_result <- NULL
      check <- try(load(file = paste0(output_folder_abs,"/",jobname,"_folder/results/",result_files[i])), silent = TRUE) 
      if (!inherits("try-error", check)) result_list[[i]] <- cluster_result
    }
    
    if (length(result_files) != num_nodes +1) {
      cat("\n\tNot all results ready\n")
    }
    setwd(original_wd)
    # results_cluster <- Sys.glob(paste0(wd_path, "*.RData")) %>% map_dfr(load)
    .GlobalEnv$cluster_result <- result_list
  }
  
  
  
  # purge function
  out[[3]] <- function (purge_local = FALSE) {
    # remove files remote
    system(
      paste0(
        ssh_command, machine, " rm -rf ", jobname, "_folder"
      )
    )
    # also remove local files if want so
    if (purge_local) {
      system(
        paste0("rm -rf ", output_folder_abs, "/", jobname,"_folder")
      )
    }
    setwd(original_wd)
  }
  
  # if recover == T, stop here
  if(recover) {
    setwd(original_wd)
    return(out)
  } 
  
  
  
  
  
  
  
  
  
  # - calculation - #
  # create wd for this run
  system(
    paste0("rm -rf ",jobname,"_folder")
  )
  system(
    paste0("mkdir ",jobname,"_folder/")
  )
  
  
  
  # export current workspace
  save.image(file = paste0(wd_path,jobname, "_workspace.RData"), compress = FALSE)
  
  
  # WRITE R
  
  
  # generate list of currently loaded packages
  package_list <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  package_list <- package_list[!is.na(package_list)]
  package_list <- paste(paste0("try(library(", package_list, "))"), collapse = "\n")
  if (compile || link) {
    objfns <- 'obj.fns <- ls()[sapply(ls(), function(nm) inherits(get(nm, envir=.GlobalEnv), c("obsfn", "parfn", "prdfn")))]'
    setmn <- sprintf('for (o in obj.fns) eval(parse(text=paste0("modelname(", o, ") <- \'%s\'")))\n', paste0(jobname, "_shared_object"))
    load_so <- paste0("dyn.load('",jobname,"_shared_object.so')")
  } else {
    objfns <- ""
    setmn <-""
    load_so <- ""
  }
  
  
  
  # generate parameter lists
  if (!is.null(var_values)) {
    var_list <- paste(
      lapply(
        seq(1,length(var_values)),
        function(i) {
          if (is.character(var_values[[i]])) {
            paste0("var_values_", i, "=c('", paste(var_values[[i]], collapse="','"),"')")
          } else {
            paste0("var_values_", i, "=c(", paste(var_values[[i]], collapse=","),")")
          }
          
          
        }
      ),
      collapse = "\n"
    )
    # cat(variable_list)
    
    # List of all names of parameters that will be changes between runs
    var_names <- paste(lapply(seq(1,length(var_values)), function(i) paste0("var_",i)))
    
    # Variables per run
    var_per_run <- paste(
      lapply(
        seq(1, length(var_values)),
        function(i) {
          paste0("var_", i, "=var_values_",i,"[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]")
        }
      ),
      collapse = "\n"
    )
  } else {
    var_list <- ""
    var_per_run <- ""
  }
  
  # cat(var_per_run)
  
  
  # define fixed pars
  fixedpars <- paste(
    "node_ID = Sys.getenv('SLURM_ARRAY_TASK_ID')",
    "job_ID = Sys.getenv('SLURM_JOB_ID')",
    paste0("jobname = ", "'",jobname,"'"),
    sep = "\n"
  )
  
  
  # WRITE R
  expr <- as.expression(substitute(...))
  # expr <- as.expression(substitute(called_function))
  cat(
    paste(
      "#!/usr/bin/env Rscript",
      "",
      "# Load packages",
      package_list,
      "try(library(tidyverse))",
      "",
      "# Load environment",
      paste0("load('",jobname,"_workspace.RData')"),
      "",
      if (resetSeeds == TRUE & exists(".Random.seed")) {
        paste0("# remove random seeds\nrm(.Random.seed)\nset.seed(as.numeric(Sys.getenv('SLURM_JOB_ID')))")
        
      },
      "",
      "# load shared object if precompiled",
      objfns,
      setmn,
      load_so,
      "",
      "# List of variablevalues",
      var_list,
      "",
      "# Define variable values per run",
      var_per_run,
      "",
      "# Fixed parameters",
      fixedpars,
      "",
      "",
      "",
      "# Paste function call",
      paste0("cluster_result <- try(", as.character(expr),")"),
      sep = "\n",
      "save(cluster_result, file = paste0(jobname,'_', node_ID, '_result.RData'))" #'_',job_ID, 
    ),
    file = paste0(wd_path, jobname,".R")
  )
  
  
  
  
  # WRITE BASH
  cat(
    paste(
      "#!/bin/bash",
      "",
      "# Job name",
      paste0("#SBATCH --job-name=",jobname),
      "# Define format of output, deactivated",
      paste0("#SBATCH --output=",jobname,"_%j-%a.out"),
      "# Define format of errorfile, deactivated",
      paste0("#SBATCH --error=",jobname,"_%j-%a.err"),
      "# Define partition",
      paste0("#SBATCH --partition=", partition),
      "# Define number of nodes per task",
      paste0("#SBATCH --nodes=", nodes),
      "# Define number of cores per node",
      paste0("#SBATCH --ntasks-per-node=",cores),
      "# Define walltime",
      paste0("#SBATCH --time=",walltime),
      "# Define of repetition",
      paste0("#SBATCH -a 0-", num_nodes),
      "# memory per CPU core",
      paste0("#SBATCH --mem-per-cpu=", mem_per_core, "gb"),
      "",
      "",
      "# Load compiler modules",
      "module load compiler/gnu/13.3",
      "# Load R modules",
      "module load math/R",
      # paste0("export OPENBLAS_NUM_THREADS=",cores),
      paste0("export OMP_NUM_THREADS=","1"), # paste0("export OMP_NUM_THREADS=",cores),
      paste0("export MKL_NUM_THREADS=", "1"), # paste0("export MKL_NUM_THREADS=",cores),
      "",
      "# Run R script",
      paste0("Rscript ", jobname, ".R"),
      sep = "\n" 
    ),
    file = paste0(wd_path,jobname,".sh")
  )
  
  
  if (compile) {
    # --- FULL RECOMPILATION (.cpp/.c -> .so) ---
    compile_files <- Sys.glob(paste0("*.cpp"))
    compile_files <- append(compile_files, Sys.glob(paste0("*.c")))
    compile_files <- paste(compile_files, collapse = " ")
    
    tar_locale <- paste0("tar -I 'zstd -T0' -cf - ", compile_files, " ", wd_path, "*")
    tar_remote <- paste0("tar -C ./ -I zstd -xf - ; mv -t ./", jobname, "_folder ", compile_files, "; ")
    
    sourcefiles <- paste(
      c(list.files(pattern = glob2rx("*.c")), list.files(pattern = glob2rx("*.cpp"))),
      collapse = " "
    )
    
    compile_remote <- paste0(
      "module load math/R; MAKEFLAGS='-j", cores, "' R CMD SHLIB ",
      sourcefiles, " -o ", jobname, "_shared_object.so; "
    )
    
  } else if (link) {
    # --- LINK ONLY (.o -> .so) ---
    object_files <- Sys.glob("*.o")
    if (length(object_files) == 0)
      stop("No .o files found for linking! You must compile first.")
    
    # Remove any old .so files before linking
    # unlink(list.files(pattern = "(\\.so)$"))
    
    compile_files <- paste(object_files, collapse = " ")
    tar_locale <- paste0("tar -I 'zstd -T0' -cf - ", compile_files, " ", wd_path, "*")
    tar_remote <- paste0("tar -C ./ -I zstd -xf - ; mv -t ./", jobname, "_folder ", compile_files, "; ")
    
    compile_remote <- paste0(
      "module load math/R; R CMD SHLIB ", paste(object_files, collapse = " "),
      " -o ", jobname, "_shared_object.so; "
    )
    
  } else {
    # --- NO BUILD ACTION (.so/.o already available) ---
    compile_files <- Sys.glob(paste0("*.so"))
    compile_files <- append(compile_files, Sys.glob(paste0("*.o")))
    compile_files <- paste(compile_files, collapse = " ")
    
    tar_locale <- paste0("tar -I 'zstd -T0' -cf - ", compile_files, " ", wd_path, "*")
    tar_remote <- paste0("tar -C ./ -I zstd -xf - ; mv -t ./", jobname, "_folder ", compile_files, "; ")
    compile_remote <- ""
  }
  
  ##
  # transfer and run files
  system(
    paste0(
      tar_locale, # compress all files in the local working dir
      " | ", ssh_command, machine, # pipe to ssh session on remote
      " 'if [ -d ", jobname, "_folder ]; then rm -Rf ", jobname,"_folder; fi ;", # remove folder if it exists
      " mkdir -p ", jobname,"_folder; ", # create new wd on remote
      tar_remote, # uncompress files in wd on remote, if necessary move files
      "cd ", jobname, "_folder; ", # change in said wd
      compile_remote, # compile files if said so, if not nothing happen
      "sbatch ", jobname, ".sh'" # start bash script
    )
  )
  
  setwd(original_wd)
  return(out)
}



#' Generate parameter list for distributed profile calculation
#' 
#' @description Generates list of `WhichPar` entries to facillitate distribute
#' profile calculation.
#' @details Lists to split the parameters for which the profiles are calculated
#' on the different nodes.
#' 
#' @param parameters list of parameters 
#' @param fits_per_node numerical, number of parameters that will be send to each node.
#' @param side determine if both sides are calculated (default) or if the profiles are split in 'left' and 'right' for calculation
#' 
#' @return List with two arrays: `from` contains the number of the starting
#' parameter, while `to` stores the respective upper end of the parameter list
#' per node.
#' @examples
#' \dontrun{
#' parameter_list <- setNames(1:10, letters[1:10])
#' var_list <- profileParsPerNode(parameter_list, 4)
#' }
#' 
#' @export
profileParsPerNode <- function(parameters, fits_per_node, side = c("both", "split")[1]) {
  # sanitize side input: must be either "left", "right" or "both"
  if (!(side %in% c("both", "split"))) {
    stop("'side' must be either 'both' or 'split'")
  }
  
  # get the number of parameters
  n_pars <- length(parameters)
  
  # Get number of fits per node
  fits_per_node <- fits_per_node
  
  # determine the number of nodes necessary
  no_nodes <- 1:ceiling(n_pars/fits_per_node)
  
  # generate the lists which parameters are send to which node
  pars_from <- fits_per_node
  pars_to_vec <- fits_per_node
  while (pars_from < (n_pars)) {
    pars_from <- pars_from + fits_per_node
    pars_to_vec <- c(pars_to_vec, pars_from)
  }
  pars_to_vec[length(pars_to_vec)] <- n_pars
  
  pars_from_vec <- c(1, pars_to_vec+1)
  pars_from_vec <- head(pars_from_vec, -1)
  
  # adjust for sides 
  if (side == "both") {
    side_vec <-  rep("both", length(no_nodes))
  } else {
    # split pars_to_vec and pars_from_vec by repeating each element twice
    pars_to_vec <- rep(pars_to_vec, each = 2)
    pars_from_vec <- rep(pars_from_vec, each = 2)
    side_vec <- rep(c("left", "right"), length(no_nodes))
  }
  
  out <- list(from=pars_from_vec, to=pars_to_vec, side = side_vec)
  
  return(out)
}



## Use Julia to calculate steady states -----------------------------------------



## .sanitizeCores (moved from tools.R) ----------------------------------------

.sanitizeCores <- function(cores)  {
  
  max.cores <- parallel::detectCores()
  min(max.cores, cores)
 #  
 # if (Sys.info()[['sysname']] == "Windows") cores <- 1
 # return(cores)
  
}




## .parallelLapply (moved from tools.R) --------------------------------------

# Cross-platform parallel-apply. Unix forks via doParallel; Windows uses a
# PSOCK cluster with explicit library-path + variable export. Driven through
# foreach::%dopar% so the caller body is the same on both platforms.
.parallelLapply <- function(X, FUN, cores = 1L,
                            extraExports = character(0),
                            envir = parent.frame()) {

  cores <- max(1L, as.integer(cores))
  if (cores == 1L)
    return(lapply(X, FUN))

  is_windows <- Sys.info()[["sysname"]] == "Windows"

  if (is_windows) {
    cluster <- parallel::makeCluster(cores)
    on.exit({
      parallel::stopCluster(cluster)
      doParallel::stopImplicitCluster()
    }, add = TRUE)
    doParallel::registerDoParallel(cl = cluster)
    parallel::clusterCall(cl = cluster, function(x) .libPaths(x), .libPaths())
    parallel::clusterEvalQ(cluster, suppressMessages(library(dMod)))
    if (length(extraExports) > 0L)
      parallel::clusterExport(cluster, envir = envir,
                              varlist = extraExports)
  } else {
    doParallel::registerDoParallel(cores = cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
  }

  i <- NULL  # silence R CMD check NSE warning
  loaded_packages <- .packages()
  out <- foreach::foreach(i = seq_along(X),
                          .packages = loaded_packages) %dopar% {
    FUN(X[[i]])
  }
  out
}



`%dopar%` <- foreach::`%dopar%`
