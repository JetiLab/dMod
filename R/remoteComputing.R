#' Detect number of free cores (on UNIX)
#' 
#' @description Read \code{/proc/loadavg} and subtract from the number of cores
#' @param machine character, e.g. "user@@localhost".
#' @export 
detectFreeCores <- function(machine = NULL) {
  if (!is.null(machine)) {
    output <- lapply(machine, function(m) {
      occupied <- as.numeric(strsplit(system(paste("ssh", m, "cat /proc/loadavg"), intern = TRUE), split = " ", fixed = TRUE)[[1]][1])  
      nCores <- as.numeric(system(paste("ssh", m, "nproc --all"), intern = TRUE))
      free <- max(c(0, round(nCores - occupied)))
      list(free, nCores, occupied)
    })
    freeCores <- unlist(lapply(output, function(o) o[[1]]))
    attr(freeCores, "ncores") <- unlist(lapply(output, function(o) o[[2]]))
    attr(freeCores, "used") <- unlist(lapply(output, function(o) o[[3]]))
  } else {
    occupied <- as.numeric(strsplit(system("cat /proc/loadavg", intern = TRUE), split = " ", fixed = TRUE)[[1]][1])      
    nCores <- as.numeric(system("nproc --all", intern = TRUE))
    freeCores <- max(c(0, round(nCores - occupied)))
    attr(freeCores, "ncores") <- nCores
    attr(freeCores, "used") <- occupied
  }
  return(freeCores)
}


#' Run an R expression in the background (only on UNIX)
#' 
#' @description Generate an R code of the expression that is copied via \code{scp}
#' to any machine (ssh-key needed). Then collect the results.
#' @details \code{runbg()} generates a workspace from the \code{input} argument
#' and copies the workspace and all C files or .so files to the remote machines via
#' \code{scp}. This will only work if *an ssh-key had been generated and added
#' to the authorized keys on the remote machine*. The
#' code snippet, i.e. the \code{...} argument, can include several intermediate results
#' but only the last call which is not redirected into a variable is returned via the
#' variable \code{.runbgOutput}, see example below.
#' @param ... Some R code
#' @param filename Character, defining the filename of the temporary file. Random
#' file name ist chosen if NULL.
#' @param machine Character vector, e.g. \code{"localhost"} or \code{"knecht1.fdm.uni-freiburg.de"}
#' or \code{c(localhost, localhost)}.
#' @param input Character vector, the objects in the workspace that are stored
#' into an R data file and copied to the remove machine.
#' @param compile Logical. If \code{TRUE}, C files are copied and compiled on the remote machine.
#' Otherwise, the .so files are copied.
#' @param wait Logical. Wait until executed. If \code{TRUE}, the code checks if the result file
#' is already present in which case it is loaded. If not present, \code{runbg()} starts, produces
#' the result and loads it as \code{.runbgOutput} directly into the workspace. If \code{wait = FALSE},
#' \code{runbg()} starts in the background and the result is only loaded into the workspace
#' when the \code{get()} function is called, see Value section. 
#' @param recover Logical. This option is useful to recover the three functions check(), get() and purge(), 
#' e.g. when a session has crashed. Then, the three functions are recreated without restarting the job.
#' They can then be used to get the results of a job wihtout having to do it manually. 
#' Requires the correct filename, so if the previous runbg was run with filename = NULL, you have to 
#' specify the tmp_filename manually.
#' @param walltime Optional. Maximum runtime in the format "HH:MM:SS". If exceeded, the job will be terminated.
#' @return List of functions \code{check()}, \code{get()}, \code{purge()} and \code{terminate()}. 
#' \code{check()} checks if the result is ready.
#' \code{get()} copies the result file
#' to the working directory and loads it into the workspace as an object called \code{.runbgOutput}. 
#' This object is a list named according to the machines that contains the results returned by each
#' machine.
#' \code{purge()} deletes the temporary folder
#' from the working directory and the remote machines.
#' \code{terminate()} kills all running processes associated with this job on the remote machines.
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
#' #' Recover a runbg job with the option "recover"
#' out_job1 <- runbg({
#'          M <- matrix(rnorm(1e2), 10, 10)
#'          solve(M)
#'          }, machine = c("localhost", "localhost"), filename = "job1")
#' Sys.sleep(1)
#' remove(out_job1)
#' try(out_job1$check())
#' out_job1 <- runbg({
#'   "This code is not run"
#' }, machine = c("localhost", "localhost"), filename = "job1", recover = T)
#' out_job1$get()
#' result <- .runbgOutput
#' print(result)
#' out_job1$purge()
#' }
runbg <- function(..., machine = "localhost", filename = NULL, input = ls(.GlobalEnv), compile = FALSE, wait = FALSE, recover = F, walltime = NULL) {
  
  
  expr <- as.expression(substitute(...))
  nmachines <- length(machine)
  
  # Set file name
  if (is.null(filename))
    filename <- paste0("tmp_", paste(sample(c(0:9, letters), 5, replace = TRUE), collapse = ""))
  
  filename0 <- filename
  filename <- paste(filename, 1:nmachines, sep = "_")
  
  # Initialize output
  out <- structure(vector("list", 4), names = c("check", "get", "purge", "terminate"))
  
  # Check
  out[[1]] <- function() {
    
    check.out <- sapply(1:nmachines, function(m) length(suppressWarnings(
      system(paste0("ssh ", machine[m], " ls ", filename[m], "_folder/ | grep -x ", filename[m], "_result.RData"), 
             intern = TRUE))))
    
    if (all(check.out) > 0) {
      cat("Result is ready!\n")
      return(TRUE)
    }
    else if (any(check.out) > 0) {
      cat("Result from machines", paste(which(check.out > 0), collapse = ", "), "are ready.")
      return(FALSE)
    }
    else if (all(check.out) == 0) {
      cat("Not ready!\n") 
      return(FALSE)
    }
      
  }
  
  # Get
  out[[2]] <- function() {
    
    result <- structure(vector(mode = "list", length = nmachines), names = machine)
    for (m in 1:nmachines) {
      .runbgOutput <- NULL
      system(paste0("scp ", machine[m], ":", filename[m], "_folder/", filename[m], "_result.RData ./"), ignore.stdout = TRUE, ignore.stderr = TRUE)
      check <- try(load(file = paste0(filename[m], "_result.RData")), silent = TRUE) 
      if (!inherits("try-error", check)) result[[m]] <- .runbgOutput
    }
    
    .GlobalEnv$.runbgOutput <- result
    
  }
  
  # Purge
  out[[3]] <- function() {
    
    for (m in 1:nmachines) {
      # Check if folder exists before trying to remove it
      folder_exists <- suppressWarnings(
        system(paste0("ssh ", machine[m], " '[ -d ", filename[m], "_folder ] && echo 1 || echo 0'"), 
               intern = TRUE)
      )
      if (folder_exists == "1") {
        system(paste0("ssh ", machine[m], " rm -r ", filename[m], "_folder"))
      }
      
      # Check if .Rout file exists before trying to remove it
      rout_exists <- suppressWarnings(
        system(paste0("ssh ", machine[m], " '[ -f ", filename[m], ".Rout ] && echo 1 || echo 0'"), 
               intern = TRUE)
      )
      if (rout_exists == "1") {
        system(paste0("ssh ", machine[m], " rm ", filename[m], ".Rout"))
      }
    }
    
    # Check if local files exist before trying to remove them
    local_files <- list.files(pattern = paste0(filename0, ".*"))
    if (length(local_files) > 0) {
      system(paste0("rm ", filename0, "*"))
    }
  }
  
  # Add terminate function
  out[[4]] <- function() {
    for (m in 1:nmachines) {
      # Get process IDs for the R processes running our job
      # Use ps to get more detailed process info including state
      pids <- suppressWarnings(
        system(paste0("ssh ", machine[m], 
               " 'ps aux | grep \"", filename[m], "\" | grep -v grep'"), 
               intern = TRUE)
      )
      
      if (length(pids) > 0) {
        # Extract PIDs of only running (R state) processes
        running_pids <- sapply(strsplit(pids, "\\s+"), function(x) {
          # Check if process state contains 'R' (running)
          if (grepl("R", x[8])) x[2] else NULL
        })
        running_pids <- running_pids[!sapply(running_pids, is.null)]
        
        if (length(running_pids) > 0) {
          # Kill only the running processes
          system(paste0("ssh ", machine[m], " 'kill ", paste(running_pids, collapse = " "), "'"))
          cat("Terminated", length(running_pids), "running processes on", machine[m], "\n")
        } else {
          cat("No running processes found on", machine[m], "\n")
        }
      }
    }
  }
  
  
  # Recover the three functions check, get, purge without letting the job being evaluated again.
  if (recover) return(out)
  
  
  
  # Check if filenames exist and load last result (only if wait == TRUE)
  resultfile <- paste(filename, "result.RData", sep = "_")
  if (all(file.exists(resultfile)) & wait) {
    
    for (m in 1:nmachines) {
      
      result <- structure(vector(mode = "list", length = nmachines), names = machine)
      load(file = resultfile[m])
      result[[m]] <- .runbgOutput
      
    }
    .GlobalEnv$.runbgOutput <- result
    return(out)
  }
  
  # Save current workspace
  save(list = input, file = paste0(filename0, ".RData"), envir = .GlobalEnv)
  
  # Get loaded packages
  pack <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  pack <- pack[!is.na(pack)]
  pack <- paste(paste0("try(library(", pack, "))"), collapse = "\n")
  
  # Define outputs
  output <- ".runbgOutput"
  
  compile.line <- NULL
  if (compile)
    compile.line <- paste(
      "cfiles <- list.files(pattern = '.c$')",
      "cppfiles <- list.files(pattern = '.cpp$')",
      "filelist <- paste(paste(cfiles, collapse = ' '), paste(cppfiles, collapse = ' '))",
      paste0("filename0 <- '", filename0, "'"),
      "system(paste0('R CMD SHLIB ', filelist, ' -o ', filename0, '.so'))",
      sep = "\n")
   
  # Write program into character
  program <- lapply(1:nmachines, function(m) paste(
    pack,
    paste0("setwd('~/", filename[m], "_folder')"),
    "rm(list = ls())",
    if (!is.null(walltime)) {
      paste0("Sys.setenv(R_TIMEOUT='", walltime, "')")
    },
    compile.line,
    paste0("load('", filename0, ".RData')"),
    "files <- list.files(pattern = '.so')",
    "for (f in files) dyn.load(f)",
    paste0(".node <- ", m),
    if (!is.null(walltime)) {
      paste0(".runbgOutput <- try(tools::pskill(Sys.getpid(), tools::SIGALRM, ", walltime, "); ", as.character(expr), ")")
    } else {
      paste0(".runbgOutput <- try(", as.character(expr), ")")
    },
    paste0("save(", output ,", file = '", filename[m], "_result.RData')"),
    sep = "\n"
  ))
  
  # Copy files to the temporary folder on the remote machine, write code and run
  for (m in 1:nmachines) {
    
    # Write program code into file
    cat(program[[m]], file = paste0(filename[m], ".R"))
    
    # Create a working directory on the remote machine
    system(paste0("ssh ", machine[m], " mkdir -p ", filename[m], "_folder/"), 
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Remove any old files in the working directory
    system(paste0("ssh ", machine[m], " rm -r ", filename[m], "_folder/*"), 
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Copy RData files to the working directory on the remote machine
    system(paste0("scp ", getwd(), "/", filename0, ".RData* ", machine[m], ":", filename[m], "_folder/"))
    
    # Copy the R scripts to the working directory on the remote machine
    system(paste0("scp ", getwd(), "/", filename[m], ".R* ", machine[m], ":", filename[m], "_folder/"))
    
    # Copy compiled files if `compile` is set to TRUE, otherwise copy shared object files
    if (compile) {
      system(paste0("scp ", getwd(), "/*.c ", getwd(), "/*.cpp ", machine[m], ":", filename[m], "_folder/"))
    } else {
      system(paste0("scp ", getwd(), "/*.so ", machine[m], ":", filename[m], "_folder/"))
    }
    
    # Execute the R script on the remote machine with the necessary environment variables
    # export OMP_NUM_THREADS=1 & MKL_NUM_THREADS=1 ensure that each job uses only one thread
    system(paste0(
      "ssh ", machine[m], 
      " 'export OMP_NUM_THREADS=1 && export MKL_NUM_THREADS=1 && ", 
      "R CMD BATCH --vanilla ", filename[m], "_folder/", filename[m], ".R'"
    ), intern = FALSE, wait = wait)
  }
  
  
  if (wait) {
    out$get()
    out$purge()
  } else {
    return(out)
  }
  
  
}


#' Run an R expression on the bwForCluster
#' 
#' @description Generate an R code of the expression that is copied via \code{scp}
#' to the bwForCluster (ssh-key needed). Then collect the results.
#' @details \code{runbg()} generates a workspace from the \code{input} argument
#' and copies the workspace and all C files or .so files to the remote machines via
#' \code{scp}. This will only work if *an ssh-key had been generated and added
#' to the authorized keys on the remote machine*. The
#' code snippet, i.e. the \code{...} argument, can include several intermediate results
#' but only the last call which is not redirected into a variable is returned via the
#' variable \code{.runbgOutput}, see example below.
#' @param ... Some R code
#' @param machine e.g. \code{fr_dk846@@bwfor.cluster.uni-mannheim.de}
#' @param filename Character, defining the filename of the temporary file. Random
#' file name ist chosen if NULL.
#' @param nodes Number of nodes, e.g. 10
#' @param cores Number of cores, e.g. 16
#' @param walltime estimated runtime in the format \code{hh:mm:ss}, e.g. \code{01:30:00}.
#' Jobs with a walltime up to 30 min are sent to a quick queue. When the walltime
#' is exceeded, all jobs are automatically killed by the queue.
#' @param input Character vector, the objects in the workspace that are stored
#' into an R data file and copied to the remove machine.
#' @param compile Logical. If \code{TRUE}, C files are copied and compiled on the remote machine.
#' Otherwise, the .so files are copied.
#' @param recover Logical, If \code{TRUE}, the scripts will not be started again.
#' Can be used to get back the check and get functions for an already started process, e.g. after local session has aborted.
#' @return List of functions \code{check()}, \code{get()} and \code{purge()}. 
#' \code{check()} checks, if the result is ready.
#' \code{get()} copies the result file
#' to the working directory and loads it into the workspace as an object called \code{.runbgOutput}. 
#' This object is a list named according to the machines that contains the results returned by each
#' machine.
#' \code{purge()} deletes the temporary folder
#' from the working directory and the remote machines.
#' @examples
#' \dontrun{
#' out_job1 <- runbg({
#'    mstrust(obj, center, fits = 10, cores = 2)
#'  }, 
#'  machine = "bwfor", nodes = 2, cores = "2:best", 
#'  walltime = "00:01:00", 
#'  filename = "job1")
#' out_job1$check()          
#' out_job1$get()
#' out_job1$purge()
#' result <- .runbgOutput
#' print(result)
#' }
#' 
#' @export
runbg_bwfor <- function(..., machine, filename = NULL, nodes = 1, cores = 1, walltime = "01:00:00", input = ls(.GlobalEnv), compile = TRUE, recover = F) {
  
  
  expr <- as.expression(substitute(...))
  
  # Set file name
  if (is.null(filename))
    filename <- paste0("tmp_", paste(sample(c(0:9, letters), 5, replace = TRUE), collapse = ""))
  
  filename0 <- filename
  filename <- paste(filename, 1:nodes, sep = "_")
  
  # Initialize output
  out <- structure(vector("list", 3), names = c("check", "get", "purge"))
  out[[1]] <- function() {
    
    check.out <- length(suppressWarnings(
      system(paste0("ssh ", machine, " ls ", filename0, "_folder/ | grep result.RData"), 
             intern = TRUE)))
    
    if (check.out == nodes) {
      cat("Result is ready!\n")
      return(TRUE)
    }
    else if (check.out  < nodes) {
      cat("Result from", check.out, "out of", nodes, "nodes are ready.")
      return(FALSE)
    }
    
    
  }
  
  out[[2]] <- function() {
    
    result <- structure(vector(mode = "list", length = nodes))
    system(paste0("scp ", machine, ":", filename0, "_folder/*", "_result.RData ./"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    for (m in 1:nodes) {
      .runbgOutput <- NULL
      check <- try(load(file = paste0(filename[m], "_result.RData")), silent = TRUE) 
      if (!inherits("try-error", check)) result[[m]] <- .runbgOutput
    }
    
    .GlobalEnv$.runbgOutput <- result
    
  }
  
  out[[3]] <- function() {
    
    system(paste0("rm ", filename0, "*"), wait = TRUE)
    system(paste0("ssh ", machine, " rm -r ", filename0, "*"), wait = TRUE)
    
  }
  if(recover) return(out)
  
  # Save current workspace
  save(list = input, file = paste0(filename0, ".RData"), envir = .GlobalEnv)
  
  # Get loaded packages
  pack <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  pack <- pack[!is.na(pack)]
  pack <- paste(paste0("try(library(", pack, "))"), collapse = "\n")
  
  # Define outputs
  output <- ".runbgOutput"

 

  # Write program into character
  program <- lapply(1:nodes, function(m) {
    paste(
      pack,
      paste0("setwd('~/", filename0, "_folder')"),
      "rm(list = ls())",
      "library(doParallel)",
      "procs <- as.numeric(Sys.getenv('MOAB_PROCCOUNT'))",
      "registerDoParallel(cores=procs)",
      paste0("load('", filename0, ".RData')"),
      "files <- list.files(pattern = '.so')",
      "for (f in files) dyn.load(f)",
      paste0(".node <- ", m),
      paste0(".runbgOutput <- try(", as.character(expr), ")"),
      
      paste0("save(", output ,", file = '", filename[m], "_result.RData')"),
      sep = "\n"
    )
  })
  
  # Write program code into file
  for (m in 1:nodes) cat(program[[m]], file = paste0(filename[m], ".R"))
  
  # Write job file to be called by msub
  job <- lapply(1:nodes, function(m) {
    paste(
      "#!/bin/sh", 
      "########## Begin MOAB/Slurm header ##########",
      "#",
      "# Give job a reasonable name",
      paste0("#MOAB -N ", filename[m]),
      "#",
      "# Request number of nodes and CPU cores per node for job",
      paste0("#MOAB -l nodes=1:ppn=", cores),
      "#",
      "# Estimated wallclock time for job",
      paste0("#MOAB -l walltime=", walltime),
      "#",
      "# Write standard output and errors in same file",
      "#MOAB -j oe ",
      "#",
      "########### End MOAB header ##########",
      "",
      "# Setup R Environment",
      "module load math/R",
      "export OPENBLAS_NUM_THREADS=1",
      "# Start program",
      paste0("R CMD BATCH --no-save --no-restore --slave ", filename0, "_folder/", filename[m], ".R"),
      sep = "\n"
    )
  })
  
  # Write job file to file
  for (m in 1:nodes) cat(job[[m]], file = paste0(filename[m], ".moab"))
  
  # Copy files to temporal folder
  system(paste0("ssh ", machine, " mkdir ", filename0, "_folder/"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(paste0("ssh ", machine, " rm ", filename0, "_folder/*"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(paste0("scp ", getwd(), "/", filename0, ".RData* ", machine, ":", filename0, "_folder/"))
  system(paste0("scp ", getwd(), "/", filename0, "*.R* ", machine, ":", filename0, "_folder/"))
  system(paste0("scp ", getwd(), "/", filename0, "*.moab ", machine, ":"))
  if (compile) {
    
    sourcefiles <- paste(
      paste0(
        filename0, "_folder/", 
        c(list.files(pattern = glob2rx("*.c")), list.files(pattern = glob2rx("*.cpp")))
      ), 
      collapse = " "
    )
    
    system(paste0("scp ", getwd(), "/*.c ", getwd(), "/*.cpp ", machine, ":", filename0, "_folder/"))
    system(paste0("ssh ", machine, " 'module load math/R; R CMD SHLIB ", sourcefiles, " -o ", filename0, "_folder/", filename0, ".so'"))
  } else {
    system(paste0("scp ", getwd(), "/*.so ", machine, ":", filename0, "_folder/"))
  }
  
  
  # Run in background
  for (m in 1:nodes) system(paste0("ssh ", machine, " msub ", filename[m], ".moab"), intern = FALSE)
  
  
  return(out)
  
  
}


#' Run an R expression on the bwForCluster via sshpass
#' 
#' @description Generate an R code of the expression that is copied via \code{scp}
#' to the bwForCluster. Then collect the results. ssh-key not needed. Password can be provided via an additional argument.
#' sshpass needs to be installed on your local machine.
#' @details \code{runbg()} generates a workspace from the \code{input} argument
#' and copies the workspace and all C files or .so files to the remote machines via
#' \code{scp}. This will only work if *an ssh-key had been generated and added
#' to the authorized keys on the remote machine*. The
#' code snippet, i.e. the \code{...} argument, can include several intermediate results
#' but only the last call which is not redirected into a variable is returned via the
#' variable \code{.runbgOutput}, see example below.
#' @param ... Some R code
#' @param machine e.g. \code{fr_dk846@@bwfor.cluster.uni-mannheim.de}
#' @param filename Character, defining the filename of the temporary file. Random
#' file name ist chosen if NULL.
#' @param nodes Number of nodes, e.g. 10
#' @param cores Number of cores, e.g. 16
#' @param walltime estimated runtime in the format \code{hh:mm:ss}, e.g. \code{01:30:00}.
#' Jobs with a walltime up to 30 min are sent to a quick queue. When the walltime
#' is exceeded, all jobs are automatically killed by the queue.
#' @param input Character vector, the objects in the workspace that are stored
#' into an R data file and copied to the remove machine.
#' @param compile Logical. If \code{TRUE}, C files are copied and compiled on the remote machine.
#' Otherwise, the .so files are copied.
#' @param recover Logical, If \code{TRUE}, the scripts will not be started again.
#' Can be used to get back the check and get functions for an already started process, e.g. after local session has aborted.
#' @param password Your ssh password in plain text (yes, no joke unfortunately), the password is handed over to sshpass for automatic login on the cluster.
#' @return List of functions \code{check()}, \code{get()} and \code{purge()}. 
#' \code{check()} checks, if the result is ready.
#' \code{get()} copies the result file
#' to the working directory and loads it into the workspace as an object called \code{.runbgOutput}. 
#' This object is a list named according to the machines that contains the results returned by each
#' machine.
#' \code{purge()} deletes the temporary folder
#' from the working directory and the remote machines.
#' @examples
#' \dontrun{
#' out_job1 <- runbg({
#'    mstrust(obj, center, fits = 10, cores = 2)
#'  }, 
#'  machine = "bwfor", nodes = 2, cores = "2:best", 
#'  walltime = "00:01:00", 
#'  filename = "job1")
#' out_job1$check()          
#' out_job1$get()
#' out_job1$purge()
#' result <- .runbgOutput
#' print(result)
#' }
#' 
#' @export
runbg_bwfor_sshpass <- function(..., machine, filename = NULL, nodes = 1, cores = 1, walltime = "01:00:00", input = ls(.GlobalEnv), compile = TRUE, recover = F, password="'begin__end'") {
  
  
  expr <- as.expression(substitute(...))
  
  # Set file name
  if (is.null(filename))
    filename <- paste0("tmp_", paste(sample(c(0:9, letters), 5, replace = TRUE), collapse = ""))
  
  filename0 <- filename
  filename <- paste(filename, 1:nodes, sep = "_")
  
  # Initialize output
  out <- structure(vector("list", 3), names = c("check", "get", "purge"))
  out[[1]] <- function() {
    
    check.out <- length(suppressWarnings(
      system(paste0("sshpass -p ",password, " ssh ", machine, " ls ", filename0, "_folder/ | grep result.RData"), 
             intern = TRUE)))
    
    if (check.out == nodes) {
      cat("Result is ready!\n")
      return(TRUE)
    }
    else if (check.out  < nodes) {
      cat("Result from", check.out, "out of", nodes, "nodes are ready.")
      return(FALSE)
    }
    
    
  }
  
  out[[2]] <- function() {
    
    result <- structure(vector(mode = "list", length = nodes))
    system(paste0("sshpass -p ",password, " scp ", machine, ":", filename0, "_folder/*", "_result.RData ./"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    for (m in 1:nodes) {
      .runbgOutput <- NULL
      check <- try(load(file = paste0(filename[m], "_result.RData")), silent = TRUE) 
      if (!inherits("try-error", check)) result[[m]] <- .runbgOutput
    }
    
    .GlobalEnv$.runbgOutput <- result
    
  }
  
  out[[3]] <- function() {
    
    system(paste0("rm ", filename0, "*"), wait = TRUE)
    system(paste0("sshpass -p ",password, " ssh ", machine, " rm -r ", filename0, "*"), wait = TRUE)
    
  }
  if(recover) return(out)
  
  # Save current workspace
  save(list = input, file = paste0(filename0, ".RData"), envir = .GlobalEnv)
  
  # Get loaded packages
  pack <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  pack <- pack[!is.na(pack)]
  pack <- paste(paste0("try(library(", pack, "))"), collapse = "\n")
  
  # Define outputs
  output <- ".runbgOutput"
  
  
  
  # Write program into character
  program <- lapply(1:nodes, function(m) {
    paste(
      pack,
      paste0("setwd('~/", filename0, "_folder')"),
      "rm(list = ls())",
      "library(doParallel)",
      "procs <- as.numeric(Sys.getenv('MOAB_PROCCOUNT'))",
      "registerDoParallel(cores=procs)",
      paste0("load('", filename0, ".RData')"),
      "files <- list.files(pattern = '.so')",
      "for (f in files) dyn.load(f)",
      paste0(".node <- ", m),
      paste0(".runbgOutput <- try(", as.character(expr), ")"),
      
      paste0("save(", output ,", file = '", filename[m], "_result.RData')"),
      sep = "\n"
    )
  })
  
  # Write program code into file
  for (m in 1:nodes) cat(program[[m]], file = paste0(filename[m], ".R"))
  
  # Write job file to be called by msub
  job <- lapply(1:nodes, function(m) {
    paste(
      "#!/bin/sh", 
      "########## Begin MOAB/Slurm header ##########",
      "#",
      "# Give job a reasonable name",
      paste0("#MOAB -N ", filename[m]),
      "#",
      "# Request number of nodes and CPU cores per node for job",
      paste0("#MOAB -l nodes=1:ppn=", cores),
      "#",
      "# Estimated wallclock time for job",
      paste0("#MOAB -l walltime=", walltime),
      "#",
      "# Write standard output and errors in same file",
      "#MOAB -j oe ",
      "#",
      "########### End MOAB header ##########",
      "",
      "# Setup R Environment",
      "module load math/R",
      "export OPENBLAS_NUM_THREADS=1",
      "# Start program",
      paste0("R CMD BATCH --no-save --no-restore --slave ", filename0, "_folder/", filename[m], ".R"),
      sep = "\n"
    )
  })
  
  # Write job file to file
  for (m in 1:nodes) cat(job[[m]], file = paste0(filename[m], ".moab"))
  
  # Copy files to temporal folder
  system(paste0("sshpass -p ",password, " ssh ", machine, " mkdir ", filename0, "_folder/"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(paste0("sshpass -p ",password, " ssh ", machine, " rm ", filename0, "_folder/*"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(paste0("sshpass -p ",password, " scp ", getwd(), "/", filename0, ".RData* ", machine, ":", filename0, "_folder/"))
  system(paste0("sshpass -p ",password, " scp ", getwd(), "/", filename0, "*.R* ", machine, ":", filename0, "_folder/"))
  system(paste0("sshpass -p ",password, " scp ", getwd(), "/", filename0, "*.moab ", machine, ":"))
  if (compile) {
    
    sourcefiles <- paste(
      paste0(
        filename0, "_folder/", 
        c(list.files(pattern = glob2rx("*.c")), list.files(pattern = glob2rx("*.cpp")))
      ), 
      collapse = " "
    )
    
    system(paste0("sshpass -p ",password, " scp ", getwd(), "/*.c ", getwd(), "/*.cpp ", machine, ":", filename0, "_folder/"))
    system(paste0("sshpass -p ",password, " ssh ", machine, " 'module load math/R; R CMD SHLIB ", sourcefiles, " -o ", filename0, "_folder/", filename0, ".so'"))
  } else {
    system(paste0("sshpass -p ",password, " scp ", getwd(), "/*.so ", machine, ":", filename0, "_folder/"))
  }
  
  
  # Run in background
  for (m in 1:nodes) system(paste0("sshpass -p ",password, " ssh ", machine, " msub ", filename[m], ".moab"), intern = FALSE)
  
  
  return(out)
  
  
}



#' Run an R expression on the bwForCluster via sshpass and slurm
#' 
#' @description Generate an R code of the expression that is copied via \code{scp}
#' to the bwForCluster. Then collect the results. ssh-key not needed. Password can be provided via an additional argument.
#' sshpass needs to be installed on your local machine.
#' @details \code{runbg()} generates a workspace from the \code{input} argument
#' and copies the workspace and all C files or .so files to the remote machines via
#' \code{scp}. This will only work if *an ssh-key had been generated and added
#' to the authorized keys on the remote machine*. The
#' code snippet, i.e. the \code{...} argument, can include several intermediate results
#' but only the last call which is not redirected into a variable is returned via the
#' variable \code{.runbgOutput}, see example below.
#' @param ... Some R code
#' @param machine e.g. \code{fr_dk846@@bwfor.cluster.uni-mannheim.de}
#' @param filename Character, defining the filename of the temporary file. Random
#' file name ist chosen if NULL. Must not contain the string "Minus".
#' @param nodes Number of nodes, e.g. 10
#' @param cores Number of cores, e.g. 16
#' @param partition character, the partition where to start the job
#' @param walltime estimated runtime in the format \code{hh:mm:ss}, e.g. \code{01:30:00}.
#' Jobs with a walltime up to 30 min are sent to a quick queue. When the walltime
#' is exceeded, all jobs are automatically killed by the queue.
#' @param input Character vector, the objects in the workspace that are stored
#' into an R data file and copied to the remove machine.
#' @param compile Logical. If \code{TRUE}, C files are copied and compiled on the remote machine.
#' Otherwise, the .so files are copied.
#' @param recover Logical, If \code{TRUE}, the scripts will not be started again.
#' Can be used to get back the check and get functions for an already started process, e.g. after local session has aborted.
#' @param password Your ssh password in plain text (yes, no joke unfortunately), the password is handed over to sshpass for automatic login on the cluster. If NULL, the standard ssh/scp is used and you will be asked for your password multiple times while uploading the scripts.
#' @return List of functions \code{check()}, \code{get()} and \code{purge()}. 
#' \code{check()} checks, if the result is ready.
#' \code{get()} copies the result file
#' to the working directory and loads it into the workspace as an object called \code{.runbgOutput}. 
#' This object is a list named according to the machines that contains the results returned by each
#' machine.
#' \code{purge()} deletes the temporary folder
#' from the working directory and the remote machines.
#' @examples
#' \dontrun{
#' out_job1 <- runbg({
#'    mstrust(obj, center, fits = 10, cores = 2)
#'  }, 
#'  machine = "bwfor", nodes = 2, cores = "2:best", 
#'  walltime = "00:01:00", 
#'  filename = "job1")
#' out_job1$check()          
#' out_job1$get()
#' out_job1$purge()
#' result <- .runbgOutput
#' print(result)
#' }
#' 
#' @export
runbg_bwfor_slurm <- function(..., machine, filename = NULL, nodes = 1, cores = 1, partition = "single", walltime = "01:00:00", input = ls(.GlobalEnv), compile = TRUE, recover = F, password="'begin__end'") {
  
  if(is.null(password)){
    ssh_command <- "ssh "
    scp_command <- "scp "
  } else {
    ssh_command <- paste0("sshpass -p ", password, " ssh ")
    scp_command <- paste0("sshpass -p ", password, " scp ")
  }
  
  
  
  expr <- as.expression(substitute(...))
  
  # Set file name
  if (is.null(filename))
    filename <- paste0("tmp_", paste(sample(c(0:9, letters), 5, replace = TRUE), collapse = ""))
  
  filename0 <- filename
  filename <- paste(filename, 1:nodes, sep = "_")
  
  # Initialize output
  out <- structure(vector("list", 3), names = c("check", "get", "purge"))
  out[[1]] <- function() {
    
    check.out <- length(suppressWarnings(
      system(paste0(ssh_command, machine, " ls ", filename0, "_folder/ | grep result.RData"), 
             intern = TRUE)))
    
    if (check.out == nodes) {
      cat("Result is ready!\n")
      return(TRUE)
    }
    else if (check.out  < nodes) {
      cat("Result from", check.out, "out of", nodes, "nodes are ready.")
      return(FALSE)
    }
    
    
  }
  
  out[[2]] <- function() {
    
    result <- structure(vector(mode = "list", length = nodes))
    system(paste0(scp_command, machine, ":", filename0, "_folder/*", "_result.RData ./"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    for (m in 1:nodes) {
      .runbgOutput <- NULL
      check <- try(load(file = paste0(filename[m], "_result.RData")), silent = TRUE) 
      if (!inherits("try-error", check)) result[[m]] <- .runbgOutput
    }
    
    .GlobalEnv$.runbgOutput <- result
    
  }
  
  out[[3]] <- function() {
    
    system(paste0("rm ", filename0, "*"), wait = TRUE)
    system(paste0(ssh_command, machine, " rm -r ", filename0, "*"), wait = TRUE)
    
  }
  if(recover) return(out)
  
  # Save current workspace
  save(list = input, file = paste0(filename0, ".RData"), envir = .GlobalEnv)
  
  # Get loaded packages
  pack <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  pack <- pack[!is.na(pack)]
  pack <- paste(paste0("try(library(", pack, "))"), collapse = "\n")
  
  # Define outputs
  output <- ".runbgOutput"
  
  
  
  # Write program into character
  program <- lapply(1:nodes, function(m) {
    paste(
      pack,
      paste0("setwd('~/", filename0, "_folder')"),
      "rm(list = ls())",
      "library(doParallel)",
      "procs <- as.numeric(Sys.getenv('SLURM_NTASKS'))",
      "registerDoParallel(cores=procs)",
      paste0("load('", filename0, ".RData')"),
      "files <- list.files(pattern = '.so')",
      "for (f in files) dyn.load(f)",
      paste0(".node <- ", m),
      paste0(".runbgOutput <- try(", as.character(expr), ")"),
      
      paste0("save(", output ,", file = '", filename[m], "_result.RData')"),
      sep = "\n"
    )
  })
  
  # Write program code into file
  for (m in 1:nodes) cat(program[[m]], file = paste0(filename[m], ".R"))
  
  # Write job file to be called by msub
  job <- lapply(1:nodes, function(m) {
    paste(
      "#!/bin/sh", 
      "########## Begin MOAB/Slurm header ##########",
      "#",
      "# Give job a reasonable name",
      paste0("#SBATCH -J ", filename[m]),
      "#",
      "# Request number of nodes and CPU cores per node for job",
      paste0("#SBATCH -n ", cores),
      "#",
      "# Estimated wallclock time for job",
      paste0("#SBATCH -t ", walltime),
      "#",
      "# Request correct partition",
      paste0("#SBATCH --partition ", partition),
      "#",
      "########### End MOAB header ##########",
      "",
      "# Setup R Environment",
      "module load math/R",
      "export OPENBLAS_NUM_THREADS=1",
      "# Start program",
      paste0("R CMD BATCH --no-save --no-restore --slave ", filename0, "_folder/", filename[m], ".R"),
      sep = "\n"
    )
  })
  
  # Write job file to file
  for (m in 1:nodes) cat(job[[m]], file = paste0(filename[m], ".sh"))
  
  # Copy files to temporal folder
  system(paste0(ssh_command, machine, " mkdir ", filename0, "_folder/"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(paste0(ssh_command, machine, " rm ", filename0, "_folder/*"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(paste0(scp_command, getwd(), "/", filename0, ".RData* ", machine, ":", filename0, "_folder/"))
  system(paste0(scp_command, getwd(), "/", filename0, "*.R* ", machine, ":", filename0, "_folder/"))
  system(paste0(scp_command, getwd(), "/", filename0, "*.sh ", machine, ":"))
  if (compile) {
    
    sourcefiles <- paste(
      paste0(
        filename0, "_folder/", 
        c(list.files(pattern = glob2rx("*.c")), list.files(pattern = glob2rx("*.cpp")))
      ), 
      collapse = " "
    )
    
    system(paste0(scp_command, getwd(), "/*.c ", getwd(), "/*.cpp ", machine, ":", filename0, "_folder/"))
    system(paste0(ssh_command, machine, " 'module load math/R; R CMD SHLIB ", sourcefiles, " -o ", filename0, "_folder/", filename0, ".so'"))
  } else {
    system(paste0(scp_command, getwd(), "/*.so ", machine, ":", filename0, "_folder/"))
  }
  
  
  # Run in background
  for (m in 1:nodes) system(paste0(ssh_command, machine, " sbatch ", filename[m], ".sh"), intern = FALSE)
  
  
  return(out)
  
  
}



#' Remote install dMod to a ssh-reachable host
#' 
#' @description Install your local dMod version to a remote host via ssh.
#' @param sshtarget The ssh host url.
#' @param source If type = local, source must point to the source directory of
#'   your dMod version. This is most probably you local dMod git repository.
#' @param type Which dMod to install. At the moment, only your local version is
#'   supported.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @importFrom utils packageVersion
#' @export
runbgInstall <- function(sshtarget, source = NULL, type = "local") {
  
  if (type == "local") {
    # Build dMod package
    if (is.null(source)) {
      stop("dMod source location not specified.")
    }
    cat("* Preparing local dMod version for remote installation:\n")
    system(eval(paste("R CMD build --no-build-vignettes", source)))
    
    # Figure out package name
    dModPkg <- paste0("dMod_", packageVersion("dMod"), ".tar.gz")
    
    # Install to remote host
    cat(paste("* Installing to remote host", sshtarget, ":\n"))
    system(eval(paste0("scp ", dModPkg, " ", sshtarget, ":~/")))
    system(eval(paste("ssh", sshtarget, "R CMD INSTALL", dModPkg)))
    
    unlink(dModPkg)
  }
  
}


