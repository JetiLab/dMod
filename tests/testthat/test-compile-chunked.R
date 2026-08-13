## Compiling a model with thousands of conditions: sources are always compiled
## one file per command, but naming every object on the final SHLIB command
## line overruns the OS per-argument limit. compile() then routes the link
## through a chunked static archive instead -- these tests pin that route down.

test_that(".compileChunks caps both count and command length", {

  files <- sprintf("/a/b/file%04d.o", 1:250)

  chunks <- .compileChunks(files, maxN = 100L)
  expect_equal(lengths(chunks), c(100L, 100L, 50L))
  expect_identical(unlist(chunks, use.names = FALSE), files)

  ## A single over-long name still forms a chunk of its own rather than
  ## producing an empty one.
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

  ## The cached key vector has to stay in step with the list it annotates,
  ## otherwise a later merge would dedup against the wrong entries.
  expect_equal(attr(m, "srckeys"), c("x.cpp", "y.cpp", "z.cpp"))
  expect_equal(.compileInfoKeys(m), attr(m, "srckeys"))
  expect_equal(.mergeCompileInfo(m, m), m)

  expect_null(.mergeCompileInfo(NULL, NULL))
  expect_null(.mergeCompileInfo(list(list(srcfile = character(0))), NULL))

})


test_that("compile() links many conditions through a chunked archive", {

  withr::local_dir(tempdir())
  ## Force the archive route with a handful of conditions instead of the
  ## thousand-odd it would take to overrun a real command line.
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
