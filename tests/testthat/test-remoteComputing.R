test_that(".saveWorkspace round-trips through zstd", {

  skip_if(!nzchar(Sys.which("zstd")), "zstd not available")

  dir <- withr::local_tempdir()
  env <- new.env()
  env$mat <- matrix(seq_len(1e4), 100)
  env$chr <- rep(letters, 100)
  env$fun <- local({ z <- 42; function() z })

  target <- file.path(dir, "ws.RData")
  out    <- dMod:::.saveWorkspace(c("mat", "chr", "fun"), target, envir = env)

  expect_identical(out, paste0(target, ".zst"))
  expect_true(file.exists(out))
  expect_false(file.exists(target))

  ## The compressed workspace is what the transfer ships; a repeated structure
  ## like one parfn per condition is exactly what makes this worth doing.
  save(list = ls(env), file = paste0(target, ".raw"), envir = env, compress = FALSE)
  expect_lt(file.size(out), file.size(paste0(target, ".raw")))

  ## Same expansion the remote command chain performs before the job starts.
  expect_equal(system2("zstd", c("-dqf", shQuote(out))), 0L)
  back <- new.env()
  load(target, envir = back)
  expect_equal(back$mat, env$mat)
  expect_equal(back$chr, env$chr)
  expect_equal(back$fun(), 42)
})

test_that(".saveWorkspace falls back to an uncompressed write without zstd", {

  dir <- withr::local_tempdir()
  env <- new.env()
  env$val <- 1:10

  ## Sys.which() finds nothing on an empty PATH.
  withr::local_envvar(PATH = "")
  target <- file.path(dir, "ws.RData")
  out    <- dMod:::.saveWorkspace("val", target, envir = env)

  expect_identical(out, target)
  back <- new.env()
  load(target, envir = back)
  expect_equal(back$val, 1:10)
})

test_that(".remoteBuildScript bundles sources only where it can", {

  dir <- withr::local_tempdir()
  src <- file.path(dir, sprintf("m%02d.cpp", 1:20))
  for (f in src) writeLines("int f(void){return 1;}", f)

  ## Below the argument limit SHLIB gets the sources themselves, so bundles
  ## would only be compiled twice.
  expect_no_match(dMod:::.remoteBuildScript(files = src, output = "o.so", cores = 1,
                                            filelist = "f.txt"),
                  "dmod_bundle_", fixed = TRUE)

  testthat::local_mocked_bindings(.compileCmdLimit = function() 10L, .package = "dMod")
  chunked <- dMod:::.remoteBuildScript(files = src, output = "o.so", cores = 1,
                                       filelist = "f.txt", bundle = 5)
  expect_match(chunked, "BUNDLE=5")
  expect_match(chunked, "dmod_bundle_", fixed = TRUE)
  ## The archive is built from the work list, which now holds the bundles.
  expect_match(chunked, "sed 's/\\.[^.]*$/.o/' \"$WORK\"", fixed = TRUE)

  expect_no_match(dMod:::.remoteBuildScript(files = src, output = "o.so", cores = 1,
                                            filelist = "f.txt", bundle = 1),
                  "dmod_bundle_", fixed = TRUE)

  ## Object files are already compiled; there is nothing to bundle.
  obj <- sub("\\.cpp$", ".o", src)
  for (f in obj) writeLines("x", f)
  expect_no_match(dMod:::.remoteBuildScript(files = obj, output = "o.so", cores = 1,
                                            filelist = "f.txt",
                                            link = TRUE, bundle = 5),
                  "dmod_bundle_", fixed = TRUE)
})
