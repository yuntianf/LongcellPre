test_that("native routines are registered in LongcellPre without a backend import", {
  ns <- asNamespace("LongcellPre")
  expect_false("Longcellsrc" %in% names(getNamespaceImports(ns)))
  calls <- getDLLRegisteredRoutines(getLoadedDLLs()[["LongcellPre"]])$.Call
  expect_length(calls, 22L)
  expect_true(all(startsWith(names(calls), "_LongcellPre_")))
  expect_true(all(vapply(calls, function(x) x$dll[["name"]] == "LongcellPre", logical(1))))
  expect_equal(isoform2sites("1,10|21,30", "|", ","), c(1, 10, 21, 30))
})

test_that("separate R workers load the bundled native routines", {
  previous <- future::plan()
  on.exit(future::plan(previous), add = TRUE)
  future::plan(future::multisession, workers = 2)
  results <- future.apply::future_lapply(1:2, function(i) {
    sites <- getFromNamespace("isoform2sites", "LongcellPre")
    list(pid = Sys.getpid(), sites = sites("1,10|21,30", "|", ","),
         dll = getLoadedDLLs()[["LongcellPre"]][["name"]],
         imports = names(getNamespaceImports("LongcellPre")))
  }, future.packages = "LongcellPre", future.globals = FALSE, future.seed = TRUE)
  expect_true(all(vapply(results, function(x) x$pid != Sys.getpid(), logical(1))))
  expect_true(all(vapply(results, function(x) identical(x$sites, c(1L, 10L, 21L, 30L)), logical(1))))
  expect_true(all(vapply(results, function(x) identical(x$dll, "LongcellPre"), logical(1))))
  expect_false(any(vapply(results, function(x) "Longcellsrc" %in% x$imports, logical(1))))
})
