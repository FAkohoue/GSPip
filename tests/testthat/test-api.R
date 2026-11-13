test_that("the complete programmatic API is exported", {
  expected <- c(
    "build_grm", "gspip_script", "optimize_training_population",
    "run_gblup_single", "run_gblup_gxe"
  )
  expect_true(all(expected %in% getNamespaceExports("GSPip")))
  expect_true(all(vapply(expected, function(x) is.function(getExportedValue("GSPip", x)),
                         logical(1))))
})

test_that("analysis functions expose required file arguments", {
  expect_true(all(c("grm_file", "output_dir") %in%
                    names(formals(optimize_training_population))))
  expect_true(all(c("grm_file", "pheno_file", "output_dir") %in%
                    names(formals(run_gblup_single))))
  expect_true(all(c("grm_file", "pheno_file", "output_dir") %in%
                    names(formals(run_gblup_gxe))))
})

test_that("loading the package does not execute workflow templates", {
  expect_false(exists("grm_file", envir = asNamespace("GSPip"), inherits = FALSE))
  expect_false(exists("pheno_file", envir = asNamespace("GSPip"), inherits = FALSE))
  expect_false(exists("output_dir", envir = asNamespace("GSPip"), inherits = FALSE))
})
