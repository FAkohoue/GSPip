test_that("build_grm passes NA calls to AGHmatrix when imputation is none", {
  path <- tempfile(fileext = ".csv")
  out <- tempfile("gspip-grm-")
  dir.create(out)
  on.exit(unlink(c(path, out), recursive = TRUE, force = TRUE), add = TRUE)
  utils::write.csv(data.frame(
    SNP = paste0("m", 1:6), CHR = c(1, 1, 1, 2, 2, 2), POS = seq(10, 60, 10),
    REF = c("A", "C", "G", "A", "C", "G"),
    ALT = c("G", "T", "A", "T", "G", "C"),
    S1 = c(0, 0, 2, 1, 0, 2),
    S2 = c(0, 1, 2, 0, NA, 1),
    S3 = c(2, 1, 0, 2, 2, 0),
    S4 = c(2, 2, 0, 1, 1, 1),
    check.names = FALSE
  ), path, row.names = FALSE, na = "NA")

  result <- build_grm(
    geno_file = path,
    output_dir = out,
    run_label = "none_test",
    qc_call_rate = 0,
    imputation_method = "none"
  )

  expect_s3_class(result, "gspip_grm")
  expect_equal(sum(is.na(result$markers)), 1)
  expect_identical(dim(result$grm), c(4L, 4L))
  expect_true(file.exists(file.path(out, "none_test_grm.csv")))
})

