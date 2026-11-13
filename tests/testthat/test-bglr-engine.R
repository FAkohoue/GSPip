make_test_grm <- function(ids) {
  n <- length(ids)
  G <- matrix(0.15, n, n, dimnames = list(ids, ids))
  diag(G) <- 1
  G
}

write_test_grm <- function(G, path) {
  utils::write.csv(
    data.frame(GenoID = rownames(G), G, check.names = FALSE),
    path, row.names = FALSE
  )
}

test_that("weighted BGLR RKHS effects are returned to the response scale", {
  converted <- GSPip:::.gspip_bglr_component(
    list(u = c(2, 6), varU = 1.25),
    label = "test",
    bglr_weights = c(2, 3)
  )
  expect_equal(converted$u, c(1, 2))
  expect_equal(converted$varU, 1.25)
})

test_that("BGLR controls reject undersized posterior samples", {
  skip_if_not_installed("BGLR")
  expect_error(
    GSPip::run_gblup_single(
      grm_file = "missing.csv",
      pheno_file = "missing.csv",
      output_dir = tempdir(),
      engine = "bglr",
      strata_col = NULL,
      weight_col = NULL,
      fix_stage1_error = FALSE,
      nonmarker_effects = FALSE,
      bglr_n_iter = 100L,
      bglr_burn_in = 50L,
      bglr_thin = 2L
    ),
    "at least 100 posterior samples"
  )
})

test_that("single-environment BGLR backend completes a small CV run", {
  skip_if_not_installed("BGLR")
  ids <- sprintf("G%02d", 1:12)
  G <- make_test_grm(ids)
  grm_path <- tempfile(fileext = ".csv")
  pheno_path <- tempfile(fileext = ".csv")
  out <- tempfile("gspip-bglr-single-")
  dir.create(out)
  on.exit(unlink(c(grm_path, pheno_path, out), recursive = TRUE, force = TRUE),
          add = TRUE)
  write_test_grm(G, grm_path)
  utils::write.csv(
    data.frame(
      GenoID = ids,
      Trait = seq_along(ids) + rep(c(-0.2, 0.2), 6),
      stringsAsFactors = FALSE
    ),
    pheno_path, row.names = FALSE
  )

  result <- run_gblup_single(
    grm_file = grm_path,
    pheno_file = pheno_path,
    output_dir = out,
    engine = "bglr",
    run_label = "small",
    trait_col = "Trait",
    strata_col = NULL,
    weight_col = NULL,
    fix_stage1_error = FALSE,
    nonmarker_effects = FALSE,
    k_folds = 2L,
    n_iterations = 1L,
    within_stratum_pa = FALSE,
    bglr_n_iter = 300L,
    bglr_burn_in = 100L,
    bglr_thin = 2L,
    bglr_seed = 17L
  )

  expect_s3_class(result, "gspip_gblup_single")
  expect_identical(result$engine, "bglr")
  expect_true(all(is.finite(result$cv_predictions[, 1])))
  expect_true("Vg_posterior" %in% result$metrics$Metric)
  expect_true(file.exists(file.path(out, "small_metrics.csv")))
})

test_that("G x E BGLR backend completes a small CV1 run", {
  skip_if_not_installed("BGLR")
  ids <- sprintf("G%02d", 1:10)
  G <- make_test_grm(ids)
  grm_path <- tempfile(fileext = ".csv")
  pheno_path <- tempfile(fileext = ".csv")
  out <- tempfile("gspip-bglr-gxe-")
  dir.create(out)
  on.exit(unlink(c(grm_path, pheno_path, out), recursive = TRUE, force = TRUE),
          add = TRUE)
  write_test_grm(G, grm_path)
  pheno <- expand.grid(
    GenoID = ids, Env = c("E1", "E2"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  pheno$Trait <- match(pheno$GenoID, ids) +
    ifelse(pheno$Env == "E2", 1.5, 0) +
    rep(c(-0.15, 0.15), length.out = nrow(pheno))
  pheno$Weight <- 4
  utils::write.csv(pheno, pheno_path, row.names = FALSE)

  result <- run_gblup_gxe(
    grm_file = grm_path,
    pheno_file = pheno_path,
    output_dir = out,
    engine = "bglr",
    run_label = "small_gxe",
    trait_col = "Trait",
    env_col = "Env",
    strata_col = NULL,
    weight_col = "Weight",
    fix_stage1_error = FALSE,
    nonmarker_effects = FALSE,
    cv_scheme = "CV1",
    k_folds = 2L,
    n_iterations = 1L,
    within_stratum_pa = FALSE,
    bglr_n_iter = 300L,
    bglr_burn_in = 100L,
    bglr_thin = 2L,
    bglr_seed = 19L
  )

  expect_s3_class(result, "gspip_gblup_gxe")
  expect_identical(result$engine, "bglr")
  expect_true(all(is.finite(result$predictions_by_environment$Mean)))
  expect_true("Vge_posterior" %in% result$metrics$Metric)
  expect_true(file.exists(file.path(out, "small_gxe_metrics.csv")))
})
