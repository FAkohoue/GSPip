test_that("packaged example objects are complete and aligned", {
  e <- new.env(parent = emptyenv())
  utils::data("gs_example_all", package = "GSPip", envir = e)
  expect_true(exists("gs_example_all", envir = e, inherits = FALSE))
  x <- e$gs_example_all

  expect_identical(
    names(x),
    c("geno", "marker_map", "marker_weights", "blues_single", "blues_gxe",
      "candidate_ids", "target_ids", "grm", "truth")
  )
  expect_identical(dim(x$geno), c(120L, 240L))
  expect_identical(dim(x$marker_map), c(240L, 6L))
  expect_identical(dim(x$blues_single), c(120L, 5L))
  expect_identical(dim(x$blues_gxe), c(360L, 6L))
  expect_identical(dim(x$grm), c(120L, 120L))
  expect_identical(dim(x$truth), c(240L, 6L))

  expect_identical(rownames(x$geno), x$blues_single$GenoID)
  expect_identical(colnames(x$geno), x$marker_map$SNP)
  expect_identical(rownames(x$grm), rownames(x$geno))
  expect_identical(colnames(x$grm), rownames(x$geno))
  expect_equal(x$grm, t(x$grm), tolerance = 1e-12)

  expect_length(intersect(x$candidate_ids$GenoID, x$target_ids$GenoID), 0L)
  expect_setequal(c(x$candidate_ids$GenoID, x$target_ids$GenoID), rownames(x$geno))

  valid_precision <- function(d) {
    is.na(d$NBL) |
      (!is.na(d$Weight) & is.finite(d$Weight) & d$Weight > 0) |
      (!is.na(d$SE) & is.finite(d$SE) & d$SE > 0)
  }
  expect_true(all(valid_precision(x$blues_single)))
  expect_true(all(valid_precision(x$blues_gxe)))
  expect_false(anyDuplicated(x$blues_single$GenoID) > 0L)
  expect_false(anyDuplicated(paste(x$blues_gxe$GenoID, x$blues_gxe$Env,
                                   sep = "\r")) > 0L)
  expect_true(all(tapply(x$blues_gxe$Population, x$blues_gxe$GenoID,
                         function(z) length(unique(z))) == 1L))
})

test_that("all packaged flat files and manifest are present", {
  extdata <- system.file("extdata", package = "GSPip")
  expect_true(nzchar(extdata))
  manifest_path <- file.path(extdata, "example_manifest.csv")
  expect_true(file.exists(manifest_path))
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)

  expect_identical(names(manifest), c("File", "Rows", "Purpose"))
  expect_equal(nrow(manifest), 20L)
  expect_false(anyDuplicated(manifest$File) > 0L)

  flat_files <- basename(manifest$File[grepl("inst[/\\\\]extdata", manifest$File)])
  expect_length(flat_files, 10L)
  expect_true(all(file.exists(file.path(extdata, flat_files))))
})

test_that("numeric, HapMap, and VCF examples decode identically", {
  extdata <- system.file("extdata", package = "GSPip")
  numeric_result <- GSPip:::read_geno(
    file.path(extdata, "example_genotypes_numeric.csv"), format = "numeric",
    meta_cols = c("CHR", "POS", "REF", "ALT")
  )
  hapmap_result <- GSPip:::read_geno(
    file.path(extdata, "example_genotypes.hmp.txt"), format = "hapmap"
  )
  vcf_result <- GSPip:::read_geno(
    file.path(extdata, "example_genotypes.vcf"), format = "vcf"
  )

  expect_identical(dim(numeric_result$M), c(240L, 120L))
  expect_equal(numeric_result$M, hapmap_result$M)
  expect_equal(numeric_result$M, vcf_result$M)

  e <- new.env(parent = emptyenv())
  utils::data("gs_example_geno", package = "GSPip", envir = e)
  expect_equal(t(numeric_result$M), e$gs_example_geno)
})

test_that("each documented package dataset can be loaded", {
  objects <- paste0("gs_example_", c(
    "geno", "marker_map", "marker_weights", "blues_single", "blues_gxe",
    "candidate_ids", "target_ids", "grm", "truth", "all"
  ))
  for (object_name in objects) {
    e <- new.env(parent = emptyenv())
    utils::data(list = object_name, package = "GSPip", envir = e)
    expect_true(exists(object_name, envir = e, inherits = FALSE),
                info = object_name)
  }
})
