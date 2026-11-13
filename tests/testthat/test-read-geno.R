test_that("numeric input preserves missing calls and marker weights", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(data.frame(
    SNP = c("m1", "m2", "m3"), CHR = c(1, 1, 2), POS = c(10, 20, 30),
    REF = c("A", "C", "G"), ALT = c("G", "T", "A"),
    MarkerWeight = c(1, 0.5, 2),
    S1 = c(0, 1, NA), S2 = c(1, 2, 0), S3 = c(2, 0, 1),
    check.names = FALSE
  ), path, row.names = FALSE, na = "NA")

  x <- GSPip:::read_geno(path, format = "numeric", weight_col = "MarkerWeight")
  expect_identical(dim(x$M), c(3L, 3L))
  expect_identical(rownames(x$M), c("m1", "m2", "m3"))
  expect_true(is.na(x$M[3, 1]))
  expect_identical(unname(x$weights), c(1, 0.5, 2))
})

test_that("HapMap input decodes allele and IUPAC calls", {
  path <- tempfile(fileext = ".hmp.txt")
  on.exit(unlink(path), add = TRUE)
  x <- data.frame(
    `rs#` = c("h1", "h2", "h3"), alleles = c("A/G", "C/T", "G/A"),
    chrom = c(1, 1, 2), pos = c(10, 20, 30), strand = "+", `assembly#` = "v1",
    center = "x", protLSID = "x", assayLSID = "x", panelLSID = "x", QCcode = "x",
    S1 = c("AA", "CT", "NN"), S2 = c("R", "TT", "GG"), S3 = c("GG", "CC", "GA"),
    check.names = FALSE
  )
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NN")

  result <- GSPip:::read_geno(path, format = "hapmap")
  expect_identical(dim(result$M), c(3L, 3L))
  expect_identical(unname(result$M[1, ]), c(0, 1, 2))
  expect_true(is.na(result$M[3, 1]))
})

test_that("VCF input locates GT outside the first FORMAT position", {
  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", "S1", "S2", "S3"), collapse = "\t"),
    paste(c("1", "10", "v1", "A", "G", ".", "PASS", ".", "DP:GT",
            "12:0/0", "9:0/1", "8:1/1"), collapse = "\t"),
    paste(c("1", "20", "v2", "C", "T", ".", "PASS", ".", "DP:GT",
            "10:0|1", "11:1|1", "4:./."), collapse = "\t"),
    paste(c("2", "30", ".", "G", "A", ".", "PASS", ".", "DP:GT",
            "7:1/1", "5:0/0", "6:0/1"), collapse = "\t"),
    paste(c("2", "40", "multi", "A", "C,G", ".", "PASS", ".", "DP:GT",
            "7:0/2", "5:0/0", "6:0/1"), collapse = "\t")
  ), path)

  result <- GSPip:::read_geno(path, format = "vcf", multiallelic_action = "remove")
  expect_identical(dim(result$M), c(3L, 3L))
  expect_identical(unname(result$M[1, ]), c(0, 1, 2))
  expect_identical(unname(result$M[2, 1:2]), c(1, 2))
  expect_true(is.na(result$M[2, 3]))
  expect_true("2:30:G:A" %in% rownames(result$M))
  expect_identical(result$removed_multiallelic, "multi")
})

test_that("installed workflow templates can be located and copied", {
  source_path <- gspip_script("grm_builder")
  expect_true(file.exists(source_path))

  destination <- tempfile(fileext = ".R")
  on.exit(unlink(destination), add = TRUE)
  copied <- gspip_script("grm_builder", copy_to = destination)
  expect_true(file.exists(copied))
  expect_match(readLines(copied, n = 5L), "GSPip GRM workflow", all = FALSE)
})

