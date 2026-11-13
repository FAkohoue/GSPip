## data-raw/generate_example_data.R
## -----------------------------------------------------------------------------
## Generate reproducible example inputs for the genomic-selection workflow.
##
## Run from the repository root with:
##   Rscript data-raw/generate_example_data.R
## or:
##   source("data-raw/generate_example_data.R")
##
## The script uses base R only. Package datasets are written to data/ and flat
## example inputs to inst/extdata/. See R/example_data.R for documentation.
## -----------------------------------------------------------------------------

#===============================================================================
# SETTINGS
#===============================================================================
data_dir         <- "data"
extdata_dir      <- file.path("inst", "extdata")
overwrite_outputs <- FALSE
# For non-interactive regeneration, GS_EXAMPLE_OVERWRITE=true overrides the
# setting above. The default remains protective.
if (identical(tolower(Sys.getenv("GS_EXAMPLE_OVERWRITE", unset = "false")), "true")) {
  overwrite_outputs <- TRUE
}
simulation_seed  <- 20261106L

n_genotypes      <- 120L
n_chromosomes    <- 3L
blocks_per_chr   <- 4L
markers_per_block <- 20L
n_populations    <- 3L
n_environments   <- 3L
n_qtl            <- 18L
missing_call_rate <- 0.01
target_fraction  <- 0.25

#===============================================================================
# CHECK SETTINGS AND OUTPUTS
#===============================================================================
whole_number <- function(x, lower = 1L) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x >= lower && x == floor(x)
}

integer_settings <- list(
  simulation_seed = c(simulation_seed, 0L),
  n_genotypes = c(n_genotypes, 12L),
  n_chromosomes = c(n_chromosomes, 1L),
  blocks_per_chr = c(blocks_per_chr, 1L),
  markers_per_block = c(markers_per_block, 2L),
  n_populations = c(n_populations, 2L),
  n_environments = c(n_environments, 2L),
  n_qtl = c(n_qtl, 1L)
)
bad_integer <- names(integer_settings)[!vapply(integer_settings, function(z) {
  whole_number(z[1L], z[2L])
}, logical(1))]
if (length(bad_integer)) {
  stop("Invalid whole-number setting(s): ", paste(bad_integer, collapse = ", "),
       call. = FALSE)
}
if (n_populations > n_genotypes) {
  stop("n_populations cannot exceed n_genotypes.", call. = FALSE)
}
if (!is.numeric(missing_call_rate) || length(missing_call_rate) != 1L ||
    !is.finite(missing_call_rate) || missing_call_rate < 0 || missing_call_rate >= 0.05) {
  stop("missing_call_rate must be at least 0 and below 0.05 so examples pass default QC.",
       call. = FALSE)
}
if (!is.numeric(target_fraction) || length(target_fraction) != 1L ||
    !is.finite(target_fraction) || target_fraction <= 0 || target_fraction >= 1) {
  stop("target_fraction must be strictly between 0 and 1.", call. = FALSE)
}

flat_files <- c(
  "example_genotypes_numeric.csv",
  "example_genotypes.hmp.txt",
  "example_genotypes.vcf",
  "example_marker_weights.csv",
  "example_BLUEs_single.csv",
  "example_BLUEs_gxe.csv",
  "example_candidate_ids.csv",
  "example_target_ids.csv",
  "example_grm.csv",
  "example_simulation_truth.csv"
)
data_files <- c(
  "gs_example_geno.rda",
  "gs_example_marker_map.rda",
  "gs_example_marker_weights.rda",
  "gs_example_blues_single.rda",
  "gs_example_blues_gxe.rda",
  "gs_example_candidate_ids.rda",
  "gs_example_target_ids.rda",
  "gs_example_grm.rda",
  "gs_example_truth.rda",
  "gs_example_all.rda"
)
manifest_file <- "example_manifest.csv"
output_paths <- c(file.path(extdata_dir, flat_files),
                  file.path(data_dir, data_files),
                  file.path(extdata_dir, manifest_file))

if (!file.exists("DESCRIPTION") ||
    !any(grepl("^Package:[[:space:]]+GSPip[[:space:]]*$", readLines("DESCRIPTION", warn = FALSE)))) {
  stop("Run this script from the GSPip package root.", call. = FALSE)
}
for (path in c(data_dir, extdata_dir)) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output directory: ", path, call. = FALSE)
  }
  write_probe <- tempfile("example-data-write-probe-", tmpdir = path)
  if (!file.create(write_probe)) stop("Output directory is not writable: ", path, call. = FALSE)
  unlink(write_probe)
}
existing <- output_paths[file.exists(output_paths)]
if (!overwrite_outputs && length(existing)) {
  stop(
    "Example output(s) already exist: ", paste(existing, collapse = ", "),
    ". Remove them or set overwrite_outputs <- TRUE.",
    call. = FALSE
  )
}

set.seed(simulation_seed)

#===============================================================================
# HELPERS
#===============================================================================
standardize <- function(x) {
  z <- as.numeric(scale(x))
  if (any(!is.finite(z))) stop("Cannot standardize a constant simulated value.", call. = FALSE)
  z
}

write_csv <- function(x, filename) {
  utils::write.csv(
    x, file.path(extdata_dir, filename), row.names = FALSE,
    quote = FALSE, na = ""
  )
}

dosage_to_hapmap <- function(x, ref, alt) {
  out <- rep("NN", length(x))
  out[!is.na(x) & x == 0] <- paste0(ref, ref)
  out[!is.na(x) & x == 1] <- paste0(ref, alt)
  out[!is.na(x) & x == 2] <- paste0(alt, alt)
  out
}

dosage_to_vcf <- function(x) {
  out <- rep("./.", length(x))
  out[!is.na(x) & x == 0] <- "0/0"
  out[!is.na(x) & x == 1] <- "0/1"
  out[!is.na(x) & x == 2] <- "1/1"
  out
}

# Simulate one LD block from founder haplotypes. Populations use different
# founder frequencies, which creates population structure without making the
# populations genetically isolated.
simulate_block <- function(population, n_markers, n_founders = 8L,
                           mutation_rate = 0.015) {
  n <- length(population)
  founder_p <- stats::runif(n_markers, 0.15, 0.85)
  founders <- vapply(founder_p, function(p) {
    stats::rbinom(n_founders, 1L, p)
  }, numeric(n_founders))
  founders <- matrix(founders, nrow = n_founders, ncol = n_markers)

  # Avoid monomorphic founder columns before sampling diploid individuals.
  mono <- which(colSums(founders) %in% c(0L, n_founders))
  if (length(mono)) {
    for (j in mono) {
      founders[, j] <- sample(c(rep(0L, n_founders / 2L),
                                rep(1L, n_founders - n_founders / 2L)))
    }
  }

  pop_weights <- matrix(
    stats::rgamma(n_populations * n_founders, shape = 1.5, rate = 1),
    nrow = n_populations, ncol = n_founders
  )
  # Shared mass ensures that every founder can occur in every population.
  pop_weights <- pop_weights + 0.20
  pop_weights <- pop_weights / rowSums(pop_weights)

  hap1 <- hap2 <- integer(n)
  for (p in seq_len(n_populations)) {
    idx <- which(population == paste0("P", p))
    hap1[idx] <- sample.int(n_founders, length(idx), replace = TRUE,
                            prob = pop_weights[p, ])
    hap2[idx] <- sample.int(n_founders, length(idx), replace = TRUE,
                            prob = pop_weights[p, ])
  }
  dosage <- founders[hap1, , drop = FALSE] + founders[hap2, , drop = FALSE]

  mutate <- matrix(stats::runif(n * n_markers) < mutation_rate,
                   nrow = n, ncol = n_markers)
  direction <- matrix(sample(c(-1L, 1L), n * n_markers, replace = TRUE),
                      nrow = n, ncol = n_markers)
  dosage[mutate] <- pmin(2L, pmax(0L, dosage[mutate] + direction[mutate]))
  storage.mode(dosage) <- "integer"
  dosage
}

#===============================================================================
# GENOTYPES AND MARKER MAP
#===============================================================================
genotype_ids <- sprintf("G%03d", seq_len(n_genotypes))
population <- rep(sprintf("P%d", seq_len(n_populations)), length.out = n_genotypes)
population <- sample(population, n_genotypes, replace = FALSE)

allele_pairs <- matrix(c(
  "A", "C", "A", "G", "A", "T", "C", "G", "C", "T", "G", "T"
), ncol = 2L, byrow = TRUE, dimnames = list(NULL, c("REF", "ALT")))

genotype_parts <- list()
map_parts <- list()
marker_counter <- 0L

for (chr in seq_len(n_chromosomes)) {
  chr_position <- 10000L
  for (block in seq_len(blocks_per_chr)) {
    block_genotype <- simulate_block(population, markers_per_block)
    marker_index <- marker_counter + seq_len(markers_per_block)
    marker_ids <- sprintf("S%d_%04d", chr, seq_len(markers_per_block) +
                            (block - 1L) * markers_per_block)
    positions <- chr_position + cumsum(sample(500L:1500L, markers_per_block,
                                               replace = TRUE))
    pair_index <- ((marker_index - 1L) %% nrow(allele_pairs)) + 1L

    genotype_parts[[length(genotype_parts) + 1L]] <- block_genotype
    map_parts[[length(map_parts) + 1L]] <- data.frame(
      SNP = marker_ids,
      CHR = chr,
      POS = positions,
      REF = allele_pairs[pair_index, "REF"],
      ALT = allele_pairs[pair_index, "ALT"],
      Block = sprintf("C%dB%d", chr, block),
      stringsAsFactors = FALSE
    )
    marker_counter <- marker_counter + markers_per_block
    chr_position <- max(positions) + 100000L
  }
}

geno_complete <- do.call(cbind, genotype_parts)
marker_map <- do.call(rbind, map_parts)
rownames(geno_complete) <- genotype_ids
colnames(geno_complete) <- marker_map$SNP

# Repair a rare monomorphic marker so all generated markers are useful.
alt_frequency <- colMeans(geno_complete) / 2
monomorphic <- which(alt_frequency <= 0 | alt_frequency >= 1)
if (length(monomorphic)) {
  for (j in monomorphic) {
    p <- stats::runif(1L, 0.15, 0.45)
    geno_complete[, j] <- stats::rbinom(n_genotypes, 2L, p)
  }
}

# Add a controlled amount of missingness. Each marker receives at most the
# configured fraction, keeping it safely above grm_builder.R's default call-rate
# threshold of 0.95.
geno_observed <- geno_complete
n_missing_per_marker <- floor(n_genotypes * missing_call_rate)
if (n_missing_per_marker > 0L) {
  for (j in seq_len(ncol(geno_observed))) {
    geno_observed[sample.int(n_genotypes, n_missing_per_marker), j] <- NA_integer_
  }
}

stopifnot(
  identical(rownames(geno_observed), genotype_ids),
  identical(colnames(geno_observed), marker_map$SNP),
  all(geno_observed %in% c(0L, 1L, 2L, NA_integer_)),
  !anyDuplicated(marker_map$SNP),
  min(colMeans(!is.na(geno_observed))) > 0.95
)

# Numeric dosage: markers in rows, metadata first, then genotype columns.
numeric_genotypes <- cbind(
  marker_map[c("SNP", "CHR", "POS", "REF", "ALT")],
  as.data.frame(t(geno_observed), check.names = FALSE)
)
write_csv(numeric_genotypes, "example_genotypes_numeric.csv")

gs_example_geno <- geno_observed
gs_example_marker_map <- marker_map

# HapMap: 11 standard columns followed by genotype calls.
hapmap_calls <- matrix("NN", nrow = nrow(marker_map), ncol = n_genotypes,
                       dimnames = list(NULL, genotype_ids))
for (j in seq_len(nrow(marker_map))) {
  hapmap_calls[j, ] <- dosage_to_hapmap(
    geno_observed[, j], marker_map$REF[j], marker_map$ALT[j]
  )
}
hapmap <- cbind(
  data.frame(
    `rs#` = marker_map$SNP,
    alleles = paste(marker_map$REF, marker_map$ALT, sep = "/"),
    chrom = marker_map$CHR,
    pos = marker_map$POS,
    strand = "+",
    `assembly#` = "NA",
    center = "NA",
    protLSID = "NA",
    assayLSID = "NA",
    panelLSID = "NA",
    QCcode = "NA",
    check.names = FALSE,
    stringsAsFactors = FALSE
  ),
  as.data.frame(hapmap_calls, check.names = FALSE, stringsAsFactors = FALSE)
)
utils::write.table(
  hapmap, file.path(extdata_dir, "example_genotypes.hmp.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)

# VCF 4.2 with diploid GT calls.
vcf_header <- c(
  "##fileformat=VCFv4.2",
  "##source=GS-script-example-generator",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
          "INFO", "FORMAT", genotype_ids), collapse = "\t")
)
vcf_records <- character(nrow(marker_map))
for (j in seq_len(nrow(marker_map))) {
  vcf_records[j] <- paste(c(
    marker_map$CHR[j], marker_map$POS[j], marker_map$SNP[j],
    marker_map$REF[j], marker_map$ALT[j], ".", "PASS", ".", "GT",
    dosage_to_vcf(geno_observed[, j])
  ), collapse = "\t")
}
writeLines(c(vcf_header, vcf_records),
           file.path(extdata_dir, "example_genotypes.vcf"))

# External marker weights, deliberately non-uniform but unrelated to phenotype.
marker_weights <- data.frame(
  MarkerID = marker_map$SNP,
  Weight = round(stats::runif(nrow(marker_map), 0.5, 1.5), 6),
  stringsAsFactors = FALSE
)
write_csv(marker_weights, "example_marker_weights.csv")
gs_example_marker_weights <- marker_weights

#===============================================================================
# PHENOTYPES
#===============================================================================
# Phenotypes use the complete latent dosages. Missing observed genotype calls
# are a measurement/QC feature and do not change the simulated genetic truth.
qtl_index <- sort(sample.int(ncol(geno_complete), n_qtl))
qtl_effect <- stats::rnorm(n_qtl, mean = 0, sd = 0.35)
additive_value <- standardize(geno_complete[, qtl_index, drop = FALSE] %*% qtl_effect)

# A second effect vector produces genuine (but moderate) genotype-by-environment
# reranking rather than only environment mean shifts.
gxe_qtl_index <- sort(sample.int(ncol(geno_complete), n_qtl))
gxe_qtl_effect <- stats::rnorm(n_qtl, mean = 0, sd = 0.20)
gxe_value <- standardize(geno_complete[, gxe_qtl_index, drop = FALSE] %*% gxe_qtl_effect)

population_effect <- setNames(seq(-0.25, 0.25, length.out = n_populations),
                              sprintf("P%d", seq_len(n_populations)))
environment_ids <- sprintf("E%d", seq_len(n_environments))
environment_mean <- seq(-0.6, 0.6, length.out = n_environments)
environment_loading <- seq(-0.45, 0.45, length.out = n_environments)

gxe_rows <- vector("list", n_environments)
for (e in seq_len(n_environments)) {
  se <- stats::runif(n_genotypes, 0.25, 0.55)
  blue <- 10 + environment_mean[e] +
    1.20 * additive_value +
    environment_loading[e] * gxe_value +
    unname(population_effect[population]) +
    stats::rnorm(n_genotypes, 0, se)
  gxe_rows[[e]] <- data.frame(
    GenoID = genotype_ids,
    Env = environment_ids[e],
    NBL = round(blue, 4),
    Population = population,
    SE = round(se, 4),
    Weight = round(1 / se^2, 6),
    stringsAsFactors = FALSE
  )
}
blues_gxe <- do.call(rbind, gxe_rows)

# Demonstrate both supported precision fallbacks: some observed rows have only
# SE, while others have only Weight. Every observed row retains at least one.
se_only_rows <- sample.int(nrow(blues_gxe), floor(0.05 * nrow(blues_gxe)))
remaining_rows <- setdiff(seq_len(nrow(blues_gxe)), se_only_rows)
weight_only_rows <- sample(remaining_rows, floor(0.05 * nrow(blues_gxe)))
blues_gxe$Weight[se_only_rows] <- NA_real_
blues_gxe$SE[weight_only_rows] <- NA_real_

# A small number of missing cells exercise prediction-only G x E rows. Each
# genotype remains observed in at least two environments.
prediction_only_rows <- unlist(lapply(seq_len(n_genotypes), function(i) {
  if (i %% 20L == 0L) i + (sample.int(n_environments, 1L) - 1L) * n_genotypes else integer(0)
}), use.names = FALSE)
blues_gxe[prediction_only_rows, c("NBL", "SE", "Weight")] <- NA_real_

stopifnot(
  !anyDuplicated(paste(blues_gxe$GenoID, blues_gxe$Env, sep = "\r")),
  all(tapply(blues_gxe$Population, blues_gxe$GenoID,
             function(x) length(unique(x))) == 1L),
  all(is.na(blues_gxe$NBL) |
        (is.finite(blues_gxe$SE) & blues_gxe$SE > 0) |
        (is.finite(blues_gxe$Weight) & blues_gxe$Weight > 0))
)
write_csv(blues_gxe, "example_BLUEs_gxe.csv")
gs_example_blues_gxe <- blues_gxe

# One BLUE per genotype for gblup_cv_single.R. Precision is the inverse-variance
# combination of the available environment records.
single_rows <- lapply(genotype_ids, function(id) {
  rows <- blues_gxe[blues_gxe$GenoID == id & !is.na(blues_gxe$NBL), ]
  effective_weight <- ifelse(!is.na(rows$Weight), rows$Weight, 1 / rows$SE^2)
  combined_weight <- sum(effective_weight)
  data.frame(
    GenoID = id,
    NBL = round(stats::weighted.mean(rows$NBL, effective_weight), 4),
    Population = rows$Population[1L],
    SE = round(sqrt(1 / combined_weight), 4),
    Weight = round(combined_weight, 6),
    stringsAsFactors = FALSE
  )
})
blues_single <- do.call(rbind, single_rows)

# Also exercise SE fallback in the single-environment input.
single_se_only <- seq(10L, n_genotypes, by = 20L)
blues_single$Weight[single_se_only] <- NA_real_
write_csv(blues_single, "example_BLUEs_single.csv")
gs_example_blues_single <- blues_single

#===============================================================================
# TRAINING-POPULATION DESIGN INPUTS AND REFERENCE GRM
#===============================================================================
target_n <- max(1L, round(n_genotypes * target_fraction))
# Sample targets across populations so the example is not population-confounded.
target_ids <- unlist(lapply(split(genotype_ids, population), function(ids) {
  sample(ids, max(1L, round(length(ids) * target_fraction)))
}), use.names = FALSE)
if (length(target_ids) > target_n) target_ids <- sample(target_ids, target_n)
if (length(target_ids) < target_n) {
  target_ids <- c(target_ids, sample(setdiff(genotype_ids, target_ids),
                                    target_n - length(target_ids)))
}
target_ids <- genotype_ids[genotype_ids %in% target_ids]
candidate_ids <- setdiff(genotype_ids, target_ids)
gs_example_candidate_ids <- data.frame(GenoID = candidate_ids,
                                       stringsAsFactors = FALSE)
gs_example_target_ids <- data.frame(GenoID = target_ids,
                                    stringsAsFactors = FALSE)
write_csv(gs_example_candidate_ids, "example_candidate_ids.csv")
write_csv(gs_example_target_ids, "example_target_ids.csv")

# Mean-imputed VanRaden relationship matrix for immediate examples. The main
# workflow can instead build its own GRM from any of the genotype files above.
geno_for_grm <- geno_observed
for (j in seq_len(ncol(geno_for_grm))) {
  geno_for_grm[is.na(geno_for_grm[, j]), j] <- mean(geno_for_grm[, j], na.rm = TRUE)
}
allele_frequency <- colMeans(geno_for_grm) / 2
centered <- sweep(geno_for_grm, 2L, 2 * allele_frequency, FUN = "-")
denominator <- 2 * sum(allele_frequency * (1 - allele_frequency))
if (!is.finite(denominator) || denominator <= 0) {
  stop("The simulated markers cannot define a VanRaden GRM.", call. = FALSE)
}
example_grm <- tcrossprod(centered) / denominator
# A tiny diagonal ridge makes the illustrative matrix safely positive definite.
diag(example_grm) <- diag(example_grm) + 1e-6
rownames(example_grm) <- colnames(example_grm) <- genotype_ids
gs_example_grm <- example_grm
grm_output <- data.frame(GenoID = genotype_ids, example_grm,
                         check.names = FALSE, stringsAsFactors = FALSE)
write_csv(grm_output, "example_grm.csv")

# Simulation truth is for checking the examples, not for fitting the workflow.
truth <- merge(
  data.frame(MarkerID = marker_map$SNP, AdditiveEffect = 0,
             GxEEffect = 0, stringsAsFactors = FALSE),
  marker_map[c("SNP", "CHR", "POS", "Block")],
  by.x = "MarkerID", by.y = "SNP", sort = FALSE
)
truth$AdditiveEffect[match(marker_map$SNP[qtl_index], truth$MarkerID)] <- qtl_effect
truth$GxEEffect[match(marker_map$SNP[gxe_qtl_index], truth$MarkerID)] <- gxe_qtl_effect
write_csv(truth, "example_simulation_truth.csv")
gs_example_truth <- truth

#===============================================================================
# R-NATIVE OBJECTS, MANIFEST, AND FINAL CHECKS
#===============================================================================
# Individual .rda files are installed as package datasets and loaded with data().
save(gs_example_geno,
     file = file.path(data_dir, "gs_example_geno.rda"), compress = "xz")
save(gs_example_marker_map,
     file = file.path(data_dir, "gs_example_marker_map.rda"), compress = "xz")
save(gs_example_marker_weights,
     file = file.path(data_dir, "gs_example_marker_weights.rda"), compress = "xz")
save(gs_example_blues_single,
     file = file.path(data_dir, "gs_example_blues_single.rda"), compress = "xz")
save(gs_example_blues_gxe,
     file = file.path(data_dir, "gs_example_blues_gxe.rda"), compress = "xz")
save(gs_example_candidate_ids,
     file = file.path(data_dir, "gs_example_candidate_ids.rda"), compress = "xz")
save(gs_example_target_ids,
     file = file.path(data_dir, "gs_example_target_ids.rda"), compress = "xz")
save(gs_example_grm,
     file = file.path(data_dir, "gs_example_grm.rda"), compress = "xz")
save(gs_example_truth,
     file = file.path(data_dir, "gs_example_truth.rda"), compress = "xz")

# One-file convenience bundle. Loading this file creates gs_example_all, a
# named list containing every documented R object.
gs_example_all <- list(
  geno = gs_example_geno,
  marker_map = gs_example_marker_map,
  marker_weights = gs_example_marker_weights,
  blues_single = gs_example_blues_single,
  blues_gxe = gs_example_blues_gxe,
  candidate_ids = gs_example_candidate_ids,
  target_ids = gs_example_target_ids,
  grm = gs_example_grm,
  truth = gs_example_truth
)
save(gs_example_all,
     file = file.path(data_dir, "gs_example_all.rda"), compress = "xz")

manifest <- data.frame(
  File = c(file.path("inst", "extdata", flat_files),
           file.path("data", data_files)),
  Rows = c(
    nrow(numeric_genotypes), nrow(hapmap), nrow(marker_map), nrow(marker_weights),
    nrow(blues_single), nrow(blues_gxe), length(candidate_ids), length(target_ids),
    nrow(example_grm), nrow(truth),
    nrow(gs_example_geno), nrow(gs_example_marker_map),
    nrow(gs_example_marker_weights), nrow(gs_example_blues_single),
    nrow(gs_example_blues_gxe), nrow(gs_example_candidate_ids),
    nrow(gs_example_target_ids), nrow(gs_example_grm), nrow(gs_example_truth),
    length(gs_example_all)
  ),
  Purpose = c(
    "Numeric genotype input for grm_builder.R",
    "HapMap genotype input for grm_builder.R",
    "VCF genotype input for grm_builder.R",
    "Optional external marker weights",
    "One BLUE per genotype for gblup_cv_single.R",
    "One BLUE per genotype x environment for gblup_cv_gxe.R",
    "Optional STPGA candidate IDs",
    "STPGA targets or fixed-split targets",
    "Ready-to-use illustrative GRM",
    "Known simulated QTL effects (validation only)",
    "R-native genotype matrix",
    "R-native marker map",
    "R-native marker-weight table",
    "R-native single-environment BLUE table",
    "R-native G x E BLUE table",
    "R-native candidate-ID table",
    "R-native target-ID table",
    "R-native genomic relationship matrix",
    "R-native simulation truth",
    "Convenience list containing every R-native object"
  ),
  stringsAsFactors = FALSE
)
write_csv(manifest, "example_manifest.csv")

missing_outputs <- output_paths[!file.exists(output_paths)]
if (length(missing_outputs)) {
  stop("Failed to create: ", paste(missing_outputs, collapse = ", "), call. = FALSE)
}

message(
  "Example generation complete: ", n_genotypes, " genotypes, ",
  nrow(marker_map), " markers, ", n_environments, " environments."
)
message("Package datasets written to: ", normalizePath(data_dir, winslash = "/"))
message("Flat examples written to: ", normalizePath(extdata_dir, winslash = "/"))
message("See inst/extdata/example_manifest.csv for the file-to-workflow mapping.")
