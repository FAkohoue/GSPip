#' Build a genomic relationship matrix
#'
#' Reads numeric, HapMap, or VCF genotypes; applies sample and marker QC;
#' optionally imputes missing calls; builds and checks a genomic relationship
#' matrix; and writes auditable result files.
#'
#' @param geno_file Path to the genotype file.
#' @param output_dir Directory for result files.
#' @param input_format One of "auto", "numeric", "hapmap", or "vcf".
#' @param run_label Prefix used for every output file.
#' @param marker_col,chrom_col,position_col,ref_col,alt_col Column names for
#'   numeric input.
#' @param meta_cols Additional non-genotype columns in numeric input.
#' @param marker_weight_col Optional inline marker-weight column for numeric input.
#' @param marker_weights_file Optional external marker-weight table.
#' @param marker_weights_id_col,marker_weights_value_col Columns in the external
#'   marker-weight table.
#' @param multiallelic_action Either "remove" or "error".
#' @param grm_method Either "VanRaden" or "Su".
#' @param qc_call_rate Minimum marker call rate.
#' @param qc_maf Minimum marker minor allele frequency.
#' @param qc_sweep_sample Maximum missing fraction per genotype.
#' @param qc_heterozygosity Maximum observed marker heterozygosity.
#' @param imputation_method One of "none", "wright", "random_hwe",
#'   "rrblup_em", or "beagle". With "none", missing calls are passed
#'   to AGHmatrix and do not cause GSPip to stop.
#' @param qc_seed Seed for stochastic imputation.
#' @param rrblup_em_tol,rrblup_cores rrBLUP EM controls.
#' @param beagle_jar,java_path,beagle_memory,beagle_threads,beagle_seed Beagle
#'   executable controls.
#' @param beagle_ne,beagle_ref,beagle_map,beagle_extra_args Optional Beagle
#'   population, reference, map, and additional arguments.
#' @param write_imputed_genotypes Write the final marker matrix when `TRUE`.
#' @param n_clusters Optional k-means cluster count.
#' @param cluster_n_pcs Number of leading GRM PCs used for clustering.
#' @param cluster_seed Clustering seed.
#' @param overwrite_outputs Permit replacement of files with the same run label.
#'
#' @return Invisibly, a \code{gspip_grm} list containing the GRM, retained
#'   marker matrix and metadata, QC tables, optional clusters, and output paths.
#' @export
build_grm <- function(
    geno_file,
    output_dir,
    input_format = "auto",
    run_label = "GRM",
    marker_col = "SNP",
    chrom_col = "CHR",
    position_col = "POS",
    ref_col = "REF",
    alt_col = "ALT",
    meta_cols = c("CHR", "POS", "REF", "ALT"),
    marker_weight_col = NULL,
    marker_weights_file = NULL,
    marker_weights_id_col = "MarkerID",
    marker_weights_value_col = "Weight",
    multiallelic_action = "remove",
    grm_method = "VanRaden",
    qc_call_rate = 0.95,
    qc_maf = 0,
    qc_sweep_sample = 1,
    qc_heterozygosity = 1,
    imputation_method = "wright",
    qc_seed = 123,
    rrblup_em_tol = 0.02,
    rrblup_cores = 1,
    beagle_jar = NULL,
    java_path = "java",
    beagle_memory = "4g",
    beagle_threads = 1,
    beagle_seed = 123,
    beagle_ne = NULL,
    beagle_ref = NULL,
    beagle_map = NULL,
    beagle_extra_args = character(0),
    write_imputed_genotypes = FALSE,
    n_clusters = NULL,
    cluster_n_pcs = 2,
    cluster_seed = 123,
    overwrite_outputs = FALSE
) {
  required <- c("data.table", "AGHmatrix", "ASRgenomics")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Please install required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (identical(imputation_method, "rrblup_em") &&
      !requireNamespace("rrBLUP", quietly = TRUE)) {
    stop("imputation_method = 'rrblup_em' requires the rrBLUP package.", call. = FALSE)
  }
  log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
  
  show_ids <- function(ids, n = 5L) {
    paste0(paste(utils::head(ids, n), collapse = ", "), if (length(ids) > n) ", ..." else "")
  }
  
  output_path <- function(suffix) file.path(output_dir, paste0(run_label, "_", suffix))
  
  is_nonempty_string <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }
  
  is_whole_number <- function(x, lower = 0) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= lower && x == floor(x)
  }
  
  check_probability <- function(x, name, upper = 1) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > upper) {
      stop(name, " must be between 0 and ", upper, ".", call. = FALSE)
    }
  }
  
  check_settings <- function() {
    if (!is_nonempty_string(geno_file) || !is_nonempty_string(output_dir)) {
      stop("geno_file and output_dir must be non-empty paths.", call. = FALSE)
    }
    if (!is_nonempty_string(run_label) || grepl('[<>:"/\\\\|?*]', run_label)) {
      stop("run_label contains an invalid file-name character.", call. = FALSE)
    }
    column_settings <- list(marker_col = marker_col, chrom_col = chrom_col,
                            position_col = position_col, ref_col = ref_col, alt_col = alt_col,
                            marker_weights_id_col = marker_weights_id_col,
                            marker_weights_value_col = marker_weights_value_col)
    if (any(!vapply(column_settings, is_nonempty_string, logical(1))) ||
        !is.character(meta_cols) || anyNA(meta_cols) || any(!nzchar(meta_cols))) {
      stop("Column-name settings and meta_cols must contain non-empty names.", call. = FALSE)
    }
    optional_strings <- list(marker_weight_col = marker_weight_col,
                             marker_weights_file = marker_weights_file,
                             beagle_ref = beagle_ref, beagle_map = beagle_map)
    invalid_optional <- names(optional_strings)[vapply(optional_strings, function(x) {
      !(is.null(x) || is_nonempty_string(x))
    }, logical(1))]
    if (length(invalid_optional)) {
      stop("These settings must be NULL or one non-empty string: ",
           paste(invalid_optional, collapse = ", "), call. = FALSE)
    }
    data_columns <- c(marker_col, chrom_col, position_col, ref_col, alt_col, marker_weight_col)
    if (anyDuplicated(data_columns)) {
      stop("marker_col, coordinate/allele columns, and marker_weight_col must be distinct.", call. = FALSE)
    }
    if (!is_nonempty_string(input_format) || !(input_format %in% c("auto", "numeric", "hapmap", "vcf"))) {
      stop("input_format must be 'auto', 'numeric', 'hapmap', or 'vcf'.", call. = FALSE)
    }
    if (!is_nonempty_string(multiallelic_action) || !(multiallelic_action %in% c("remove", "error"))) {
      stop("multiallelic_action must be 'remove' or 'error'.", call. = FALSE)
    }
    if (!is_nonempty_string(grm_method) || !(grm_method %in% c("VanRaden", "Su"))) {
      stop("grm_method must be 'VanRaden' or 'Su'.", call. = FALSE)
    }
    if (!is_nonempty_string(imputation_method) ||
        !(imputation_method %in% c("none", "wright", "random_hwe", "rrblup_em", "beagle"))) {
      stop("Unknown imputation_method. See the SETTINGS comments.", call. = FALSE)
    }
    check_probability(qc_call_rate, "qc_call_rate")
    check_probability(qc_maf, "qc_maf", 0.5)
    check_probability(qc_sweep_sample, "qc_sweep_sample")
    check_probability(qc_heterozygosity, "qc_heterozygosity")
    integer_settings <- list(qc_seed = qc_seed, rrblup_cores = rrblup_cores,
                             beagle_threads = beagle_threads, beagle_seed = beagle_seed,
                             cluster_seed = cluster_seed)
    bad_integer <- names(integer_settings)[!vapply(integer_settings, is_whole_number,
                                                    logical(1), lower = 0)]
    if (length(bad_integer)) {
      stop("Invalid non-negative whole-number setting(s): ", paste(bad_integer, collapse = ", "),
           call. = FALSE)
    }
    if (rrblup_cores < 1L || beagle_threads < 1L) {
      stop("rrblup_cores and beagle_threads must be at least 1.", call. = FALSE)
    }
    if (!is.numeric(rrblup_em_tol) || length(rrblup_em_tol) != 1L ||
        !is.finite(rrblup_em_tol) || rrblup_em_tol <= 0) {
      stop("rrblup_em_tol must be positive.", call. = FALSE)
    }
    if (!is.null(beagle_ne) && (!is.numeric(beagle_ne) || length(beagle_ne) != 1L ||
                                !is.finite(beagle_ne) || beagle_ne <= 0)) {
      stop("beagle_ne must be NULL or a positive number.", call. = FALSE)
    }
    if (!is.character(beagle_extra_args) || anyNA(beagle_extra_args)) {
      stop("beagle_extra_args must be a character vector without NA.", call. = FALSE)
    }
    reserved_beagle_args <- grepl("^(gt|out|nthreads|seed|ref|map|ne)=", beagle_extra_args,
                                  ignore.case = TRUE)
    if (any(reserved_beagle_args)) {
      stop("beagle_extra_args cannot override gt, out, nthreads, seed, ref, map, or ne; use the named settings.",
           call. = FALSE)
    }
    if (!is.null(n_clusters) && !is_whole_number(n_clusters, 2)) {
      stop("n_clusters must be NULL or a whole number of at least 2.", call. = FALSE)
    }
    if (!is_whole_number(cluster_n_pcs, 1)) stop("cluster_n_pcs must be positive.", call. = FALSE)
    scalar_logical <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
    if (!scalar_logical(write_imputed_genotypes) || !scalar_logical(overwrite_outputs)) {
      stop("write_imputed_genotypes and overwrite_outputs must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.null(marker_weight_col) && !is.null(marker_weights_file)) {
      stop("Use either marker_weight_col or marker_weights_file, not both.", call. = FALSE)
    }
    if ((!is.null(marker_weight_col) || !is.null(marker_weights_file)) && grm_method != "VanRaden") {
      stop("AGHmatrix marker weights are supported here only for the VanRaden GRM.", call. = FALSE)
    }
    if (imputation_method == "rrblup_em" && grm_method != "VanRaden") {
      stop("rrblup_em produces additive dosages and is supported here only with VanRaden.", call. = FALSE)
    }
    if (imputation_method == "beagle") {
      if (is.null(beagle_jar) || !is_nonempty_string(beagle_jar) || !file.exists(beagle_jar)) {
        stop("For Beagle imputation, beagle_jar must point to an existing Beagle .jar file.", call. = FALSE)
      }
      if (!is_nonempty_string(beagle_memory) || !grepl("^[1-9][0-9]*[kKmMgG]$", beagle_memory)) {
        stop("beagle_memory must look like '4g' or '800m'.", call. = FALSE)
      }
      if (!is_nonempty_string(java_path)) stop("java_path must be a non-empty command or path.", call. = FALSE)
      if (!is.null(beagle_ref) && !file.exists(beagle_ref)) {
        stop("beagle_ref not found: ", beagle_ref, call. = FALSE)
      }
      if (!is.null(beagle_map) && !file.exists(beagle_map)) {
        stop("beagle_map not found: ", beagle_map, call. = FALSE)
      }
    }
  }
  
  read_external_weights <- function(path, marker_ids) {
    if (!file.exists(path)) stop("Marker-weights file not found: ", path, call. = FALSE)
    x <- as.data.frame(fread(path, header = TRUE, check.names = FALSE))
    if (!all(c(marker_weights_id_col, marker_weights_value_col) %in% names(x))) {
      stop("Marker-weights file must contain columns '", marker_weights_id_col, "' and '",
           marker_weights_value_col, "'.", call. = FALSE)
    }
    ids <- as.character(x[[marker_weights_id_col]])
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
      stop("Marker IDs in the weights file must be non-missing and unique.", call. = FALSE)
    }
    weights <- suppressWarnings(as.numeric(x[[marker_weights_value_col]]))
    if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
      stop("External marker weights must be finite and non-negative.", call. = FALSE)
    }
    names(weights) <- ids
    missing_ids <- setdiff(marker_ids, ids)
    if (length(missing_ids)) stop("Weights are missing for marker(s): ", show_ids(missing_ids), call. = FALSE)
    weights[marker_ids]
  }
  
  #===============================================================================
  # INTERNAL QC
  #===============================================================================
  perform_qc <- function(M) {
    sample_observed <- colSums(!is.na(M))
    sample_missing <- 1 - sample_observed / nrow(M)
    keep_samples <- sample_missing <= qc_sweep_sample
    sample_qc <- data.frame(
      GenoID = colnames(M), MarkersObserved = sample_observed,
      MissingFraction = sample_missing, Kept = keep_samples,
      FailureReason = ifelse(keep_samples, "", "missing_fraction"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    if (sum(keep_samples) < 2L) {
      stop("Fewer than two genotypes remain after the sample-missingness filter.", call. = FALSE)
    }
    M2 <- M[, keep_samples, drop = FALSE]
    n_observed <- rowSums(!is.na(M2))
    call_rate <- n_observed / ncol(M2)
    allele_frequency <- rowMeans(M2, na.rm = TRUE) / 2
    allele_frequency[n_observed == 0L] <- NA_real_
    maf <- pmin(allele_frequency, 1 - allele_frequency)
    heterozygosity <- rowSums(M2 == 1, na.rm = TRUE) / n_observed
    heterozygosity[n_observed == 0L] <- NA_real_
    polymorphic <- !is.na(allele_frequency) & allele_frequency > 0 & allele_frequency < 1
    pass_observed <- n_observed > 0L
    pass_call_rate <- !is.na(call_rate) & call_rate >= qc_call_rate
    pass_maf <- !is.na(maf) & maf >= qc_maf
    pass_heterozygosity <- !is.na(heterozygosity) & heterozygosity <= qc_heterozygosity
    keep_markers <- pass_observed & pass_call_rate & pass_maf & pass_heterozygosity & polymorphic
    reasons <- vapply(seq_len(nrow(M2)), function(i) {
      failed <- c(if (!pass_observed[i]) "no_observed_calls",
                  if (!pass_call_rate[i]) "call_rate",
                  if (!pass_maf[i]) "maf",
                  if (!pass_heterozygosity[i]) "heterozygosity",
                  if (!polymorphic[i] && pass_observed[i]) "monomorphic")
      paste(failed, collapse = ";")
    }, character(1))
    marker_qc <- data.frame(
      MarkerID = rownames(M2), CallsObserved = n_observed, CallRate = call_rate,
      AlleleFrequency = allele_frequency, MAF = maf,
      ObservedHeterozygosity = heterozygosity, Polymorphic = polymorphic,
      Kept = keep_markers, FailureReason = reasons,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    if (sum(keep_markers) < 2L) stop("Fewer than two polymorphic markers remain after QC.", call. = FALSE)
    list(M = M2[keep_markers, , drop = FALSE], sample_qc = sample_qc,
         marker_qc = marker_qc, kept_samples = colnames(M2),
         kept_markers = rownames(M2)[keep_markers])
  }
  
  #===============================================================================
  # IMPUTATION
  #===============================================================================
  impute_random_hwe <- function(X) {
    p <- colMeans(X, na.rm = TRUE) / 2
    for (j in seq_len(ncol(X))) {
      missing <- which(is.na(X[, j]))
      if (length(missing)) X[missing, j] <- stats::rbinom(length(missing), 2L, p[j])
    }
    X
  }
  
  impute_wright <- function(X) {
    p <- colMeans(X, na.rm = TRUE) / 2
    expected_heterozygosity <- 2 * p * (1 - p)
    F_values <- rep(0, nrow(X))
    for (i in seq_len(nrow(X))) {
      observed <- !is.na(X[i, ])
      denominator <- sum(expected_heterozygosity[observed])
      if (denominator > 0) F_values[i] <- 1 - sum(X[i, observed] == 1) / denominator
    }
    F_values <- pmax(-0.99, pmin(0.99, F_values))
    for (j in seq_len(ncol(X))) {
      missing <- which(is.na(X[, j]))
      if (!length(missing)) next
      q <- 1 - p[j]
      Fm <- F_values[missing]
      probabilities <- cbind(q^2 + p[j] * q * Fm,
                             2 * p[j] * q * (1 - Fm),
                             p[j]^2 + p[j] * q * Fm)
      probabilities <- pmax(probabilities, 0)
      probabilities <- probabilities / rowSums(probabilities)
      X[missing, j] <- vapply(seq_along(missing), function(k) {
        sample.int(3L, 1L, prob = probabilities[k, ]) - 1L
      }, integer(1))
    }
    X
  }
  
  impute_rrblup_em <- function(X) {
    result <- tryCatch(
      rrBLUP::A.mat(X - 1, min.MAF = -1, max.missing = 1,
                    impute.method = "EM", tol = rrblup_em_tol,
                    n.core = rrblup_cores, return.imputed = TRUE),
      error = function(e) stop("rrBLUP EM imputation failed: ", conditionMessage(e), call. = FALSE)
    )
    if (!is.list(result) || is.null(result$imputed)) {
      stop("rrBLUP did not return an imputed marker matrix.", call. = FALSE)
    }
    out <- as.matrix(result$imputed) + 1
    if (!identical(dim(out), dim(X))) stop("rrBLUP changed the marker-matrix dimensions.", call. = FALSE)
    dimnames(out) <- dimnames(X)
    if (anyNA(out) || any(!is.finite(out))) stop("rrBLUP EM left non-finite values.", call. = FALSE)
    out
  }
  
  validate_beagle_metadata <- function(info, genotype_ids) {
    required <- c("MarkerID", "CHROM", "POS", "REF", "ALT")
    if (!all(required %in% names(info))) stop("Internal Beagle marker metadata are incomplete.", call. = FALSE)
    bad <- is.na(info$CHROM) | !nzchar(as.character(info$CHROM)) |
      is.na(info$POS) | !is.finite(info$POS) | info$POS <= 0 | info$POS != floor(info$POS) |
      is.na(info$REF) | !nzchar(as.character(info$REF)) |
      is.na(info$ALT) | !nzchar(as.character(info$ALT)) | grepl(",", info$ALT, fixed = TRUE)
    if (any(bad)) {
      stop("Beagle requires real CHROM, positive integer POS, REF, and one ALT allele for every retained marker. ",
           "Invalid marker(s): ", show_ids(info$MarkerID[bad]), call. = FALSE)
    }
    if (any(grepl("[\\t\\r\\n]", genotype_ids))) {
      stop("VCF sample IDs cannot contain tabs or line breaks.", call. = FALSE)
    }
  }
  
  write_beagle_vcf <- function(path, X, info) {
    validate_beagle_metadata(info, rownames(X))
    calls <- matrix("./.", nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))
    calls[!is.na(X) & X == 0] <- "0/0"
    calls[!is.na(X) & X == 1] <- "0/1"
    calls[!is.na(X) & X == 2] <- "1/1"
    table <- data.frame(
      `#CHROM` = as.character(info$CHROM),
      POS = format(info$POS, scientific = FALSE, trim = TRUE),
      ID = as.character(info$MarkerID), REF = as.character(info$REF),
      ALT = as.character(info$ALT), QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT",
      check.names = FALSE, stringsAsFactors = FALSE
    )
    table <- cbind(table, as.data.frame(t(calls), check.names = FALSE, stringsAsFactors = FALSE))
    writeLines(c("##fileformat=VCFv4.2", "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>"),
               path, useBytes = TRUE)
    fwrite(table, path, sep = "\t", quote = FALSE, na = ".", append = TRUE, col.names = TRUE)
  }
  
  run_beagle_imputation <- function(X, info) {
    work_dir <- tempfile(pattern = paste0(run_label, "_beagle_"), tmpdir = output_dir)
    if (!dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create temporary Beagle directory.", call. = FALSE)
    }
    on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
    input_vcf <- file.path(work_dir, "input.vcf")
    output_prefix <- file.path(work_dir, "imputed")
    # Beagle requires chromosome blocks with non-decreasing physical positions.
    marker_order <- order(as.character(info$CHROM), info$POS, info$MarkerID)
    write_beagle_vcf(input_vcf, X[, marker_order, drop = FALSE],
                     info[marker_order, , drop = FALSE])
    quoted <- function(x) shQuote(normalizePath(x, winslash = "/", mustWork = TRUE))
    args <- c(paste0("-Xmx", beagle_memory), "-jar", quoted(beagle_jar),
              shQuote(paste0("gt=", normalizePath(input_vcf, winslash = "/", mustWork = TRUE))),
              shQuote(paste0("out=", normalizePath(work_dir, winslash = "/", mustWork = TRUE),
                             "/imputed")),
              paste0("nthreads=", beagle_threads), paste0("seed=", beagle_seed))
    if (!is.null(beagle_ne)) args <- c(args, paste0("ne=", format(beagle_ne, scientific = FALSE)))
    if (!is.null(beagle_ref)) {
      if (!file.exists(beagle_ref)) stop("beagle_ref not found: ", beagle_ref, call. = FALSE)
      args <- c(args, shQuote(paste0("ref=", normalizePath(beagle_ref, winslash = "/",
                                                            mustWork = TRUE))))
    }
    if (!is.null(beagle_map)) {
      if (!file.exists(beagle_map)) stop("beagle_map not found: ", beagle_map, call. = FALSE)
      args <- c(args, shQuote(paste0("map=", normalizePath(beagle_map, winslash = "/",
                                                            mustWork = TRUE))))
    }
    args <- c(args, beagle_extra_args)
    log_msg("Running Beagle; see ", output_path("beagle.log"))
    command_output <- tryCatch(system2(java_path, args = args, stdout = TRUE, stderr = TRUE),
                               error = function(e) structure(conditionMessage(e), status = 1L))
    status <- attr(command_output, "status")
    beagle_log <- paste0(output_prefix, ".log")
    log_lines <- c(paste("Java command:", java_path, paste(args, collapse = " ")),
                   "", as.character(command_output))
    if (file.exists(beagle_log)) log_lines <- c(log_lines, "", readLines(beagle_log, warn = FALSE))
    writeLines(log_lines, output_path("beagle.log"))
    output_vcf <- paste0(output_prefix, ".vcf.gz")
    if ((!is.null(status) && status != 0L) || !file.exists(output_vcf)) {
      stop("Beagle failed (status ", if (is.null(status)) "unknown" else status,
           "). Review ", output_path("beagle.log"), ".", call. = FALSE)
    }
    result <- read_geno(output_vcf, format = "vcf",
                        multiallelic_action = multiallelic_action)
    missing_markers <- setdiff(colnames(X), rownames(result$M))
    missing_samples <- setdiff(rownames(X), colnames(result$M))
    if (length(missing_markers) || length(missing_samples)) {
      stop("Beagle output does not contain every original marker and genotype.", call. = FALSE)
    }
    out <- t(result$M[colnames(X), rownames(X), drop = FALSE])
    dimnames(out) <- dimnames(X)
    if (anyNA(out)) stop("Beagle left missing calls at the original markers.", call. = FALSE)
    out
  }
  
  impute_markers <- function(M, info) {
    X <- t(M)
    missing_before <- sum(is.na(X))
    if (missing_before == 0L || imputation_method == "none") {
      return(list(M = M, n_imputed = 0L, fractional = FALSE))
    }
    set.seed(qc_seed)
    X <- switch(imputation_method,
                wright = impute_wright(X),
                random_hwe = impute_random_hwe(X),
                rrblup_em = impute_rrblup_em(X),
                beagle = run_beagle_imputation(X, info))
    list(M = t(X), n_imputed = missing_before - sum(is.na(X)),
         fractional = identical(imputation_method, "rrblup_em"))
  }
  
  #===============================================================================
  # GRM HELPERS
  #===============================================================================
  check_inverse <- function(G) {
    result <- tryCatch(suppressMessages(suppressWarnings(G.inverse(G = G, sparseform = FALSE))),
                       error = function(e) e)
    if (inherits(result, "error")) {
      return(list(invertible = FALSE, note = gsub("\\s+", " ", conditionMessage(result))))
    }
    if (is.null(result$Ginv) || any(!is.finite(result$Ginv))) {
      return(list(invertible = FALSE, note = "no finite inverse returned"))
    }
    list(invertible = TRUE, note = "inverted")
  }
  
  clean_grm <- function(G, ids) {
    G <- as.matrix(G)
    G <- (G + t(G)) / 2
    dimnames(G) <- list(ids, ids)
    G
  }
  
  record_check <- function(step, check) {
    inverse_log <<- c(inverse_log, sprintf("G.inverse, %s: %s (%s)", step,
                                           if (check$invertible) "invertible" else "not invertible",
                                           check$note))
    log_msg(tail(inverse_log, 1L))
    check$invertible
  }
  
  #===============================================================================
  # RUN
  #===============================================================================
  check_settings()
  if (!file.exists(geno_file)) stop("Genotype file not found: ", geno_file, call. = FALSE)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output_dir: ", output_dir, call. = FALSE)
  }
  if (file.access(output_dir, 2) != 0) stop("output_dir is not writable: ", output_dir, call. = FALSE)
  
  output_suffixes <- c("grm.csv", "clusters.csv", "marker_weights_used.csv", "qc_report.txt",
                       "marker_qc.csv", "sample_qc.csv", "marker_map.csv",
                       "imputed_genotypes.csv", "beagle.log", "grm_info.txt")
  existing_outputs <- output_suffixes[file.exists(vapply(output_suffixes, output_path, character(1)))]
  if (!overwrite_outputs && length(existing_outputs)) {
    stop("Output file(s) already exist for run_label '", run_label, "': ",
         paste(existing_outputs, collapse = ", "),
         ". Change run_label or set overwrite_outputs <- TRUE.", call. = FALSE)
  }
  
  log_msg("Reading genotypes with the internal read_geno() function")
  marker_input <- read_geno(
    path = geno_file,
    format = input_format,
    marker_col = marker_col,
    chrom_col = chrom_col,
    position_col = position_col,
    ref_col = ref_col,
    alt_col = alt_col,
    meta_cols = meta_cols,
    weight_col = marker_weight_col,
    multiallelic_action = multiallelic_action
  )
  resolved_format <- marker_input$format
  M_input <- marker_input$M
  marker_info <- marker_input$info
  marker_weights <- marker_input$weights
  if (!is.null(marker_weights_file)) {
    marker_weights <- read_external_weights(marker_weights_file, rownames(M_input))
  }
  n_markers_input <- nrow(M_input)
  n_genotypes_input <- ncol(M_input)
  n_missing_input <- sum(is.na(M_input))
  if (n_genotypes_input < 2L) stop("At least two genotypes are required.", call. = FALSE)
  
  log_msg("Applying sample and marker QC")
  qc <- perform_qc(M_input)
  M_qc <- qc$M
  marker_info <- marker_info[match(rownames(M_qc), marker_info$MarkerID), , drop = FALSE]
  if (anyNA(marker_info$MarkerID)) stop("Could not align retained markers to their metadata.", call. = FALSE)
  weights_used <- if (is.null(marker_weights)) NULL else unname(marker_weights[rownames(M_qc)])
  if (!is.null(weights_used) && (length(weights_used) != nrow(M_qc) || anyNA(weights_used))) {
    stop("Could not align marker weights to retained markers.", call. = FALSE)
  }
  if (!is.null(weights_used) && !any(weights_used > 0)) {
    stop("All retained marker weights are zero.", call. = FALSE)
  }
  
  marker_qc_output <- qc$marker_qc
  if (length(marker_input$removed_multiallelic)) {
    removed_rows <- data.frame(
      MarkerID = marker_input$removed_multiallelic, CallsObserved = NA_integer_,
      CallRate = NA_real_, AlleleFrequency = NA_real_, MAF = NA_real_,
      ObservedHeterozygosity = NA_real_, Polymorphic = NA, Kept = FALSE,
      FailureReason = "multiallelic_or_unsupported",
      stringsAsFactors = FALSE, check.names = FALSE
    )
    marker_qc_output <- rbind(marker_qc_output, removed_rows)
  }
  fwrite(qc$sample_qc, output_path("sample_qc.csv"))
  fwrite(marker_qc_output, output_path("marker_qc.csv"))
  fwrite(marker_info, output_path("marker_map.csv"))
  if (length(marker_input$removed_multiallelic)) {
    log_msg(length(marker_input$removed_multiallelic), " multiallelic/unsupported marker(s) removed")
  }
  
  log_msg("Imputation method: ", imputation_method,
          if (imputation_method == "none") " (NA calls are retained for AGHmatrix)" else "")
  imputed <- impute_markers(M_qc, marker_info)
  M <- imputed$M
  if (imputation_method != "none" && anyNA(M)) {
    stop("The selected imputation method left missing genotype calls.", call. = FALSE)
  }
  
  if (write_imputed_genotypes) {
    fwrite(data.frame(MarkerID = rownames(M), M, check.names = FALSE),
           output_path("imputed_genotypes.csv"))
  }
  
  # Always pass impute.method = "none". With imputation_method = "none", NA calls
  # intentionally reach AGHmatrix; AGHmatrix handles them without this script stopping.
  log_msg("Computing the GRM: ", grm_method)
  G <- as.matrix(Gmatrix(SNPmatrix = t(M), method = grm_method,
                         integer = !imputed$fractional, impute.method = "none",
                         weights = weights_used))
  G <- (G + t(G)) / 2
  dimnames(G) <- list(colnames(M), colnames(M))
  if (any(!is.finite(G))) {
    stop("AGHmatrix returned missing or non-finite GRM values. Review marker QC and coding.", call. = FALSE)
  }
  
  # Test the computed matrix, then blend, then bend only if necessary.
  inverse_log <- character(0)
  G_computed <- G
  tuneup_applied <- "none"
  if (!record_check("computed GRM", check_inverse(G))) {
    G <- clean_grm(G.tuneup(G = G_computed, blend = TRUE, rcn = TRUE)$Gb, colnames(M))
    tuneup_applied <- "blend"
    if (!record_check("after blend", check_inverse(G))) {
      G <- clean_grm(G.tuneup(G = G_computed, bend = TRUE, rcn = TRUE)$Gb, colnames(M))
      tuneup_applied <- "bend"
      if (!record_check("after bend", check_inverse(G))) {
        stop("The GRM cannot be inverted after blending and bending: ",
             paste(inverse_log, collapse = "; "), call. = FALSE)
      }
    }
  }
  
  n_genotypes <- nrow(G)
  diag_mean <- mean(diag(G))
  offdiag_mean <- (sum(G) - sum(diag(G))) / (n_genotypes * (n_genotypes - 1))
  eig <- eigen(G, symmetric = TRUE)
  
  clusters <- NULL
  cluster_info <- character(0)
  if (!is.null(n_clusters)) {
    if (n_clusters > n_genotypes) stop("n_clusters exceeds the number of genotypes.", call. = FALSE)
    if (cluster_n_pcs > n_genotypes) stop("cluster_n_pcs exceeds the PCs available.", call. = FALSE)
    used <- seq_len(cluster_n_pcs)
    pcs <- sweep(eig$vectors[, used, drop = FALSE], 2L,
                 sqrt(pmax(eig$values[used], 0)), `*`)
    colnames(pcs) <- paste0("PC", used)
    if (nrow(unique(as.data.frame(pcs))) < n_clusters) {
      stop("n_clusters exceeds the distinct positions on the selected PCs.", call. = FALSE)
    }
    set.seed(cluster_seed)
    population <- kmeans(pcs, centers = n_clusters, nstart = 25)$cluster
    clusters <- data.frame(GenoID = rownames(G), as.data.frame(pcs, check.names = FALSE),
                           Population = population, check.names = FALSE,
                           stringsAsFactors = FALSE)
    fwrite(clusters, output_path("clusters.csv"))
    positive_sum <- sum(pmax(eig$values, 0))
    variance_percent <- if (positive_sum > 0) 100 * pmax(eig$values, 0) / positive_sum else rep(0, nrow(G))
    cluster_info <- c("", sprintf("Clusters: %d using %d PCs (cumulative variance %.1f%%), seed %s",
                                   n_clusters, cluster_n_pcs,
                                   sum(variance_percent[used]), cluster_seed),
                      capture.output(print(table(Population = population))))
  }
  
  fwrite(data.frame(GenoID = rownames(G), G, check.names = FALSE), output_path("grm.csv"))
  if (!is.null(weights_used)) {
    fwrite(data.frame(MarkerID = rownames(M), Weight = weights_used,
                      stringsAsFactors = FALSE), output_path("marker_weights_used.csv"))
  }
  
  qc_lines <- c(
    "QC is calculated from observed calls before imputation.",
    sprintf("Genotypes: %d input; %d retained; %d removed", n_genotypes_input,
            ncol(M_qc), n_genotypes_input - ncol(M_qc)),
    sprintf("Markers: %d biallelic input; %d retained; %d removed by QC", n_markers_input,
            nrow(M_qc), n_markers_input - nrow(M_qc)),
    sprintf("Unsupported/multiallelic records removed before QC: %d",
            length(marker_input$removed_multiallelic)),
    sprintf("Thresholds: marker call rate >= %s; MAF >= %s; sample missingness <= %s; marker heterozygosity <= %s",
            qc_call_rate, qc_maf, qc_sweep_sample, qc_heterozygosity),
    sprintf("Missing calls: %d in biallelic input; %d after QC; %d imputed; %d passed to AGHmatrix",
            n_missing_input, sum(is.na(M_qc)), imputed$n_imputed, sum(is.na(M)))
  )
  writeLines(qc_lines, output_path("qc_report.txt"))
  
  writeLines(c(
    paste("Run:", run_label),
    paste("Date:", format(Sys.time())),
    paste("Genotype file:", normalizePath(geno_file)),
    paste("Input format:", resolved_format),
    paste("GRM method:", grm_method),
    paste("Marker weights:", if (is.null(weights_used)) "equal (AGHmatrix default)" else "supplied"),
    paste("Imputation:", imputation_method),
    paste("Imputation seed:", if (imputation_method %in% c("wright", "random_hwe")) qc_seed else "not used"),
    if (imputation_method == "rrblup_em") {
      paste("rrBLUP EM tolerance:", rrblup_em_tol, "| cores:", rrblup_cores)
    },
    if (imputation_method == "beagle") {
      paste("Beagle seed:", beagle_seed, "| threads:", beagle_threads,
            "| memory:", beagle_memory, "| ne:", if (is.null(beagle_ne)) "Beagle default" else beagle_ne)
    },
    if (imputation_method == "beagle") {
      paste("Beagle reference:", if (is.null(beagle_ref)) "none" else normalizePath(beagle_ref),
            "| map:", if (is.null(beagle_map)) "none (1 cM/Mb)" else normalizePath(beagle_map))
    },
    paste("Calls passed as NA to AGHmatrix:", sum(is.na(M))),
    inverse_log,
    paste("Tune-up applied:", tuneup_applied),
    qc_lines,
    "",
    sprintf("Mean diagonal: %.4f", diag_mean),
    sprintf("Mean off-diagonal: %.4f", offdiag_mean),
    sprintf("Difference (Vg factor in Cullis heritability): %.4f", diag_mean - offdiag_mean),
    sprintf("Smallest eigenvalue: %.4g", min(eig$values)),
    cluster_info,
    "", capture.output(sessionInfo())
  ), output_path("grm_info.txt"))
  
  log_msg(sprintf("Finished: %d genotypes, %d markers; %d calls remain NA; tune-up %s",
                  n_genotypes, nrow(M), sum(is.na(M)), tuneup_applied))
  log_msg("Results are in: ", normalizePath(output_dir))
  result <- list(
    grm = G,
    markers = M,
    marker_map = marker_info,
    marker_weights = weights_used,
    sample_qc = qc$sample_qc,
    marker_qc = marker_qc_output,
    clusters = clusters,
    input_format = resolved_format,
    imputation_method = imputation_method,
    tuneup_applied = tuneup_applied,
    output_files = stats::setNames(
      vapply(output_suffixes, output_path, character(1)), output_suffixes
    )
  )
  result$output_files <- result$output_files[file.exists(result$output_files)]
  class(result) <- c("gspip_grm", "list")
  invisible(result)
}

#' @export
print.gspip_grm <- function(x, ...) {
  cat("GSPip genomic relationship matrix\\n")
  cat("  Genotypes :", nrow(x$grm), "\\n")
  cat("  Markers   :", nrow(x$markers), "\\n")
  cat("  Input     :", x$input_format, "\\n")
  cat("  Imputation:", x$imputation_method, "\\n")
  cat("  Tune-up   :", x$tuneup_applied, "\\n")
  invisible(x)
}

