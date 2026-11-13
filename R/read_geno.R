# Internal unified genotype reader.
#
# This file intentionally has no @export tag. In a package, read_geno() is an
# internal implementation detail. The standalone GRM script sources this file
# so that it uses exactly the same reader as the future package interface.
#
# Returned orientation is markers x samples, matching grm_builder.R:
#   $M                     numeric 0/1/2/NA matrix
#   $info                  MarkerID/CHROM/POS/REF/ALT data.frame
#   $weights               named numeric vector or NULL
#   $removed_multiallelic  character vector
#   $format                resolved input format

.gs_read_geno_api <- 1L

.geno_nonempty_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.geno_show_ids <- function(x, n = 5L) {
  paste0(paste(utils::head(x, n), collapse = ", "), if (length(x) > n) ", ..." else "")
}

.detect_genotype_format <- function(path, format) {
  choices <- c("auto", "numeric", "hapmap", "vcf")
  if (!.geno_nonempty_string(format) || !(format %in% choices)) {
    stop("format must be 'auto', 'numeric', 'hapmap', or 'vcf'.", call. = FALSE)
  }
  if (format != "auto") return(format)
  lower <- tolower(path)
  if (grepl("\\.vcf(\\.gz)?$", lower)) return("vcf")
  if (grepl("(\\.hmp(\\.txt)?|\\.hapmap(\\.txt)?)$", lower)) return("hapmap")
  "numeric"
}

.validate_genotype_ids <- function(marker_ids, sample_ids) {
  if (anyNA(marker_ids) || any(!nzchar(marker_ids)) || anyDuplicated(marker_ids)) {
    stop("Marker IDs must be non-missing, non-empty, and unique.", call. = FALSE)
  }
  if (length(sample_ids) < 2L || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) {
    stop("At least two uniquely named sample columns are required.", call. = FALSE)
  }
}

.coerce_numeric_dosage <- function(x) {
  original <- as.matrix(x)
  text <- matrix(as.character(original), nrow = nrow(original), dimnames = dimnames(original))
  missing <- is.na(original) | toupper(trimws(text)) %in% c("", ".", "NA", "N/A", "-9")
  suppressWarnings(M <- matrix(as.numeric(text), nrow = nrow(original), dimnames = dimnames(original)))
  M[missing] <- NA_real_
  invalid <- !missing & is.na(M)
  if (any(invalid)) {
    stop("Numeric input contains a non-numeric genotype value: ",
         text[which(invalid)[1L]], call. = FALSE)
  }
  if (any(!is.na(M) & !(M %in% c(0, 1, 2)))) {
    stop("Numeric genotype calls must be 0, 1, 2, or missing.", call. = FALSE)
  }
  M
}

.marker_info <- function(ids, chrom = NA_character_, pos = NA_real_,
                         ref = NA_character_, alt = NA_character_) {
  data.frame(
    MarkerID = ids,
    CHROM = rep_len(chrom, length(ids)),
    POS = rep_len(pos, length(ids)),
    REF = rep_len(ref, length(ids)),
    ALT = rep_len(alt, length(ids)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.read_numeric_genotypes <- function(path, marker_col, chrom_col, position_col,
                                    ref_col, alt_col, meta_cols, weight_col) {
  x <- as.data.frame(data.table::fread(
    path, header = TRUE, check.names = FALSE,
    na.strings = c("NA", "N/A", "", ".")
  ))
  if (anyDuplicated(names(x))) stop("Duplicated input column names.", call. = FALSE)
  if (!(marker_col %in% names(x))) stop("Marker ID column not found: ", marker_col, call. = FALSE)

  ids <- as.character(x[[marker_col]])
  non_genotype <- unique(c(marker_col, meta_cols, chrom_col, position_col,
                           ref_col, alt_col, weight_col))
  sample_ids <- setdiff(names(x), non_genotype)
  .validate_genotype_ids(ids, sample_ids)

  M <- .coerce_numeric_dosage(x[sample_ids])
  dimnames(M) <- list(ids, sample_ids)
  get_column <- function(name, default) {
    if (is.null(name) || !(name %in% names(x))) rep(default, nrow(x)) else x[[name]]
  }
  info <- .marker_info(
    ids,
    as.character(get_column(chrom_col, NA_character_)),
    suppressWarnings(as.numeric(get_column(position_col, NA_real_))),
    as.character(get_column(ref_col, NA_character_)),
    as.character(get_column(alt_col, NA_character_))
  )

  weights <- NULL
  if (!is.null(weight_col)) {
    if (!(weight_col %in% names(x))) stop("Marker-weight column not found: ", weight_col, call. = FALSE)
    weights <- suppressWarnings(as.numeric(x[[weight_col]]))
    if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
      stop("Marker weights must be non-missing, finite, and non-negative.", call. = FALSE)
    }
    weights <- stats::setNames(weights, ids)
  }

  list(M = M, info = info, weights = weights,
       removed_multiallelic = character(0))
}

.split_hapmap_alleles <- function(values) {
  text <- toupper(trimws(as.character(values)))
  parts <- strsplit(text, "[/|]", perl = TRUE)
  valid <- !is.na(text) & lengths(parts) == 2L &
    vapply(parts, function(z) length(z) == 2L && all(!is.na(z)) && all(nzchar(z)), logical(1))
  ref <- alt <- rep(NA_character_, length(parts))
  if (any(valid)) {
    ref[valid] <- vapply(parts[valid], `[[`, character(1), 1L)
    alt[valid] <- vapply(parts[valid], `[[`, character(1), 2L)
  }
  list(ref = ref, alt = alt, valid = valid)
}

.decode_hapmap_calls <- function(calls, ref, alt, marker_id) {
  z <- toupper(trimws(as.character(calls)))
  z[z %in% c("", ".", "-", "--", "N", "NN", "NA", "N/N", "./.")] <- NA_character_
  out <- rep(NA_real_, length(z))

  numeric_calls <- !is.na(z) & z %in% c("0", "1", "2")
  out[numeric_calls] <- as.numeric(z[numeric_calls])
  out[z %in% c(ref, paste0(ref, ref), paste0(ref, "/", ref), paste0(ref, "|", ref))] <- 0
  out[z %in% c(alt, paste0(alt, alt), paste0(alt, "/", alt), paste0(alt, "|", alt))] <- 2

  normalized <- gsub("[/|]", "", z)
  heterozygous <- !is.na(z) &
    nchar(normalized) == nchar(ref) + nchar(alt) &
    (normalized == paste0(ref, alt) | normalized == paste0(alt, ref))
  out[heterozygous] <- 1

  if (nchar(ref) == 1L && nchar(alt) == 1L) {
    pair <- paste(sort(c(ref, alt)), collapse = "")
    iupac <- c(AC = "M", AG = "R", AT = "W", CG = "S", CT = "Y", GT = "K")
    if (pair %in% names(iupac)) out[!is.na(z) & z == unname(iupac[pair])] <- 1
  }
  bad <- !is.na(z) & is.na(out)
  if (any(bad)) {
    stop("Unsupported HapMap call '", z[which(bad)[1L]], "' at marker ", marker_id,
         ". Expected allele calls, an IUPAC heterozygote, 0/1/2, or missing.",
         call. = FALSE)
  }
  out
}

.read_hapmap_genotypes <- function(path, multiallelic_action) {
  x <- as.data.frame(data.table::fread(
    path, header = TRUE, check.names = FALSE,
    na.strings = c("NA", "N", "NN", "", ".")
  ))
  expected <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center",
                "protLSID", "assayLSID", "panelLSID", "QCcode")
  lower_names <- tolower(names(x))
  if (anyDuplicated(lower_names)) stop("HapMap has duplicated column names ignoring case.", call. = FALSE)
  lookup <- stats::setNames(names(x), lower_names)
  missing_columns <- expected[!(tolower(expected) %in% names(lookup))]
  if (length(missing_columns)) {
    stop("HapMap input is missing standard column(s): ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  actual <- unname(lookup[tolower(expected)])
  sample_ids <- setdiff(names(x), actual)
  ids <- as.character(x[[actual[1L]]])
  .validate_genotype_ids(ids, sample_ids)

  alleles <- .split_hapmap_alleles(x[[actual[2L]]])
  unsupported <- !alleles$valid
  if (any(unsupported) && multiallelic_action == "error") {
    stop("HapMap contains invalid or multiallelic allele definitions: ",
         .geno_show_ids(ids[unsupported]), call. = FALSE)
  }
  keep <- !unsupported
  if (sum(keep) < 2L) stop("Fewer than two biallelic HapMap markers remain.", call. = FALSE)

  kept_rows <- which(keep)
  M <- matrix(NA_real_, nrow = length(kept_rows), ncol = length(sample_ids),
              dimnames = list(ids[keep], sample_ids))
  for (j in seq_along(kept_rows)) {
    i <- kept_rows[j]
    M[j, ] <- .decode_hapmap_calls(
      unlist(x[i, sample_ids, drop = FALSE], use.names = FALSE),
      alleles$ref[i], alleles$alt[i], ids[i]
    )
  }
  info <- .marker_info(
    ids[keep], as.character(x[[actual[3L]]][keep]),
    suppressWarnings(as.numeric(x[[actual[4L]]][keep])),
    alleles$ref[keep], alleles$alt[keep]
  )
  list(M = M, info = info, weights = NULL,
       removed_multiallelic = ids[unsupported])
}

.read_vcf_table <- function(path) {
  tryCatch(
    as.data.frame(data.table::fread(
      path, skip = "#CHROM", sep = "\t", header = TRUE,
      check.names = FALSE, data.table = FALSE, na.strings = "."
    )),
    error = function(e) stop("Could not read VCF: ", conditionMessage(e), call. = FALSE)
  )
}

.decode_vcf_calls <- function(sample_fields, format_fields, marker_ids) {
  formats <- strsplit(as.character(format_fields), ":", fixed = TRUE)
  gt_index <- vapply(formats, function(z) match("GT", z), integer(1))
  if (anyNA(gt_index)) {
    stop("VCF FORMAT lacks GT at marker ", marker_ids[which(is.na(gt_index))[1L]], ".",
         call. = FALSE)
  }
  out <- matrix(NA_real_, nrow = nrow(sample_fields), ncol = ncol(sample_fields),
                dimnames = dimnames(sample_fields))
  for (i in seq_len(nrow(sample_fields))) {
    raw <- as.character(unlist(sample_fields[i, ], use.names = FALSE))
    calls <- vapply(strsplit(raw, ":", fixed = TRUE), function(z) {
      if (length(z) < gt_index[i]) NA_character_ else z[gt_index[i]]
    }, character(1))
    calls[calls %in% c(".", "./.", ".|.", "0/.", "./0", "1/.", "./1",
                        "0|.", ".|0", "1|.", ".|1")] <- NA_character_
    out[i, calls %in% c("0/0", "0|0")] <- 0
    out[i, calls %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
    out[i, calls %in% c("1/1", "1|1")] <- 2
    bad <- !is.na(calls) & is.na(out[i, ])
    if (any(bad)) {
      stop("Unsupported VCF GT '", calls[which(bad)[1L]], "' at marker ", marker_ids[i],
           ". Only diploid biallelic GT is supported.", call. = FALSE)
    }
  }
  out
}

.read_vcf_genotypes <- function(path, multiallelic_action) {
  x <- .read_vcf_table(path)
  required <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  if (!all(required %in% names(x))) {
    stop("VCF must contain the nine standard fixed columns through FORMAT.", call. = FALSE)
  }
  format_index <- match("FORMAT", names(x))
  if (format_index >= ncol(x)) stop("VCF has no sample columns.", call. = FALSE)
  sample_ids <- names(x)[seq.int(format_index + 1L, ncol(x))]

  ids <- as.character(x[["ID"]])
  generated <- is.na(ids) | ids == "" | ids == "."
  ids[generated] <- paste(
    as.character(x[["#CHROM"]][generated]), as.character(x[["POS"]][generated]),
    as.character(x[["REF"]][generated]), as.character(x[["ALT"]][generated]), sep = ":"
  )
  .validate_genotype_ids(ids, sample_ids)

  alt <- as.character(x[["ALT"]])
  unsupported <- is.na(alt) | alt %in% c("", ".", "*") | grepl(",", alt, fixed = TRUE)
  if (any(unsupported) && multiallelic_action == "error") {
    stop("VCF contains multiallelic or unsupported ALT records: ",
         .geno_show_ids(ids[unsupported]), call. = FALSE)
  }
  keep <- !unsupported
  if (sum(keep) < 2L) stop("Fewer than two biallelic VCF markers remain.", call. = FALSE)

  M <- .decode_vcf_calls(
    as.matrix(x[keep, sample_ids, drop = FALSE]), x[["FORMAT"]][keep], ids[keep]
  )
  rownames(M) <- ids[keep]
  colnames(M) <- sample_ids
  info <- .marker_info(
    ids[keep], as.character(x[["#CHROM"]][keep]),
    suppressWarnings(as.numeric(x[["POS"]][keep])),
    as.character(x[["REF"]][keep]), alt[keep]
  )
  list(M = M, info = info, weights = NULL,
       removed_multiallelic = ids[unsupported])
}

read_geno <- function(path, format = "auto", marker_col = "SNP",
                      chrom_col = "CHR", position_col = "POS",
                      ref_col = "REF", alt_col = "ALT",
                      meta_cols = c("CHR", "POS", "REF", "ALT"),
                      weight_col = NULL,
                      multiallelic_action = c("remove", "error")) {
  if (!.geno_nonempty_string(path) || !file.exists(path)) {
    stop("Genotype file not found: ", path, call. = FALSE)
  }
  allowed_formats <- c("auto", "numeric", "hapmap", "vcf")
  if (!.geno_nonempty_string(format) || !(format %in% allowed_formats)) {
    stop("format must be 'auto', 'numeric', 'hapmap', or 'vcf'.", call. = FALSE)
  }
  columns <- list(marker_col = marker_col, chrom_col = chrom_col,
                  position_col = position_col, ref_col = ref_col, alt_col = alt_col)
  if (any(!vapply(columns, .geno_nonempty_string, logical(1))) ||
      !is.character(meta_cols) || anyNA(meta_cols) || any(!nzchar(meta_cols))) {
    stop("Column-name arguments and meta_cols must contain non-empty names.", call. = FALSE)
  }
  if (!(is.null(weight_col) || .geno_nonempty_string(weight_col))) {
    stop("weight_col must be NULL or one non-empty column name.", call. = FALSE)
  }
  if (anyDuplicated(c(unlist(columns, use.names = FALSE), weight_col))) {
    stop("Marker, coordinate, allele, and weight columns must be distinct.", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The internal genotype reader requires the data.table package.", call. = FALSE)
  }
  multiallelic_action <- match.arg(multiallelic_action)
  resolved <- .detect_genotype_format(path, format)
  if (!is.null(weight_col) && resolved != "numeric") {
    stop("An inline marker-weight column is supported only for numeric input.", call. = FALSE)
  }

  result <- switch(
    resolved,
    numeric = .read_numeric_genotypes(
      path, marker_col, chrom_col, position_col, ref_col, alt_col, meta_cols, weight_col
    ),
    hapmap = .read_hapmap_genotypes(path, multiallelic_action),
    vcf = .read_vcf_genotypes(path, multiallelic_action)
  )
  result$format <- resolved
  class(result) <- c("gs_genotypes", "list")
  result
}

