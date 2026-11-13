################################################################################
# GSPip GRM workflow template
#
# Edit SETTINGS, then run the whole script. For programmatic use, call
# GSPip::build_grm() directly.
################################################################################

library(GSPip)

# Files
geno_file    <- "path/to/genotypes.csv"
input_format <- "auto"       # "auto", "numeric", "hapmap", or "vcf"
output_dir   <- "path/to/results"
run_label    <- "GRM"

# Numeric-table columns (ignored for HapMap/VCF)
marker_col   <- "SNP"
chrom_col    <- "CHR"
position_col <- "POS"
ref_col      <- "REF"
alt_col      <- "ALT"
meta_cols    <- c("CHR", "POS", "REF", "ALT")
marker_weight_col <- NULL

# Optional external marker weights for every supported input format
marker_weights_file      <- NULL
marker_weights_id_col    <- "MarkerID"
marker_weights_value_col <- "Weight"

multiallelic_action <- "remove" # "remove" or "error"
grm_method <- "VanRaden"        # "VanRaden" or "Su"

# QC
qc_call_rate      <- 0.95
qc_maf            <- 0
qc_sweep_sample   <- 1
qc_heterozygosity <- 1

# Imputation
imputation_method <- "wright" # "none", "wright", "random_hwe", "rrblup_em", "beagle"
# IMPORTANT: "none" keeps NA calls and passes them to AGHmatrix.
qc_seed       <- 123
rrblup_em_tol <- 0.02
rrblup_cores  <- 1

# Beagle settings (used only for Beagle imputation)
beagle_jar        <- NULL
java_path         <- "java"
beagle_memory     <- "4g"
beagle_threads    <- 1
beagle_seed       <- 123
beagle_ne         <- NULL
beagle_ref        <- NULL
beagle_map        <- NULL
beagle_extra_args <- character(0)

write_imputed_genotypes <- FALSE
n_clusters    <- NULL
cluster_n_pcs <- 2
cluster_seed  <- 123
overwrite_outputs <- FALSE

result <- build_grm(
  geno_file = geno_file,
  output_dir = output_dir,
  input_format = input_format,
  run_label = run_label,
  marker_col = marker_col,
  chrom_col = chrom_col,
  position_col = position_col,
  ref_col = ref_col,
  alt_col = alt_col,
  meta_cols = meta_cols,
  marker_weight_col = marker_weight_col,
  marker_weights_file = marker_weights_file,
  marker_weights_id_col = marker_weights_id_col,
  marker_weights_value_col = marker_weights_value_col,
  multiallelic_action = multiallelic_action,
  grm_method = grm_method,
  qc_call_rate = qc_call_rate,
  qc_maf = qc_maf,
  qc_sweep_sample = qc_sweep_sample,
  qc_heterozygosity = qc_heterozygosity,
  imputation_method = imputation_method,
  qc_seed = qc_seed,
  rrblup_em_tol = rrblup_em_tol,
  rrblup_cores = rrblup_cores,
  beagle_jar = beagle_jar,
  java_path = java_path,
  beagle_memory = beagle_memory,
  beagle_threads = beagle_threads,
  beagle_seed = beagle_seed,
  beagle_ne = beagle_ne,
  beagle_ref = beagle_ref,
  beagle_map = beagle_map,
  beagle_extra_args = beagle_extra_args,
  write_imputed_genotypes = write_imputed_genotypes,
  n_clusters = n_clusters,
  cluster_n_pcs = cluster_n_pcs,
  cluster_seed = cluster_seed,
  overwrite_outputs = overwrite_outputs
)
print(result)

