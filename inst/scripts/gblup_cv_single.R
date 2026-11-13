################################################################################
# GSPip single-environment GBLUP workflow template
#
# Edit SETTINGS, then run the whole script. The analysis itself is implemented
# by GSPip::run_gblup_single().
################################################################################

library(GSPip)

# Files
grm_file   <- "path/to/GRM.csv"
pheno_file <- "path/to/stage1_BLUEs.csv"
output_dir <- "path/to/results"
run_label  <- "GBLUP_NBL"

# Columns
id_col     <- "GenoID"
trait_col  <- "NBL"
strata_col <- "Population" # NULL for none or GRM-derived strata
weight_col <- "Weight"     # May be NULL only with the no-weight settings below
se_col     <- "SE"

# Optional GRM-derived strata
derive_strata_from_grm  <- FALSE
grm_strata_n_clusters   <- NULL
grm_strata_max_clusters <- 10L
grm_strata_n_pcs        <- 2L
grm_strata_seed         <- 123L

# Model
fix_stage1_error  <- TRUE
nonmarker_effects <- TRUE
# With neither Weight nor SE, use:
# weight_col <- NULL; fix_stage1_error <- FALSE; nonmarker_effects <- FALSE

# Evaluation
k_folds      <- 5L
n_iterations <- 50L
base_seed    <- 0L
evaluation_mode   <- "repeated_cv" # "repeated_cv" or "fixed_split"
training_ids_file <- NULL
target_ids_file   <- NULL
within_stratum_pa <- TRUE
min_stratum_n     <- 10L

# Run
overwrite_outputs <- FALSE
asreml_workspace  <- 32e07

result <- run_gblup_single(
  grm_file = grm_file,
  pheno_file = pheno_file,
  output_dir = output_dir,
  run_label = run_label,
  id_col = id_col,
  trait_col = trait_col,
  strata_col = strata_col,
  weight_col = weight_col,
  se_col = se_col,
  derive_strata_from_grm = derive_strata_from_grm,
  grm_strata_n_clusters = grm_strata_n_clusters,
  grm_strata_max_clusters = grm_strata_max_clusters,
  grm_strata_n_pcs = grm_strata_n_pcs,
  grm_strata_seed = grm_strata_seed,
  fix_stage1_error = fix_stage1_error,
  nonmarker_effects = nonmarker_effects,
  k_folds = k_folds,
  n_iterations = n_iterations,
  base_seed = base_seed,
  evaluation_mode = evaluation_mode,
  training_ids_file = training_ids_file,
  target_ids_file = target_ids_file,
  within_stratum_pa = within_stratum_pa,
  min_stratum_n = min_stratum_n,
  overwrite_outputs = overwrite_outputs,
  asreml_workspace = asreml_workspace
)

print(result$summary)
