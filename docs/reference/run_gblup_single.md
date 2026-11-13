# Run single-environment GBLUP evaluation

Fits weighted genomic BLUP models from one stage-1 BLUE per genotype and
evaluates them by repeated cross-validation or a fixed training/target
split. ASReml-R is required only when this function is called.

## Usage

``` r
run_gblup_single(
  grm_file,
  pheno_file,
  output_dir,
  run_label = "GBLUP_NBL",
  id_col = "GenoID",
  trait_col = "NBL",
  strata_col = "Population",
  weight_col = "Weight",
  se_col = "SE",
  derive_strata_from_grm = FALSE,
  grm_strata_n_clusters = NULL,
  grm_strata_max_clusters = 10L,
  grm_strata_n_pcs = 2L,
  grm_strata_seed = 123L,
  fix_stage1_error = TRUE,
  nonmarker_effects = TRUE,
  k_folds = 5L,
  n_iterations = 50L,
  base_seed = 0L,
  evaluation_mode = c("repeated_cv", "fixed_split"),
  training_ids_file = NULL,
  target_ids_file = NULL,
  within_stratum_pa = TRUE,
  min_stratum_n = 10L,
  overwrite_outputs = FALSE,
  asreml_workspace = 3.2e+08
)
```

## Arguments

- grm_file:

  Path to a square GRM CSV whose first column contains genotype IDs.

- pheno_file:

  Path to the stage-1 BLUE table.

- output_dir:

  Directory receiving analysis outputs.

- run_label:

  Prefix used for output files.

- id_col, trait_col, strata_col, weight_col, se_col:

  Input column names. \`strata_col\` and \`weight_col\` may be \`NULL\`
  where documented.

- derive_strata_from_grm:

  Derive CV strata from GRM PCs.

- grm_strata_n_clusters:

  Optional fixed cluster count; \`NULL\` selects automatically.

- grm_strata_max_clusters, grm_strata_n_pcs, grm_strata_seed:

  GRM clustering controls.

- fix_stage1_error:

  Fix the stage-1 residual scale to one.

- nonmarker_effects:

  Include the non-marker genotype effect.

- k_folds, n_iterations, base_seed:

  Cross-validation controls.

- evaluation_mode:

  Either \`"repeated_cv"\` or \`"fixed_split"\`.

- training_ids_file, target_ids_file:

  ID files required for a fixed split.

- within_stratum_pa, min_stratum_n:

  Controls within-stratum prediction ability.

- overwrite_outputs:

  Permit replacement of outputs with the same run label.

- asreml_workspace:

  Workspace passed to ASReml-R.

## Value

Invisibly, a \`gspip_gblup_single\` list containing metrics,
predictions, input data used by the model, and output paths.
