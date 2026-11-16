# Empirical validation data

These files were validated with GSPip on 2026-11-15. They contain 204
genotypes and 9,126 biallelic dosage markers coded 0/1/2. All genotype,
phenotype, population, and GRM identifiers align exactly.

## Correct settings for these phenotype data

`Pheno_3k_forgblup.csv` has no `Weight`, `SE`, or `Env` column. It can therefore
be used with the single-environment workflow for traits `NBL` and `BL` using:

```r
weight_col        <- NULL
fix_stage1_error  <- FALSE
nonmarker_effects <- FALSE
```

It must not be passed to `run_gblup_gxe()` as-is: a G x E analysis requires an
environment column and a valid Weight or SE for every observed row.

The supplied `Population` column may be used as `strata_col`. Alternatively,
population-free validation was exercised with:

```r
strata_col              <- NULL
derive_strata_from_grm  <- TRUE
grm_strata_n_clusters   <- NULL
grm_strata_max_clusters <- 10L
grm_strata_n_pcs        <- 5L
```

Automatic silhouette selection chose three genomic clusters (average
silhouette 0.397). Those clusters did not reproduce the eight supplied
population labels closely, so they should be treated as fold-balancing groups,
not renamed as biological populations.

## Validation results

- `build_grm()` retained all 9,126 markers and all 204 genotypes with no
  imputation. The resulting blended GRM matched
  `unweighted_GRM_3k_blended.csv` with upper-triangle correlation 1.000000,
  RMSE 2.36e-19, and maximum absolute difference 2.78e-17.
- The supplied GRM is symmetric and positive definite; its smallest eigenvalue
  is 0.02000000044.
- Repeated 5-fold CV (two iterations) completed for NBL and BL with both ASReml
  and BGLR. Every expected output file existed and every CV and training
  prediction was finite.
- Mean out-of-sample prediction ability was 0.637 (ASReml) and 0.621 (BGLR) for
  NBL, and 0.725 (ASReml) and 0.712 (BGLR) for BL.
- Mean ASReml and BGLR CV predictions correlated 0.993 for NBL and 0.995 for BL.
- The automatic five-PC, unknown-cluster-count BGLR run completed with finite
  predictions and selected `K = 3`.

For these data, ASReml estimated residual variance very close to zero and
reported singular/negative-definite Average Information matrix warnings.
Validation therefore used `asreml_ai_sing = TRUE`. This is an explicit,
auditable exception, not a recommended default. GSPip flags these fits with
`Boundary_fit = 1` and reports heritability, PVE, and accuracy as `NA` because
no external stage-1 heritability was supplied, rather than presenting an
unstable near-one estimate. BGLR completed without model
warnings in this validation run.
