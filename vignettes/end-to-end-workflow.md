# End-to-end genomic selection workflow

## Scope

This vignette covers the complete analysis sequence:

1. prepare genotype and phenotype inputs;
2. construct an additive or dominance genomic relationship matrix;
3. select a genomic training population with or without a predefined target set;
4. evaluate a single-record GBLUP model by repeated cross-validation;
5. evaluate a multi-environment GBLUP model under CV1 or CV2;
6. fit a fixed training set and predict a specified target population;
7. identify the output files used for reporting and downstream analysis.

All paths are relative to the repository directory. Absolute paths are also accepted.

## 1. Install the required software

Install the R packages used for data import, GRM construction, training-population optimization, and cross-validation:

```r
install.packages(c(
  "data.table",
  "AGHmatrix",
  "ASRgenomics",
  "rrBLUP",
  "STPGA",
  "cluster"
))
```

`rrBLUP` is required only for `imputation_method <- "rrblup_em"`.

GBLUP model fitting requires a licensed installation of ASReml-R. Confirm that the package loads:

```r
library(asreml)
```

Beagle imputation requires Java and a Beagle JAR file. No Beagle installation is required for the other imputation methods.

## 2. Create the analysis directories

Use separate directories for inputs and outputs:

```text
genomic-selection-workflow/
├── data/
│   ├── genotypes.csv
│   ├── stage1_blues.csv
│   ├── stage1_blues_by_environment.csv
│   ├── candidate_ids.csv
│   └── target_ids.csv
├── results/
├── R/
│   └── read_geno.R
├── grm_builder.R
├── optimize_training_population.R
├── gblup_cv_single.R
└── gblup_cv_gxe.R
```

Start R in the repository directory:

```r
setwd("path/to/genomic-selection-workflow")
```

Each analysis writes files with the prefix stored in `run_label`. A distinct label prevents accidental collisions between analyses.

## 3. Prepare genotype data

### 3.1 Numeric dosage table

Numeric input contains one marker per row and one genotype per sample column. Dosages are `0`, `1`, `2`, or missing.

```csv
SNP,CHR,POS,REF,ALT,G001,G002,G003,G004
S1_10234,1,10234,T,G,0,2,2,0
S1_18945,1,18945,G,T,0,0,0,2
S1_27656,1,27656,C,T,0,NA,2,1
S1_36367,1,36367,A,C,2,1,0,0
```

The following values are recognized as missing in numeric input: `NA`, `N/A`, an empty field, `.`, and `-9`.

Columns listed in `marker_col`, `chrom_col`, `position_col`, `ref_col`, `alt_col`, `meta_cols`, or `marker_weight_col` are treated as marker information. Every other column is treated as a genotype.

### 3.2 HapMap

HapMap input contains the 11 standard columns from `rs#` through `QCcode`, followed by sample columns. Allele letters, two-letter calls, IUPAC heterozygotes, and numeric `0/1/2` calls are accepted.

### 3.3 VCF

VCF input contains the nine fixed columns from `#CHROM` through `FORMAT`, followed by sample columns. The reader extracts diploid biallelic `GT` calls:

| GT | Dosage |
|---|---:|
| `0/0` or `0|0` | 0 |
| `0/1`, `1/0`, `0|1`, or `1|0` | 1 |
| `1/1` or `1|1` | 2 |

Multiallelic, haploid, and polyploid records are not converted to biallelic dosage. Records without a VCF ID receive the ID `CHROM:POS:REF:ALT`.

## 4. Construct an additive GRM

### 4.1 Configure the files

Open `grm_builder.R` and set:

```r
geno_file    <- "data/genotypes.csv"
input_format <- "numeric"
output_dir   <- "results"
run_label    <- "maize_additive"
```

For automatic format detection:

```r
input_format <- "auto"
```

The recognized extensions are `.vcf`, `.vcf.gz`, `.hmp`, `.hmp.txt`, `.hapmap`, and `.hapmap.txt`. Other extensions are read as numeric tables.

### 4.2 Configure numeric-table columns

```r
marker_col   <- "SNP"
chrom_col    <- "CHR"
position_col <- "POS"
ref_col      <- "REF"
alt_col      <- "ALT"
meta_cols    <- c("CHR", "POS", "REF", "ALT")
```

These settings are ignored for HapMap and VCF input.

### 4.3 Set marker quality-control thresholds

```r
qc_call_rate      <- 0.95
qc_maf            <- 0.05
qc_sweep_sample   <- 0.20
qc_heterozygosity <- 1
```

The filters are applied in the following order:

1. remove genotypes with missing fraction greater than `qc_sweep_sample`;
2. calculate marker call rate, allele frequency, minor allele frequency, and observed heterozygosity from the remaining observed calls;
3. remove markers that fail the configured thresholds;
4. remove markers with no observed calls or no allele variation;
5. impute the retained matrix with the selected method.

### 4.4 Select an imputation method

Wright imputation:

```r
imputation_method <- "wright"
qc_seed           <- 123
```

Hardy-Weinberg random imputation:

```r
imputation_method <- "random_hwe"
qc_seed           <- 123
```

EM imputation:

```r
imputation_method <- "rrblup_em"
rrblup_em_tol     <- 0.02
rrblup_cores      <- 1
```

No imputation:

```r
imputation_method <- "none"
```

With `"none"`, missing calls remain `NA` and are passed to `AGHmatrix::Gmatrix()`.

Beagle imputation:

```r
imputation_method <- "beagle"
beagle_jar        <- "software/beagle.jar"
java_path         <- "java"
beagle_memory     <- "4g"
beagle_threads    <- 4
beagle_seed       <- 123
beagle_ne         <- 100
beagle_ref        <- NULL
beagle_map        <- NULL
beagle_extra_args <- character(0)
```

Beagle requires chromosome, positive integer position, REF, and ALT values for every retained marker.

### 4.5 Select the GRM method

```r
grm_method <- "VanRaden"
```

### 4.6 Add marker weights when required

For an inline weight column in numeric input:

```r
marker_weight_col <- "MarkerWeight"
marker_weights_file <- NULL
```

For a separate marker-weight file:

```r
marker_weight_col          <- NULL
marker_weights_file        <- "data/marker_weights.csv"
marker_weights_id_col      <- "MarkerID"
marker_weights_value_col   <- "Weight"
```

The external table contains one finite, non-negative weight for every accepted biallelic marker before QC. Marker weights are supported only for the VanRaden GRM.

### 4.7 Configure optional outputs

```r
write_imputed_genotypes <- TRUE
n_clusters              <- 3
cluster_n_pcs           <- 5
cluster_seed            <- 123
overwrite_outputs       <- FALSE
```

Set `n_clusters <- NULL` when no genomic clustering output is required.

### 4.8 Run the analysis

```r
source("grm_builder.R")
```

### 4.9 Identify the GRM outputs

The principal outputs are:

| File | Content |
|---|---|
| `results/maize_additive_grm.csv` | Genomic relationship matrix |
| `results/maize_additive_sample_qc.csv` | Sample missingness and filter status |
| `results/maize_additive_marker_qc.csv` | Marker statistics and filter status |
| `results/maize_additive_marker_map.csv` | Retained marker map |
| `results/maize_additive_qc_report.txt` | QC counts and thresholds |
| `results/maize_additive_grm_info.txt` | Matrix diagnostics and inversion result |
| `results/maize_additive_clusters.csv` | GRM PCs and cluster assignments when requested |
| `results/maize_additive_imputed_genotypes.csv` | Final dosage matrix when requested |

The GRM CSV contains `GenoID` in the first column and one matrix column for each genotype.

## 5. Construct a dominance GRM

Use the genotype, QC, and output settings from Section 4, then set:

```r
grm_method <- "Su"
```

Disable marker weights:

```r
marker_weight_col  <- NULL
marker_weights_file <- NULL
```

Select `"none"`, `"wright"`, `"random_hwe"`, or `"beagle"` imputation. EM imputation is restricted to the additive VanRaden analysis.

Run:

```r
source("grm_builder.R")
```

The output structure is the same as for the additive analysis.

## 6. Prepare ID files for training-population optimization

Candidate and target files contain one unique ID per row in the column named by `id_col`:

```csv
GenoID
G001
G002
G003
```

Every ID must occur in the GRM. Candidate and target sets must be disjoint.

## 7. Optimize a training population for predefined targets

### 7.1 Configure the files

Open `optimize_training_population.R` and set:

```r
grm_file           <- "results/maize_additive_grm.csv"
candidate_ids_file <- "data/candidate_ids.csv"
target_ids_file    <- "data/target_ids.csv"
output_dir         <- "results"
run_label          <- "maize_targeted_tp"
id_col             <- "GenoID"
```

Set `candidate_ids_file <- NULL` to use every GRM genotype outside the target set as a candidate.

### 7.2 Configure the target design

```r
optimization_mode <- "targeted"
training_size     <- 200
criterion         <- "PEVMEAN"
lambda            <- 1e-6
```

The available criteria are:

```text
PEVMEAN, PEVMEAN2, PEVMAX, CDMEAN, CDMEAN2, CDMAX
```

All six criteria are minimized.

### 7.3 Configure GRM principal components

Automatic selection by cumulative positive-eigenvalue variance:

```r
stpga_n_pcs               <- NULL
stpga_pc_variance_target  <- 0.95
stpga_max_pcs             <- 100
```

Manual selection:

```r
stpga_n_pcs <- 25
```

### 7.4 Configure replicated genetic-algorithm searches

```r
ga_population         <- 100
ga_elite              <- 5
ga_iterations         <- 500
ga_min_before_stop    <- 200
ga_mutation_prob      <- 0.8
ga_mutation_intensity <- 1
ga_replicates         <- 20
ga_seed               <- 123
ga_cores              <- 1
```

On Windows, `ga_cores` must be `1`.

### 7.5 Configure the random benchmark

```r
random_replicates <- 100
random_seed       <- 456
overwrite_outputs <- FALSE
```

### 7.6 Run the optimization

```r
source("optimize_training_population.R")
```

### 7.7 Identify the optimization outputs

| File | Content |
|---|---|
| `results/maize_targeted_tp_selected_training.csv` | Selected training genotypes |
| `results/maize_targeted_tp_target_ids.csv` | Target genotypes |
| `results/maize_targeted_tp_candidate_ids.csv` | Candidate genotypes used by the search |
| `results/maize_targeted_tp_ga_replicates.csv` | Objective value for every GA replicate |
| `results/maize_targeted_tp_ga_training_sets.csv` | Membership of every GA solution |
| `results/maize_targeted_tp_ga_convergence.csv` | Objective trajectory for every GA replicate |
| `results/maize_targeted_tp_random_benchmark.csv` | Objectives for equal-size random sets |
| `results/maize_targeted_tp_random_training_sets.csv` | Membership of random sets |
| `results/maize_targeted_tp_pc_variance.csv` | GRM eigenvalues and PCs used |
| `results/maize_targeted_tp_run_info.txt` | Settings and result summary |

The selected training and target files can be supplied directly to a fixed-split GBLUP analysis.

## 8. Optimize without a predefined target set

Configure the optimizer as follows:

```r
grm_file           <- "results/maize_additive_grm.csv"
candidate_ids_file <- NULL
target_ids_file    <- NULL
output_dir         <- "results"
run_label          <- "maize_no_test_tp"
id_col             <- "GenoID"

optimization_mode <- "no_test"
training_size     <- 200
criterion         <- "PEVMEAN"
lambda            <- 1e-6
```

Use the PC, GA, benchmark, and overwrite settings from Section 7, then run:

```r
source("optimize_training_population.R")
```

The generated `_target_ids.csv` contains the complement of the selected training set within the candidate population.

## 9. Prepare single-record phenotype data

Single-record GBLUP requires one row per genotype:

```csv
GenoID,Trait,Population,SE,Weight
G001,10.66,P1,0.566,3.122
G002,10.92,P1,0.385,6.747
G003,10.38,P2,0.457,
G004,,P2,,
```

For an observed trait:

1. a positive finite value in `Weight` is used when present;
2. otherwise, a positive finite `SE` is converted to `1 / SE^2`;
3. the analysis stops if neither value is available.

Rows with a missing trait remain in the prediction table and require neither weight nor SE.

Every phenotype genotype must occur in the GRM.

## 10. Run repeated cross-validation for one record per genotype

### 10.1 Configure the files and columns

Open `gblup_cv_single.R` and set:

```r
grm_file   <- "results/maize_additive_grm.csv"
pheno_file <- "data/stage1_blues.csv"
output_dir <- "results"
run_label  <- "maize_single_cv"

id_col     <- "GenoID"
trait_col  <- "Trait"
strata_col <- "Population"
weight_col <- "Weight"
se_col     <- "SE"
```

### 10.2 Configure the weighted model

```r
fix_stage1_error  <- TRUE
nonmarker_effects <- TRUE
```

This specification fixes the stage-1 error variance for each record at `1 / weight` and includes the non-marker genotype effect.

For an unweighted single-record model:

```r
weight_col        <- NULL
fix_stage1_error  <- FALSE
nonmarker_effects <- FALSE
```

This specification gives all observed records equal weight and estimates the residual variance.

### 10.3 Configure repeated cross-validation

```r
evaluation_mode <- "repeated_cv"
k_folds         <- 5
n_iterations    <- 50
base_seed       <- 1000
```

Iteration `i` uses seed `base_seed + i`.

### 10.4 Select a stratification method

Supplied population or stratum:

```r
strata_col             <- "Population"
derive_strata_from_grm <- FALSE
```

GRM-derived strata:

```r
strata_col                  <- NULL
derive_strata_from_grm      <- TRUE
grm_strata_n_clusters       <- NULL
grm_strata_max_clusters     <- 10
grm_strata_n_pcs            <- 2
grm_strata_seed             <- 123
```

No genotype stratification:

```r
strata_col             <- NULL
derive_strata_from_grm <- FALSE
```

### 10.5 Configure stratum-specific metrics

```r
within_stratum_pa <- TRUE
min_stratum_n     <- 10
overwrite_outputs <- FALSE
```

### 10.6 Run the analysis

```r
source("gblup_cv_single.R")
```

### 10.7 Identify the principal outputs

| File | Content |
|---|---|
| `results/maize_single_cv_metrics.csv` | Fold- and iteration-level metrics |
| `results/maize_single_cv_summary.csv` | Mean and standard deviation across iterations |
| `results/maize_single_cv_cv_predictions.csv` | Out-of-sample predictions by genotype and iteration |
| `results/maize_single_cv_train_predictions.csv` | In-sample training predictions |
| `results/maize_single_cv_train_se.csv` | Standard errors of in-sample predictions |
| `results/maize_single_cv_intercepts.csv` | Model intercepts |
| `results/maize_single_cv_run_info.txt` | Settings and package versions |

Fold `0` in `_metrics.csv` contains metrics calculated from predictions pooled across folds within an iteration.

## 11. Prepare genotype-by-environment phenotype data

Multi-environment GBLUP requires one row per genotype × environment combination:

```csv
GenoID,Env,Trait,Population,SE,Weight
G001,E1,10.02,P1,0.400,6.250
G001,E2,11.68,P1,0.459,4.747
G002,E1,11.01,P1,0.448,4.982
G002,E2,,P1,,
```

Every observed trait requires a positive finite weight or a positive finite SE. An unweighted G × E model is not supported. When a stratum column is supplied, each genotype must retain the same stratum in every environment.

## 12. Run CV1 for untested genotypes

### 12.1 Configure files and columns

Open `gblup_cv_gxe.R` and set:

```r
grm_file   <- "results/maize_additive_grm.csv"
pheno_file <- "data/stage1_blues_by_environment.csv"
output_dir <- "results"
run_label  <- "maize_gxe_cv1"

id_col     <- "GenoID"
env_col    <- "Env"
trait_col  <- "Trait"
strata_col <- "Population"
weight_col <- "Weight"
se_col     <- "SE"
```

### 12.2 Configure the model

```r
fix_stage1_error  <- TRUE
nonmarker_effects <- TRUE
```

### 12.3 Configure CV1

```r
evaluation_mode <- "repeated_cv"
cv_scheme       <- "CV1"
k_folds         <- 5
n_iterations    <- 50
base_seed       <- 1000
```

CV1 assigns complete genotypes to folds. All environment records of a validation genotype are excluded from model fitting.

### 12.4 Configure accuracy

Use estimated Cullis heritability:

```r
env_h2 <- NULL
```

Or supply one heritability for every environment:

```r
env_h2 <- c(E1 = 0.62, E2 = 0.55, E3 = 0.58)
```

### 12.5 Configure stratification and run

Use one of the stratification settings in Section 10.4, then run:

```r
source("gblup_cv_gxe.R")
```

The file `_cv_predictions.csv` contains genotype-level out-of-sample predictions. The file `_cv_predictions_by_env.csv` contains the corresponding genotype × environment predictions.

## 13. Run CV2 for tested genotypes in unobserved environments

Use the G × E file, column, model, accuracy, and stratification settings from Section 12. Set:

```r
run_label       <- "maize_gxe_cv2"
evaluation_mode <- "repeated_cv"
cv_scheme       <- "CV2"
k_folds         <- 5
n_iterations    <- 50
base_seed       <- 1000
```

Run:

```r
source("gblup_cv_gxe.R")
```

CV2 assigns genotype × environment records to folds within environment and stratum. A genotype may therefore have records in both training and validation data. `Untested_val_rows` counts validation records whose genotype has no other observed record in the training data.

## 14. Predict a fixed target population with the single-record model

### 14.1 Prepare the split files

Use the optimizer outputs or prepare two CSV files with the configured `id_col`:

```text
results/maize_targeted_tp_selected_training.csv
results/maize_targeted_tp_target_ids.csv
```

Training and target IDs must be unique, disjoint, present in the GRM, and present in the phenotype table. Training genotypes require observed trait values. Target trait values may be observed or missing.

### 14.2 Configure the fixed split

In `gblup_cv_single.R`, set the files, columns, weights, and model as in Section 10, then set:

```r
run_label          <- "maize_single_fixed"
evaluation_mode    <- "fixed_split"
training_ids_file  <- "results/maize_targeted_tp_selected_training.csv"
target_ids_file    <- "results/maize_targeted_tp_target_ids.csv"
```

`k_folds`, `n_iterations`, and `base_seed` do not define a fixed split.

### 14.3 Run the model

```r
source("gblup_cv_single.R")
```

The file `_split_ids.csv` records the training or target role of every genotype. Metrics are calculated only for target genotypes with observed trait values. Predictions are produced for targets with missing trait values.

## 15. Predict a fixed target population with the G × E model

In `gblup_cv_gxe.R`, set the files, columns, weights, model, and accuracy settings as in Section 12. Then set:

```r
run_label          <- "maize_gxe_fixed"
evaluation_mode    <- "fixed_split"
cv_scheme          <- "CV1"
training_ids_file  <- "results/maize_targeted_tp_selected_training.csv"
target_ids_file    <- "results/maize_targeted_tp_target_ids.csv"
```

Run:

```r
source("gblup_cv_gxe.R")
```

Fixed-split G × E analysis is defined at genotype level and therefore requires CV1. All environment records of a target genotype are excluded from fitting.

## 16. Interpret GBLUP outputs

The GBLUP entry points write the following common outputs:

| Suffix | Content |
|---|---|
| `_metrics.csv` | Metrics for each iteration and fold |
| `_summary.csv` | Across-iteration mean and standard deviation |
| `_cv_predictions.csv` | Out-of-sample genomic predictions |
| `_cv_predictions_by_env.csv` | G × E predictions for individual environment records |
| `_train_predictions.csv` | In-sample fitted values for training genotypes |
| `_train_se.csv` | Standard errors of in-sample predictions |
| `_intercepts.csv` | Model intercepts |
| `_split_ids.csv` | Training and target membership in fixed-split mode |
| `_derived_strata.csv` | GRM PC coordinates and assigned strata |
| `_derived_strata_selection.csv` | Silhouette statistics for candidate cluster counts |
| `_derived_strata_pc_variance.csv` | GRM eigenvalues and PCs used for clustering |
| `_run_info.txt` | Analysis settings and package versions |

The principal metrics are:

| Metric | Definition |
|---|---|
| `PredictionAbility` | Correlation between observed values and out-of-sample predictions |
| `Accuracy` | Prediction ability divided by the square root of heritability |
| `RMSE` | Root mean squared out-of-sample prediction error |
| `Relative_RMSE` | RMSE as a percentage of the absolute trait mean |
| `Standardized_RMSE` | RMSE divided by the observed trait standard deviation |
| `Heritability` | Cullis heritability based on prediction-error variance of genotype differences |
| `PVE` | Proportion of variance attributed to the genomic main effect |
| `GOF` | In-sample correlation between observed and fitted values |
| `Converged` | ASReml convergence indicator |

For G × E analyses, prediction ability and GOF are calculated within environment and then averaged. Stage-1 weights enter model fitting but do not weight correlations or RMSE.

## 17. Re-run an analysis

An existing output prefix causes the analysis to stop when:

```r
overwrite_outputs <- FALSE
```

Use a new `run_label` for a separate analysis. Set `overwrite_outputs <- TRUE` only when replacement of the existing files is intended.
