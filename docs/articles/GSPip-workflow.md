# Complete GSPip workflow

## Complete GSPip workflow

### Objectives

GSPip integrates genotype processing, genomic relationship matrix (GRM)
construction, training-population optimisation, and genomic best linear
unbiased prediction (GBLUP). This tutorial demonstrates the complete
workflow with the example data distributed with the package.

The analyses address five use cases: (i) construction of additive and
dominance GRMs, (ii) optimisation of a training population for
predefined or unspecified targets, (iii) repeated cross-validation with
one stage-1 best linear unbiased estimate (BLUE) per genotype, (iv)
genotype-by-environment (G × E) prediction under CV1 and CV2, and (v)
prediction of a fixed target population.

The example set contains 120 genotypes, 240 markers, three populations,
and three environments. The numerical settings provide reproducible
examples; they are not universal inferential defaults.

### 1. Install and load GSPip

Install the package from GitHub:

``` r
install.packages("pak")
pak::pak("FAkohoue/GSPip")
```

Load GSPip and install the optional R packages used in this tutorial:

``` r
library(GSPip)
install.packages(c("STPGA", "BGLR", "rrBLUP"))
```

Selection of Training Populations by Genetic Algorithm (STPGA) is
required for training-population optimisation. Bayesian Generalized
Linear Regression (BGLR) is required for Bayesian prediction. `rrBLUP`
is required for expectation-maximisation imputation. ASReml-R is
licensed software and requires a separate installation. Beagle
imputation requires Java and a local Beagle JAR file.

### 2. Locate the example inputs

``` r
numeric_file <- system.file(
  "extdata", "example_genotypes_numeric.csv", package = "GSPip"
)
hapmap_file <- system.file(
  "extdata", "example_genotypes.hmp.txt", package = "GSPip"
)
vcf_file <- system.file(
  "extdata", "example_genotypes.vcf", package = "GSPip"
)
marker_weights_file <- system.file(
  "extdata", "example_marker_weights.csv", package = "GSPip"
)
single_blues_file <- system.file(
  "extdata", "example_BLUEs_single.csv", package = "GSPip"
)
gxe_blues_file <- system.file(
  "extdata", "example_BLUEs_gxe.csv", package = "GSPip"
)
candidate_ids_file <- system.file(
  "extdata", "example_candidate_ids.csv", package = "GSPip"
)
target_ids_file <- system.file(
  "extdata", "example_target_ids.csv", package = "GSPip"
)
example_grm_file <- system.file(
  "extdata", "example_grm.csv", package = "GSPip"
)

input_paths <- c(
  numeric_file, hapmap_file, vcf_file, marker_weights_file,
  single_blues_file, gxe_blues_file, candidate_ids_file,
  target_ids_file, example_grm_file
)
stopifnot(all(file.exists(input_paths)))
```

Create one output directory:

``` r
output_dir <- file.path(getwd(), "gspip-tutorial-results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
```

Every output file begins with `run_label`. Hence, each analysis below
uses a distinct label.

### 3. Inspect the input formats

#### Numeric genotypes

The numeric file contains one marker per row. The first five columns
contain marker information; the remaining 120 columns contain genotype
dosages.

``` r
numeric_genotypes <- data.table::fread(numeric_file)
dim(numeric_genotypes)
numeric_genotypes[1:5, 1:10]
```

Dosages are coded as `0`, `1`, `2`, or missing. Every column not
assigned to marker identification, coordinates, alleles, metadata, or
weights is interpreted as a genotype.

#### HapMap and Variant Call Format genotypes

The HapMap file contains the 11 standard fields from `rs#` through
`QCcode`, followed by sample columns. GSPip accepts allele letters,
two-letter calls, International Union of Pure and Applied Chemistry
(IUPAC) heterozygote codes, and numeric `0/1/2` calls.

The Variant Call Format (VCF) file contains the nine fixed fields
through `FORMAT`, followed by sample columns. Diploid biallelic `GT`
calls are converted to alternative-allele dosage. Unsupported records
are removed when `multiallelic_action = "remove"` and rejected when
`multiallelic_action = "error"`.

#### Phenotypes and ID files

``` r
single_blues <- data.table::fread(single_blues_file)
gxe_blues <- data.table::fread(gxe_blues_file)

names(single_blues)
names(gxe_blues)
dim(single_blues) # 120 genotypes
dim(gxe_blues)    # 120 genotypes × 3 environments
```

The single-record file contains one BLUE per genotype. The G × E file
contains one BLUE per genotype and environment. Candidate and target
files contain one unique `GenoID` per row.

### 4. Construct an additive GRM

#### Step 1: run quality control, imputation, and GRM construction

The following call: (i) removes genotypes with more than 20% missing
calls, (ii) retains markers with at least 95% call rate and 5% minor
allele frequency, (iii) imputes remaining missing calls with the Wright
method, and (iv) constructs an additive VanRaden GRM.

``` r
additive_fit <- build_grm(
  geno_file = numeric_file,
  output_dir = output_dir,
  input_format = "numeric",
  run_label = "additive_grm",
  marker_col = "SNP",
  chrom_col = "CHR",
  position_col = "POS",
  ref_col = "REF",
  alt_col = "ALT",
  meta_cols = c("CHR", "POS", "REF", "ALT"),
  multiallelic_action = "remove",
  grm_method = "VanRaden",
  qc_call_rate = 0.95,
  qc_maf = 0.05,
  qc_sweep_sample = 0.20,
  qc_heterozygosity = 1,
  imputation_method = "wright",
  qc_seed = 123,
  write_imputed_genotypes = TRUE,
  n_clusters = 3,
  cluster_n_pcs = 5,
  cluster_seed = 123,
  overwrite_outputs = FALSE
)

additive_fit
additive_fit$output_files
```

#### Step 2: examine the returned object

``` r
dim(additive_fit$grm)
additive_fit$grm[1:5, 1:5]
table(additive_fit$sample_qc$Keep)
table(additive_fit$marker_qc$Keep)
head(additive_fit$clusters)
```

The returned `gspip_grm` object contains the GRM, retained marker
matrix, marker map, quality-control tables, optional clusters, and
output paths.

#### Step 3: identify the output files

| Suffix                   | Content                                      |
|--------------------------|----------------------------------------------|
| `_grm.csv`               | Genomic relationship matrix                  |
| `_sample_qc.csv`         | Sample missingness and filter status         |
| `_marker_qc.csv`         | Marker statistics and filter status          |
| `_marker_map.csv`        | Retained marker coordinates and alleles      |
| `_qc_report.txt`         | Quality-control counts and thresholds        |
| `_grm_info.txt`          | Inversion result and matrix diagnostics      |
| `_clusters.csv`          | Principal components and cluster assignments |
| `_imputed_genotypes.csv` | Final marker matrix                          |

The GRM created above can be recovered with:

``` r
grm_file <- unname(additive_fit$output_files[["grm.csv"]])
```

The packaged GRM is used in the remaining sections so that each
procedure can be run independently:

``` r
grm_file <- example_grm_file
```

#### Alternative genotype formats

``` r
hapmap_fit <- build_grm(
  geno_file = hapmap_file,
  output_dir = output_dir,
  input_format = "hapmap",
  run_label = "hapmap_grm",
  grm_method = "VanRaden",
  imputation_method = "wright"
)

vcf_fit <- build_grm(
  geno_file = vcf_file,
  output_dir = output_dir,
  input_format = "vcf",
  run_label = "vcf_grm",
  grm_method = "VanRaden",
  imputation_method = "wright"
)
```

#### External marker weights

``` r
weighted_fit <- build_grm(
  geno_file = numeric_file,
  output_dir = output_dir,
  input_format = "numeric",
  run_label = "weighted_grm",
  marker_weights_file = marker_weights_file,
  marker_weights_id_col = "MarkerID",
  marker_weights_value_col = "Weight",
  grm_method = "VanRaden",
  imputation_method = "wright"
)
```

The external table must contain one finite, non-negative weight for
every accepted biallelic marker before quality control. The retained
weights are written to `_marker_weights_used.csv`.

#### Other imputation methods

Set `imputation_method` to `"random_hwe"`, `"rrblup_em"`, or `"none"`
for random Hardy-Weinberg, expectation-maximisation, or no imputation,
respectively.

Beagle is configured within
[`build_grm()`](https://FAkohoue.github.io/GSPip/reference/build_grm.md):

``` r
beagle_fit <- build_grm(
  geno_file = vcf_file,
  output_dir = output_dir,
  input_format = "vcf",
  run_label = "beagle_grm",
  grm_method = "VanRaden",
  imputation_method = "beagle",
  beagle_jar = "path/to/beagle.jar",
  java_path = "java",
  beagle_memory = "4g",
  beagle_threads = 4,
  beagle_seed = 123,
  beagle_ne = 100
)
```

Beagle requires chromosome, positive integer position, reference allele,
and alternative allele for every retained marker.

### 5. Construct a dominance GRM

The Su dominance GRM uses the same genotype and quality-control
arguments. Marker weights and expectation-maximisation imputation are
not supported for this method.

``` r
dominance_fit <- build_grm(
  geno_file = numeric_file,
  output_dir = output_dir,
  input_format = "numeric",
  run_label = "dominance_grm",
  grm_method = "Su",
  qc_call_rate = 0.95,
  qc_maf = 0.05,
  qc_sweep_sample = 0.20,
  imputation_method = "wright",
  qc_seed = 123
)
```

### 6. Optimise a training population for predefined targets

#### Step 1: run targeted optimisation

The example contains 90 candidates and 30 targets. The following
analysis selects 40 training genotypes.

``` r
targeted_fit <- optimize_training_population(
  grm_file = example_grm_file,
  output_dir = output_dir,
  candidate_ids_file = candidate_ids_file,
  target_ids_file = target_ids_file,
  run_label = "targeted_tp",
  id_col = "GenoID",
  optimization_mode = "targeted",
  training_size = 40,
  criterion = "PEVMEAN",
  lambda = 1e-6,
  stpga_n_pcs = NULL,
  stpga_pc_variance_target = 0.95,
  stpga_max_pcs = 100,
  ga_population = 100,
  ga_elite = 5,
  ga_iterations = 500,
  ga_min_before_stop = 200,
  ga_mutation_prob = 0.8,
  ga_mutation_intensity = 1,
  ga_replicates = 20,
  ga_seed = 123,
  ga_cores = 1,
  random_replicates = 100,
  random_seed = 456,
  overwrite_outputs = FALSE
)

targeted_fit
```

The available criteria are `PEVMEAN`, `PEVMEAN2`, `PEVMAX`, `CDMEAN`,
`CDMEAN2`, and `CDMAX`. All are minimised.

#### Step 2: examine convergence and random benchmarks

``` r
targeted_fit$replicates
head(targeted_fit$convergence)
summary(targeted_fit$random_benchmark$Objective)

best_objective <- min(targeted_fit$replicates$Objective)
mean(targeted_fit$random_benchmark$Objective <= best_objective)
```

The selected set should be interpreted with the convergence
trajectories, consistency among genetic-algorithm replicates, and the
equal-size random benchmark.

#### Step 3: retain the fixed-split paths

``` r
training_ids_file <- unname(
  targeted_fit$output_files[["selected_training.csv"]]
)
fixed_target_ids_file <- unname(
  targeted_fit$output_files[["target_ids.csv"]]
)
```

### 7. Optimise without a predefined target population

In `no_test` mode, the complement of each proposed training set forms
its objective population. The complement of the best solution is written
to `_target_ids.csv`.

``` r
no_test_fit <- optimize_training_population(
  grm_file = example_grm_file,
  output_dir = output_dir,
  candidate_ids_file = NULL,
  target_ids_file = NULL,
  run_label = "no_test_tp",
  optimization_mode = "no_test",
  training_size = 40,
  criterion = "PEVMEAN",
  stpga_pc_variance_target = 0.95,
  stpga_max_pcs = 100,
  ga_population = 100,
  ga_elite = 5,
  ga_iterations = 500,
  ga_min_before_stop = 200,
  ga_replicates = 20,
  ga_seed = 123,
  ga_cores = 1,
  random_replicates = 100,
  random_seed = 456
)
```

### 8. Run repeated cross-validation with one BLUE per genotype

#### ASReml analysis

The single-record ASReml model is $`y = \mu + g + r + e`$, where $`g`$
is the genomic effect, $`r`$ is the optional non-marker genotype effect,
and $`e`$ is the stage-1 error.

``` r
single_asreml <- run_gblup_single(
  grm_file = example_grm_file,
  pheno_file = single_blues_file,
  output_dir = output_dir,
  engine = "asreml",
  run_label = "single_asreml_cv",
  id_col = "GenoID",
  trait_col = "NBL",
  strata_col = "Population",
  weight_col = "Weight",
  se_col = "SE",
  derive_strata_from_grm = FALSE,
  fix_stage1_error = TRUE,
  nonmarker_effects = TRUE,
  accuracy_heritability = NULL,
  k_folds = 5,
  n_iterations = 50,
  base_seed = 1000,
  evaluation_mode = "repeated_cv",
  within_stratum_pa = TRUE,
  min_stratum_n = 10,
  overwrite_outputs = FALSE,
  asreml_ai_sing = FALSE
)

single_asreml$summary
```

With `fix_stage1_error = TRUE`, each stage-1 error variance is fixed at
`1 / weight`. The non-marker term is then estimable separately from the
known stage-1 error.

#### BGLR analysis

BGLR uses the GRM as a reproducing kernel Hilbert space (RKHS) kernel.
It estimates the residual scale and does not include the non-marker
identity term.

``` r
single_bglr <- run_gblup_single(
  grm_file = example_grm_file,
  pheno_file = single_blues_file,
  output_dir = output_dir,
  engine = "bglr",
  run_label = "single_bglr_cv",
  id_col = "GenoID",
  trait_col = "NBL",
  strata_col = "Population",
  weight_col = "Weight",
  se_col = "SE",
  fix_stage1_error = FALSE,
  nonmarker_effects = FALSE,
  k_folds = 5,
  n_iterations = 50,
  base_seed = 1000,
  evaluation_mode = "repeated_cv",
  bglr_n_iter = 12000,
  bglr_burn_in = 2000,
  bglr_thin = 10,
  bglr_seed = 123,
  bglr_R2 = 0.5,
  bglr_verbose = FALSE
)
```

#### Equal-weight analysis

When stage-1 precision columns are unavailable, residual scale
estimation and exclusion of non-marker effects are required:

``` r
weight_col <- NULL
fix_stage1_error <- FALSE
nonmarker_effects <- FALSE
```

Stage-1 precision information remains preferable when it is available.

#### Fold stratification

Use a supplied population or stratum:

``` r
strata_col <- "Population"
derive_strata_from_grm <- FALSE
```

Derive genomic strata from GRM principal components:

``` r
strata_col <- NULL
derive_strata_from_grm <- TRUE
grm_strata_n_clusters <- NULL
grm_strata_max_clusters <- 10
grm_strata_n_pcs <- 2
grm_strata_seed <- 123
```

With `grm_strata_n_clusters = NULL`, candidate cluster counts are
compared by average silhouette width. The resulting clusters are genomic
groups, not predefined biological populations.

Run without genotype stratification:

``` r
strata_col <- NULL
derive_strata_from_grm <- FALSE
```

### 9. Run G × E cross-validation

#### CV1 with ASReml

The multi-environment ASReml model is
$`y = \mu + E + g + r + gE + rE + e`$. CV1 assigns complete genotypes to
folds. Therefore, all environment records of a validation genotype are
excluded from fitting.

``` r
gxe_cv1_asreml <- run_gblup_gxe(
  grm_file = example_grm_file,
  pheno_file = gxe_blues_file,
  output_dir = output_dir,
  engine = "asreml",
  run_label = "gxe_cv1_asreml",
  id_col = "GenoID",
  trait_col = "NBL",
  env_col = "Env",
  strata_col = "Population",
  weight_col = "Weight",
  se_col = "SE",
  derive_strata_from_grm = FALSE,
  fix_stage1_error = TRUE,
  nonmarker_effects = TRUE,
  env_h2 = NULL,
  cv_scheme = "CV1",
  k_folds = 5,
  n_iterations = 50,
  base_seed = 1000,
  evaluation_mode = "repeated_cv",
  within_stratum_pa = TRUE,
  min_stratum_n = 10,
  overwrite_outputs = FALSE,
  asreml_ai_sing = FALSE
)
```

Environment-specific stage-1 heritabilities can be supplied for the
accuracy denominator:

``` r
env_h2 <- c(E1 = 0.62, E2 = 0.55, E3 = 0.58)
```

#### CV2 with ASReml

CV2 assigns genotype × environment records to folds within environment
and stratum. A genotype may therefore have records in both training and
validation data.

``` r
gxe_cv2_asreml <- run_gblup_gxe(
  grm_file = example_grm_file,
  pheno_file = gxe_blues_file,
  output_dir = output_dir,
  engine = "asreml",
  run_label = "gxe_cv2_asreml",
  id_col = "GenoID",
  trait_col = "NBL",
  env_col = "Env",
  strata_col = "Population",
  weight_col = "Weight",
  se_col = "SE",
  fix_stage1_error = TRUE,
  nonmarker_effects = TRUE,
  cv_scheme = "CV2",
  k_folds = 5,
  n_iterations = 50,
  base_seed = 1000,
  evaluation_mode = "repeated_cv"
)
```

`Untested_val_rows` counts validation records whose genotype has no
other observed record in the training data.

#### CV1 with BGLR

The BGLR G × E model contains one genomic main-effect kernel and one
genomic interaction kernel.

``` r
gxe_cv1_bglr <- run_gblup_gxe(
  grm_file = example_grm_file,
  pheno_file = gxe_blues_file,
  output_dir = output_dir,
  engine = "bglr",
  run_label = "gxe_cv1_bglr",
  id_col = "GenoID",
  trait_col = "NBL",
  env_col = "Env",
  strata_col = "Population",
  weight_col = "Weight",
  se_col = "SE",
  fix_stage1_error = FALSE,
  nonmarker_effects = FALSE,
  cv_scheme = "CV1",
  k_folds = 5,
  n_iterations = 50,
  base_seed = 1000,
  evaluation_mode = "repeated_cv",
  bglr_n_iter = 12000,
  bglr_burn_in = 2000,
  bglr_thin = 10,
  bglr_seed = 123,
  bglr_R2 = 0.5,
  bglr_verbose = FALSE
)
```

The same arguments apply to BGLR CV2 after setting `cv_scheme = "CV2"`
and assigning a distinct `run_label`.

### 10. Predict a fixed target population

#### One BLUE per genotype

Use the training and target files created by targeted optimisation.
Alternatively, supply independently defined files with the same `GenoID`
structure.

``` r
single_fixed <- run_gblup_single(
  grm_file = example_grm_file,
  pheno_file = single_blues_file,
  output_dir = output_dir,
  engine = "asreml",
  run_label = "single_fixed",
  id_col = "GenoID",
  trait_col = "NBL",
  strata_col = "Population",
  weight_col = "Weight",
  se_col = "SE",
  fix_stage1_error = TRUE,
  nonmarker_effects = TRUE,
  evaluation_mode = "fixed_split",
  training_ids_file = training_ids_file,
  target_ids_file = fixed_target_ids_file
)
```

Training and target IDs must be unique, disjoint, present in the GRM,
and present in the phenotype table. Training genotypes require observed
trait values. Target trait values may be observed for retrospective
assessment or missing for prediction.

#### G × E prediction

Fixed G × E prediction is defined at genotype level and therefore
requires CV1.

``` r
gxe_fixed <- run_gblup_gxe(
  grm_file = example_grm_file,
  pheno_file = gxe_blues_file,
  output_dir = output_dir,
  engine = "asreml",
  run_label = "gxe_fixed",
  id_col = "GenoID",
  trait_col = "NBL",
  env_col = "Env",
  strata_col = "Population",
  weight_col = "Weight",
  se_col = "SE",
  fix_stage1_error = TRUE,
  nonmarker_effects = TRUE,
  cv_scheme = "CV1",
  evaluation_mode = "fixed_split",
  training_ids_file = training_ids_file,
  target_ids_file = fixed_target_ids_file
)
```

All environment records of a target genotype are excluded from fitting.
`_split_ids.csv` records the role of every genotype.

### 11. Interpret the GBLUP outputs

Both GBLUP functions return the complete analysis in memory and write
auditable files.

``` r
names(single_asreml)
single_asreml$summary
head(single_asreml$metrics)
head(single_asreml$cv_predictions)
single_asreml$output_files
```

| Suffix | Content |
|----|----|
| `_metrics.csv` | Fold- and iteration-level metrics |
| `_summary.csv` | Mean and standard deviation across iterations |
| `_cv_predictions.csv` | Out-of-sample genotype predictions |
| `_cv_predictions_by_env.csv` | G × E predictions by genotype and environment |
| `_train_predictions.csv` | In-sample genomic predictions |
| `_train_se.csv` | ASReml standard errors or BGLR posterior standard deviations |
| `_intercepts.csv` | Model intercepts |
| `_split_ids.csv` | Fixed training and target membership |
| `_derived_strata.csv` | Genomic principal components and assigned strata |
| `_derived_strata_selection.csv` | Silhouette statistics for candidate cluster counts |
| `_derived_strata_pc_variance.csv` | GRM eigenvalues and principal components used for clustering |
| `_run_info.txt` | Settings, input counts, and package versions |

Fold `0` in `_metrics.csv` contains statistics calculated from
predictions pooled across all folds within an iteration.

`PredictionAbility` is the correlation between observed BLUEs and
out-of-sample predictions. `Accuracy` divides prediction ability by the
square root of heritability. For G × E analyses, prediction ability is
calculated within environment and then averaged.

`RMSE` is the root mean squared error of out-of-sample predictions.
`Relative_RMSE` expresses RMSE as a percentage of the absolute trait
mean. `Standardized_RMSE` divides RMSE by the observed standard
deviation; the G × E calculation first removes environment means.

ASReml analyses report residual maximum likelihood (REML) variance
components, proportion of variance explained (PVE), and Cullis
heritability when the fitted model supports these quantities. BGLR
analyses report posterior summaries with the suffix `_posterior`,
together with the deviance information criterion (DIC) and retained
Markov chain Monte Carlo (MCMC) sample count.

`Boundary_fit = 1` identifies an unstable ASReml variance solution.
Heritability and PVE are not reported for such a fit.

### 12. Use editable workflow templates

The exported functions form the package application programming
interface (API). GSPip also installs thin wrapper templates for analyses
managed through settings files.

``` r
dir.create("analysis", showWarnings = FALSE)
gspip_script("grm_builder", copy_to = "analysis/")
gspip_script("optimize_training_population", copy_to = "analysis/")
gspip_script("gblup_cv_single", copy_to = "analysis/")
gspip_script("gblup_cv_gxe", copy_to = "analysis/")
```

Each copied file calls the corresponding exported function. Hence,
direct function calls and template-based analyses use the same
implementation and produce the same output structure.

### 13. Preserve existing results

Every analysis stops when files with the selected `run_label` already
exist and `overwrite_outputs = FALSE`. A new `run_label` creates a
separate analysis. Existing files are replaced only when
`overwrite_outputs = TRUE` is set explicitly.

Therefore, the complete workflow can be reproduced from the input files,
function arguments, output files, and package versions recorded in
`_run_info.txt`.
