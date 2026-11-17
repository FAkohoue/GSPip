# Genomic selection workflow for R

This toolkit prepares genotype data, constructs genomic relationship matrices (GRMs), optimizes genomic training populations, and evaluates GBLUP models. It supports analyses with one phenotype per genotype and multi-environment analyses with genotype-by-environment (G × E) effects.

The complete workflow can:

- read numeric dosage tables, HapMap files, and diploid biallelic VCF/VCF.GZ files;
- filter samples and markers, impute missing genotypes, and apply optional marker weights;
- construct additive VanRaden or dominance Su GRMs;
- derive genomic clusters from GRM principal components;
- select an optimized training population with STPGA;
- run repeated cross-validation for one phenotype per genotype;
- run CV1 or CV2 for genotype-by-environment data;
- fit a fixed training set and predict a specified target population;
- report variance components, heritability, prediction ability, accuracy, and prediction errors.

The repository is run directly from R. Each analysis has a settings block at the beginning of its entry-point file; configure that block and source the file.

The [end-to-end workflow vignette](vignettes/end-to-end-workflow.md) gives complete procedures for GRM construction, training-population optimization, repeated cross-validation, G × E prediction, and fixed training/target analyses.

## Functions

### Genotype import and quality control

The genotype reader converts all supported formats to a marker-by-genotype dosage matrix containing `0`, `1`, `2`, or `NA`. Quality control is calculated from observed calls before imputation.

It can:

- remove samples above a missing-data threshold;
- filter markers by call rate, minor allele frequency, heterozygosity, and polymorphism;
- remove or reject unsupported multiallelic records;
- preserve marker coordinates and alleles for downstream use;
- accept inline or external marker weights.

### Genotype imputation

Five imputation choices are available:

| Method | Behavior |
|---|---|
| `"none"` | Leaves missing calls as `NA` and passes them to `AGHmatrix::Gmatrix()` |
| `"wright"` | Draws genotypes from marker allele frequencies and genotype-specific Wright inbreeding estimates |
| `"random_hwe"` | Draws genotypes from Hardy-Weinberg frequencies at each marker |
| `"rrblup_em"` | Uses `rrBLUP::A.mat()` multivariate EM imputation; additive VanRaden GRM only |
| `"beagle"` | Runs Beagle through Java using marker coordinates and alleles |

### Genomic relationship matrices

The GRM builder creates either:

- an additive VanRaden GRM; or
- a dominance Su GRM.

After construction, the matrix is checked for symmetry and invertibility. When necessary, blending and then bending are attempted. Optional k-means clustering uses a configurable number of leading GRM principal components.

### Training-population optimization

STPGA selects a fixed-size training set using genomic principal-component predictors from the GRM. Trait values are not used during selection.

Two designs are supported:

| Design | Target population |
|---|---|
| `"targeted"` | A supplied set of target genotypes |
| `"no_test"` | The complement of each proposed training set within the candidate population |

The available minimized design criteria are `PEVMEAN`, `PEVMEAN2`, `PEVMAX`, `CDMEAN`, `CDMEAN2`, and `CDMAX`. Replicated genetic-algorithm searches are compared with equal-size random training sets.

### Single-record GBLUP

For one stage-1 BLUE per genotype, the model is:

```text
y = mu + g + r + e
```

where `g` is the genomic effect with covariance proportional to the GRM, `r` is an optional non-marker genotype effect, and `e` is the stage-1 error. The analysis supports repeated k-fold cross-validation and a fixed training/target split.

### Multi-environment GBLUP

For one stage-1 BLUE per genotype × environment combination, the model is:

```text
y = mu + E + g + r + gE + rE + e
```

`E` is a fixed environment effect, `gE` is the genomic G × E effect, and `rE` is the optional non-marker G × E effect.

- CV1 holds out all records of a genotype and evaluates prediction of untested genotypes.
- CV2 holds out genotype × environment records and evaluates prediction in environments where a genotype has not been observed.

### Fixed training/target prediction

A selected or user-defined training population can be fitted once to predict a specified target set. Target phenotypes may be present for retrospective evaluation or missing for prediction-only use. For G × E data, fixed-split prediction uses CV1 because the split is defined at genotype level.

## Installation

### 1. Install R dependencies

Install the packages used for GRM construction, optimization, and data handling:

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

GBLUP analyses require ASReml-R, which is licensed software distributed by VSNi. Install ASReml-R and confirm that it loads before running a GBLUP analysis:

```r
library(asreml)
```

Beagle imputation additionally requires Java and a Beagle JAR file. Set `beagle_jar` to the downloaded JAR path and `java_path` to the Java executable.

### 2. Download the toolkit

Clone the repository or download and extract a release. Keep `R/read_geno.R` in the `R` directory beside the analysis entry points.

Open R in the repository directory, or set the working directory explicitly:

```r
setwd("path/to/genomic-selection-workflow")
```

### 3. Prepare an output directory

Each run writes files named `<run_label>_<output>.csv` or `<run_label>_<output>.txt`. Use a new `run_label` for each analysis. Existing files are preserved unless `overwrite_outputs <- TRUE`.

## Input formats

### Genotypes

Set `input_format` to `"numeric"`, `"hapmap"`, `"vcf"`, or `"auto"`. Automatic detection recognizes `.vcf`, `.vcf.gz`, `.hmp`, `.hmp.txt`, `.hapmap`, and `.hapmap.txt`; other extensions are treated as numeric tables.

| Format | Required layout | Dosage conversion |
|---|---|---|
| Numeric table | Markers in rows, one marker-ID column, and at least two genotype columns | Values are already `0`, `1`, `2`, or missing |
| HapMap | The 11 standard columns through `QCcode`, followed by sample columns | Counts the second allele listed in `alleles`; allele letters, two-letter calls, IUPAC heterozygotes, and `0/1/2` are accepted |
| VCF or VCF.GZ | Nine fixed VCF columns through `FORMAT`, followed by sample columns with diploid `GT` | `0/0 = 0`, heterozygote `= 1`, and `1/1 = 2` |

Only biallelic diploid markers are used. Set `multiallelic_action <- "remove"` to report and exclude unsupported records or `"error"` to stop. VCF records without an `ID` receive the generated ID `CHROM:POS:REF:ALT`.

A numeric dosage table may look like this:

| SNP | CHR | POS | REF | ALT | G001 | G002 | G003 |
|---|---:|---:|---|---|---:|---:|---:|
| S1_10234 | 1 | 10234 | T | G | 0 | 2 | 2 |
| S1_18945 | 1 | 18945 | G | T | 0 | 0 | 0 |
| S1_27656 | 1 | 27656 | C | T | 0 | NA | 2 |

Every column not listed as the marker column, a coordinate/allele column, a metadata column, or a marker-weight column is interpreted as a genotype. Marker IDs and genotype column names must be unique.

### Marker weights

Numeric input may contain an inline marker-weight column:

```r
marker_weight_col <- "MarkerWeight"
```

Alternatively, all genotype formats can use an external CSV containing one row per accepted marker:

| MarkerID | Weight |
|---|---:|
| S1_10234 | 1.25 |
| S1_18945 | 0.80 |

Weights must be finite and non-negative, and at least one retained weight must be positive. Marker weights are supported only with `grm_method <- "VanRaden"`.

### Genomic relationship matrix

The GRM input used by optimization and GBLUP must be a square, symmetric CSV. The first column contains genotype IDs, and the remaining column names contain the same IDs in the same or a reorderable set.

| GenoID | G001 | G002 | G003 |
|---|---:|---:|---:|
| G001 | 1.85 | -0.34 | -0.21 |
| G002 | -0.34 | 1.91 | -0.28 |
| G003 | -0.21 | -0.28 | 1.76 |

All entries must be finite. Row and column IDs must be unique and identical as sets.

### Candidate, training, and target IDs

ID files are CSV files with one unique, non-empty ID per row. The column name is configured with `id_col`, which defaults to `GenoID`.

```csv
GenoID
G001
G002
G003
```

Every ID must occur in the GRM. Candidate and target sets must be disjoint. Training and target sets used for fixed-split evaluation must also be disjoint.

### One phenotype per genotype

The single-record GBLUP input contains one row per genotype:

| GenoID | Trait | Population | SE | Weight |
|---|---:|---|---:|---:|
| G001 | 10.66 | P1 | 0.566 | 3.122 |
| G002 | 10.92 | P1 | 0.385 | 6.747 |
| G003 | 10.38 | P2 | 0.457 | NA |
| G004 | NA | P2 | NA | NA |

For observed traits, an existing positive `Weight` is used first; otherwise the analysis calculates `1 / SE^2`. Rows with a missing trait are retained for prediction and require neither `Weight` nor `SE`.

An unweighted single-record analysis is also supported. Use:

```r
weight_col        <- NULL
fix_stage1_error  <- FALSE
nonmarker_effects <- FALSE
```

This gives observed BLUEs equal weight and estimates the residual variance.

### Genotype-by-environment phenotypes

The G × E input contains one row per genotype × environment combination:

| GenoID | Env | Trait | Population | SE | Weight |
|---|---|---:|---|---:|---:|
| G001 | E1 | 10.02 | P1 | 0.400 | 6.250 |
| G001 | E2 | 11.68 | P1 | 0.459 | 4.747 |
| G002 | E1 | 11.01 | P1 | 0.448 | 4.982 |
| G002 | E2 | NA | P1 | NA | NA |

Every observed trait must have a positive finite weight or a positive finite SE. Unweighted G × E analysis is not supported. A genotype must have the same supplied stratum in every environment.

### Stratification

Cross-validation can use one of three fold-balancing strategies.

Use a supplied population or stratum column:

```r
strata_col             <- "Population"
derive_strata_from_grm <- FALSE
```

Derive genomic strata from the GRM:

```r
strata_col                 <- NULL
derive_strata_from_grm     <- TRUE
grm_strata_n_clusters      <- NULL  # select by average silhouette width
grm_strata_max_clusters    <- 10
grm_strata_n_pcs           <- 2
grm_strata_seed            <- 123
```

Run without genotype stratification:

```r
strata_col             <- NULL
derive_strata_from_grm <- FALSE
```

## Tutorials

### Use case 1: Build an additive GRM

This workflow reads genotypes, applies sample and marker QC, imputes missing calls, builds a VanRaden GRM, and writes diagnostics.

1. Prepare a numeric, HapMap, or VCF genotype file.
2. Open `grm_builder.R` and configure the input and output settings:

   ```r
   geno_file    <- "data/genotypes.csv"
   input_format <- "numeric"
   output_dir   <- "results"
   run_label    <- "maize_grm"

   marker_col   <- "SNP"
   chrom_col    <- "CHR"
   position_col <- "POS"
   ref_col      <- "REF"
   alt_col      <- "ALT"
   meta_cols    <- c("CHR", "POS", "REF", "ALT")
   ```

3. Choose the GRM and QC settings:

   ```r
   grm_method        <- "VanRaden"
   qc_call_rate      <- 0.95
   qc_maf            <- 0.05
   qc_sweep_sample   <- 0.20
   qc_heterozygosity <- 1
   ```

4. Choose an imputation method and seed:

   ```r
   imputation_method <- "wright"
   qc_seed           <- 123
   ```

5. Optionally request genomic clusters and the final dosage matrix:

   ```r
   n_clusters              <- 3
   cluster_n_pcs           <- 5
   write_imputed_genotypes <- TRUE
   ```

6. Run the analysis:

   ```r
   source("grm_builder.R")
   ```

7. Use `results/maize_grm_grm.csv` as the GRM input for optimization or GBLUP. Review the sample and marker QC tables before using the matrix.

GRM outputs:

| Suffix | Content |
|---|---|
| `_grm.csv` | GRM with `GenoID` in the first column |
| `_sample_qc.csv` | Missingness, keep flag, and failure reason for every genotype |
| `_marker_qc.csv` | Marker statistics, keep flag, and failure reason |
| `_marker_map.csv` | Retained marker IDs, chromosome, position, REF, and ALT |
| `_qc_report.txt` | QC counts and thresholds |
| `_grm_info.txt` | Input format, imputation counts, inversion checks, tune-up, and matrix diagnostics |
| `_clusters.csv` | GRM PC coordinates and population assignment when `n_clusters` is set |
| `_marker_weights_used.csv` | Retained marker weights when weights are supplied |
| `_imputed_genotypes.csv` | Final marker-by-genotype matrix when requested |
| `_beagle.log` | Beagle command and output when Beagle is used |

### Use case 2: Build a dominance GRM

1. Configure the genotype file, output directory, QC, and imputation as in Use case 1.
2. Select the Su dominance method:

   ```r
   grm_method <- "Su"
   ```

3. Do not supply inline or external marker weights. If using EM imputation, choose another imputation method because `rrblup_em` is restricted to the additive GRM.
4. Run:

   ```r
   source("grm_builder.R")
   ```

5. Use the generated `_grm.csv` and review `_grm_info.txt` for the matrix diagnostics.

### Use case 3: Optimize a training population for known targets

This workflow selects candidate genotypes that minimize a genomic design criterion for a supplied target population.

1. Prepare a GRM and a target-ID CSV. Optionally prepare a candidate-ID CSV; when omitted, every GRM genotype outside the target set is eligible.
2. Open `optimize_training_population.R` and configure the files:

   ```r
   grm_file           <- "results/maize_grm_grm.csv"
   candidate_ids_file <- "data/candidate_ids.csv"  # or NULL
   target_ids_file    <- "data/target_ids.csv"
   output_dir         <- "results"
   run_label          <- "maize_targeted_tp"
   id_col             <- "GenoID"
   ```

3. Configure the design:

   ```r
   optimization_mode <- "targeted"
   training_size     <- 200
   criterion         <- "PEVMEAN"
   lambda            <- 1e-6
   ```

4. Select the GRM PCs automatically or set `stpga_n_pcs` to a positive integer:

   ```r
   stpga_n_pcs              <- NULL
   stpga_pc_variance_target <- 0.95
   stpga_max_pcs            <- 100
   ```

5. Configure replicated optimization and the random benchmark:

   ```r
   ga_replicates     <- 20
   ga_iterations     <- 500
   ga_seed           <- 123
   random_replicates <- 100
   random_seed       <- 456
   ```

6. Run:

   ```r
   source("optimize_training_population.R")
   ```

7. Use `_selected_training.csv` and `_target_ids.csv` in a fixed-split GBLUP analysis. Inspect `_ga_replicates.csv`, `_ga_convergence.csv`, and `_random_benchmark.csv` to assess search consistency and improvement over random sets.

Optimization outputs:

| Suffix | Content |
|---|---|
| `_selected_training.csv` | Best training set |
| `_target_ids.csv` | Target set used for downstream fixed-split prediction |
| `_candidate_ids.csv` | Candidate population used by the optimizer |
| `_ga_replicates.csv` | Seed, objective, iterations recorded, and best-replicate flag |
| `_ga_training_sets.csv` | Membership of every GA training set |
| `_ga_convergence.csv` | Objective trajectory for each GA replicate |
| `_random_benchmark.csv` | Objectives for equal-size random training sets |
| `_random_training_sets.csv` | Membership of the random training sets |
| `_pc_variance.csv` | GRM eigenvalues, variance percentages, and PCs used |
| `_run_info.txt` | Design settings and result summary |

### Use case 4: Optimize without a separate target population

1. Prepare a GRM and, optionally, a candidate-ID CSV.
2. Configure `optimize_training_population.R`:

   ```r
   grm_file           <- "results/maize_grm_grm.csv"
   candidate_ids_file <- NULL
   target_ids_file    <- NULL
   optimization_mode  <- "no_test"
   training_size      <- 200
   criterion          <- "PEVMEAN"
   output_dir         <- "results"
   run_label          <- "maize_no_test_tp"
   ```

3. Configure the PC selection, GA, and random-benchmark settings as in Use case 3.
4. Run:

   ```r
   source("optimize_training_population.R")
   ```

5. The generated `_target_ids.csv` contains the complement of the best selected training set within the candidate population.

### Use case 5: Repeated GBLUP cross-validation with one record per genotype

1. Prepare a GRM and a phenotype table containing one row per genotype.
2. Open `gblup_cv_single.R` and configure the files and columns:

   ```r
   grm_file   <- "results/maize_grm_grm.csv"
   pheno_file <- "data/stage1_blues.csv"
   output_dir <- "results"
   run_label  <- "maize_single_cv"

   id_col     <- "GenoID"
   trait_col  <- "Trait"
   strata_col <- "Population"
   weight_col <- "Weight"
   se_col     <- "SE"
   ```

3. Configure the model. With stage-1 weights, the following fixes each residual variance to `1 / weight` and includes a non-marker effect:

   ```r
   fix_stage1_error  <- TRUE
   nonmarker_effects <- TRUE
   ```

4. Configure repeated cross-validation:

   ```r
   evaluation_mode <- "repeated_cv"
   k_folds         <- 5
   n_iterations    <- 50
   base_seed       <- 1000
   ```

5. Choose one stratification configuration from the Input formats section.
6. Run:

   ```r
   source("gblup_cv_single.R")
   ```

7. Read `_summary.csv` for the across-iteration results, `_metrics.csv` for fold-level results, and `_cv_predictions.csv` for out-of-sample genomic predictions.

### Use case 6: CV1 for untested genotypes across environments

1. Prepare a GRM and a phenotype table with one row per genotype × environment combination.
2. Open `gblup_cv_gxe.R` and configure the files and columns:

   ```r
   grm_file   <- "results/maize_grm_grm.csv"
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

3. Configure the model and CV1:

   ```r
   fix_stage1_error  <- TRUE
   nonmarker_effects <- TRUE
   evaluation_mode   <- "repeated_cv"
   cv_scheme         <- "CV1"
   k_folds           <- 5
   n_iterations      <- 50
   base_seed         <- 1000
   ```

4. Optionally provide named heritabilities for environment-specific accuracy; leave `env_h2 <- NULL` to use the estimated Cullis heritability:

   ```r
   env_h2 <- c(E1 = 0.62, E2 = 0.55, E3 = 0.58)
   ```

5. Choose a stratification configuration and run:

   ```r
   source("gblup_cv_gxe.R")
   ```

6. Use `_cv_predictions.csv` for genotype-level predictions and `_cv_predictions_by_env.csv` for out-of-sample predictions in each environment.

### Use case 7: CV2 for tested genotypes in unobserved environments

1. Prepare the same GRM and genotype × environment phenotype format used for CV1.
2. Configure `gblup_cv_gxe.R` as in Use case 6, changing the scheme:

   ```r
   evaluation_mode <- "repeated_cv"
   cv_scheme       <- "CV2"
   ```

3. Set `k_folds`, `n_iterations`, `base_seed`, model settings, and stratification.
4. Run:

   ```r
   source("gblup_cv_gxe.R")
   ```

5. Review `Untested_val_rows` in `_metrics.csv`. It counts held-out records whose genotype had no other observed environment in the training data and was therefore predicted as an untested genotype.

### Use case 8: Evaluate or predict a fixed target population

Use this workflow after STPGA optimization or with user-defined training and target ID files.

1. Prepare the GRM, phenotype table, training-ID CSV, and target-ID CSV. Use the same GRM that was used for training-population optimization. Every listed genotype must occur in the phenotype table; training genotypes require observed trait values, while target trait values may be observed or missing.
2. Configure either `gblup_cv_single.R` or `gblup_cv_gxe.R`:

   ```r
   evaluation_mode   <- "fixed_split"
   training_ids_file <- "results/maize_targeted_tp_selected_training.csv"
   target_ids_file   <- "results/maize_targeted_tp_target_ids.csv"
   ```

3. For G × E data, use genotype-level CV1:

   ```r
   cv_scheme <- "CV1"
   ```

4. Configure file paths, column names, weights, and model settings as in the corresponding repeated-CV tutorial.
5. Run the selected analysis:

   ```r
   source("gblup_cv_single.R")
   # or
   source("gblup_cv_gxe.R")
   ```

6. `_split_ids.csv` records every genotype's training or target role. Target predictions are out-of-sample. When target phenotypes are present, the metrics provide a retrospective evaluation; when they are missing, predictions are still produced but cannot contribute to evaluation metrics.

## GBLUP outputs

| Suffix | Content |
|---|---|
| `_metrics.csv` | Metrics for each iteration and fold; fold 0 contains predictions pooled within an iteration |
| `_summary.csv` | Mean and standard deviation of each metric across iterations |
| `_cv_predictions.csv` | Out-of-sample genomic prediction per genotype and iteration, plus the across-iteration mean |
| `_cv_predictions_by_env.csv` | G × E only: observed and out-of-sample predicted value for each genotype × environment row |
| `_train_predictions.csv` | In-sample training predictions; not a measure of predictive performance |
| `_train_se.csv` | Standard errors for the in-sample predictions |
| `_intercepts.csv` | Intercept from each fitted model |
| `_split_ids.csv` | Fixed-split only: training or target role for each included genotype |
| `_derived_strata.csv` | GRM PC coordinates and selected stratum when GRM-derived strata are enabled |
| `_derived_strata_selection.csv` | Candidate cluster counts and silhouette diagnostics |
| `_derived_strata_pc_variance.csv` | GRM eigenvalues, variance percentages, and PCs used for clustering |
| `_run_info.txt` | Analysis settings and package versions |

## Reported metrics

| Metric | Interpretation |
|---|---|
| `Vg`, `Vge`, `Vr`, `Vre` | REML variance-component estimates from the training data |
| `Residual_scale` | Estimated residual multiplier; equals 1 when stage-1 error is fixed |
| `Ve` | Mean stage-1 error variance |
| `PVE` | Proportion of variance attributed to the genomic main effect |
| `Heritability` | Cullis heritability based on the mean prediction-error variance of genotype differences |
| `Heritability_min`, `Heritability_max` | Cullis bounds based on the largest and smallest standard error of a difference |
| `GOF` | In-sample correlation between observed and fitted values |
| `PredictionAbility` | Correlation between observed values and out-of-sample predictions |
| `Accuracy` | Prediction ability divided by the square root of heritability |
| `RMSE` | Root mean squared out-of-sample prediction error |
| `Relative_RMSE` | RMSE as a percentage of the absolute trait mean |
| `Standardized_RMSE` | RMSE divided by the observed trait standard deviation; G × E removes environment means first |
| `Converged` | `1` when ASReml reports convergence and `0` otherwise |

For G × E analyses, prediction ability and GOF are calculated within environment and then averaged. Correlations and RMSE give each phenotype record equal weight; stage-1 weights affect model fitting, not metric calculation.
