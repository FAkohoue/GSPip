# GSPip

GSPip is an R package for genomic selection analyses. It integrates: (i)
genotype import and quality control, (ii) construction of genomic
relationship matrices (GRMs), (iii) optimisation of training
populations, and (iv) genomic best linear unbiased prediction (GBLUP)
with or without genotype-by-environment (G × E) effects.

Numeric dosage tables, HapMap files, and Variant Call Format (VCF) files
are supported. Genomic prediction can be fitted with ASReml-R by
residual maximum likelihood (REML) or with Bayesian Generalized Linear
Regression (BGLR) by Bayesian reproducing kernel Hilbert space (RKHS)
regression.

## Main functions

| Function | Purpose | Principal result |
|----|----|----|
| [`build_grm()`](https://FAkohoue.github.io/GSPip/reference/build_grm.md) | Import genotypes, apply sample and marker quality control, impute missing calls, and construct an additive or dominance GRM | `gspip_grm` object and auditable GRM output files |
| [`optimize_training_population()`](https://FAkohoue.github.io/GSPip/reference/optimize_training_population.md) | Select a fixed-size training population with Selection of Training Populations by Genetic Algorithm (STPGA) | `gspip_training_optimization` object, selected IDs, convergence records, and random benchmarks |
| [`run_gblup_single()`](https://FAkohoue.github.io/GSPip/reference/run_gblup_single.md) | Evaluate one stage-1 best linear unbiased estimate (BLUE) per genotype by repeated cross-validation or a fixed split | `gspip_gblup_single` object, metrics, and genotype predictions |
| [`run_gblup_gxe()`](https://FAkohoue.github.io/GSPip/reference/run_gblup_gxe.md) | Evaluate stage-1 BLUEs across environments under CV1, CV2, or a fixed genotype-level split | `gspip_gblup_gxe` object, metrics, genotype predictions, and environment-specific predictions |
| [`gspip_script()`](https://FAkohoue.github.io/GSPip/reference/gspip_script.md) | Locate or copy an editable workflow template | Path to the installed or copied template |

The [complete workflow
tutorial](https://FAkohoue.github.io/GSPip/docs/articles/GSPip-workflow.md)
provides step-by-step analyses for every supported use case.

## Installation

Install the development version from GitHub:

``` r
install.packages("pak")
pak::pak("FAkohoue/GSPip")
```

Load the package:

``` r
library(GSPip)
```

The core installation imports `AGHmatrix`, `ASRgenomics`, `cluster`, and
`data.table`. Additional packages are required for specific analyses:

| Package or software | Use                                              |
|---------------------|--------------------------------------------------|
| `STPGA`             | Training-population optimisation                 |
| `rrBLUP`            | Expectation-maximisation genotype imputation     |
| `BGLR`              | Open-source Bayesian genomic prediction          |
| ASReml-R            | Licensed REML genomic prediction                 |
| Java and Beagle     | Linkage-disequilibrium-aware genotype imputation |

Install the optional R packages required for STPGA,
expectation-maximisation imputation, and BGLR:

``` r
install.packages(c("STPGA", "rrBLUP", "BGLR"))
```

ASReml-R requires a separate installation and a valid licence. Beagle
requires Java and a local Beagle JAR file.

## Genotype input

[`build_grm()`](https://FAkohoue.github.io/GSPip/reference/build_grm.md)
accepts three genotype formats.

| Format | Required structure | Internal dosage |
|----|----|----|
| Numeric table | Markers in rows, a marker-ID column, and at least two genotype columns | `0`, `1`, `2`, or missing |
| HapMap | Eleven standard HapMap columns through `QCcode`, followed by sample columns | Count of the second allele in `alleles` |
| VCF or VCF.GZ | Nine fixed VCF columns through `FORMAT`, followed by sample columns with diploid `GT` | `0/0 = 0`, heterozygote `= 1`, and `1/1 = 2` |

With `input_format = "auto"`, the extensions `.vcf`, `.vcf.gz`, `.hmp`,
`.hmp.txt`, `.hapmap`, and `.hapmap.txt` are detected automatically.
Other extensions are treated as numeric input.

Only diploid biallelic markers are converted to dosage. Unsupported
records can be removed or rejected with `multiallelic_action = "remove"`
or `multiallelic_action = "error"`. VCF records without an ID receive
the identifier `CHROM:POS:REF:ALT`.

### Numeric dosage table

``` text
SNP       CHR  POS    REF  ALT  G001  G002  G003
S1_10234  1    10234  T    G    0     2     2
S1_18945  1    18945  G    T    0     0     0
S1_27656  1    27656  C    T    0     NA    2
```

Every column not assigned to marker identification, coordinates,
alleles, metadata, or marker weights is interpreted as a genotype.
Marker IDs and genotype names must be unique.

### Marker weights

Marker weights may be supplied as: (i) an inline column in numeric input
or (ii) an external table for any supported genotype format.

``` text
MarkerID  Weight
S1_10234  1.25
S1_18945  0.80
```

Weights must be finite and non-negative. At least one retained weight
must be positive. Marker weights are supported for the additive VanRaden
GRM, but not for the dominance Su GRM.

### Quality control and imputation

Sample filtering precedes marker filtering. Marker call rate, minor
allele frequency, heterozygosity, and polymorphism are calculated from
observed calls before imputation.

| `imputation_method` | Method |
|----|----|
| `"none"` | Retain missing calls and pass them to [`AGHmatrix::Gmatrix()`](https://rdrr.io/pkg/AGHmatrix/man/Gmatrix.html) |
| `"wright"` | Random imputation from marker allele frequency and genotype-specific Wright inbreeding estimates |
| `"random_hwe"` | Random imputation from Hardy-Weinberg genotype frequencies |
| `"rrblup_em"` | Multivariate expectation-maximisation imputation with [`rrBLUP::A.mat()`](https://rdrr.io/pkg/rrBLUP/man/A.mat.html) |
| `"beagle"` | External Beagle imputation using chromosome, position, reference allele, and alternative allele |

## Phenotype input

### One BLUE per genotype

[`run_gblup_single()`](https://FAkohoue.github.io/GSPip/reference/run_gblup_single.md)
requires one row per genotype.

``` text
GenoID  Trait  Population  SE     Weight
G001    10.66  P1          0.566  3.122
G002    10.92  P1          0.385  6.747
G003    10.38  P2          0.457  NA
G004    NA     P2          NA     NA
```

### One BLUE per genotype and environment

[`run_gblup_gxe()`](https://FAkohoue.github.io/GSPip/reference/run_gblup_gxe.md)
requires one row per genotype × environment combination.

``` text
GenoID  Env  Trait  Population  SE     Weight
G001    E1   10.02  P1          0.400  6.250
G001    E2   11.68  P1          0.459  4.747
G002    E1   11.01  P1          0.448  4.982
G002    E2   NA     P1          NA     NA
```

For an observed BLUE, a positive finite value in `Weight` is used first.
Otherwise, a positive finite standard error (SE) is converted to
`1 / SE^2`. Rows with a missing trait remain available for prediction
and require neither value.

When both precision columns are absent, GSPip uses equal weights,
estimates the residual scale, and excludes non-marker effects.
Therefore, stage-1 precision information should be supplied when
available.

Every phenotype genotype must occur in the GRM. Genotypes present only
in the GRM receive missing values in prediction tables unless they are
included through a fixed target file.

## Genomic relationship matrix and ID files

Optimisation and GBLUP functions require a square, symmetric GRM stored
as a comma-separated value (CSV) file. The first column contains
genotype IDs; the remaining column names contain the same IDs.

``` text
GenoID  G001   G002   G003
G001    1.85  -0.34  -0.21
G002   -0.34   1.91  -0.28
G003   -0.21  -0.28   1.76
```

Candidate, training, and target files contain one unique ID per row in
the column specified by `id_col`.

``` text
GenoID
G001
G002
G003
```

All listed IDs must occur in the GRM. Candidate and target sets must be
disjoint. Fixed training and target sets must also be disjoint.

## Quick example

GSPip contains aligned example files with 120 genotypes, 240 markers,
three populations, and three environments.

Locate the numeric genotype file and create an output directory:

``` r
library(GSPip)

geno_file <- system.file(
  "extdata", "example_genotypes_numeric.csv",
  package = "GSPip"
)

output_dir <- file.path(tempdir(), "gspip-example")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
```

Construct an additive VanRaden GRM:

``` r
grm_fit <- build_grm(
  geno_file = geno_file,
  output_dir = output_dir,
  input_format = "numeric",
  run_label = "example_grm",
  grm_method = "VanRaden",
  qc_call_rate = 0.95,
  qc_maf = 0.05,
  imputation_method = "wright",
  qc_seed = 123
)

grm_fit
grm_fit$output_files
```

The returned object contains the GRM, retained marker matrix, marker
map, quality-control tables, optional clusters, and paths to all written
files.

## Training-population optimisation

[`optimize_training_population()`](https://FAkohoue.github.io/GSPip/reference/optimize_training_population.md)
uses GRM principal components as genomic predictors. Trait values are
not used during selection.

Two optimisation designs are available:

| Mode | Objective population |
|----|----|
| `"targeted"` | A supplied set of target genotypes |
| `"no_test"` | The complement of each proposed training set within the candidate population |

The available STPGA criteria are `PEVMEAN`, `PEVMEAN2`, `PEVMAX`,
`CDMEAN`, `CDMEAN2`, and `CDMAX`. All are minimised. Replicated
genetic-algorithm searches are reported with equal-size random
benchmarks.

## Genomic prediction engines

Both GBLUP functions require an explicit engine.

| Engine | Model fitting | Supported random terms |
|----|----|----|
| `engine = "asreml"` | REML with ASReml-R | Genomic effects and optional non-marker effects |
| `engine = "bglr"` | Bayesian RKHS regression with BGLR | Genomic main effect and, for G × E, a genomic interaction kernel |

BGLR estimates the residual scale and requires:

``` r
fix_stage1_error  <- FALSE
nonmarker_effects <- FALSE
```

Its posterior quantities use the suffix `_posterior`. They are not REML
or Cullis estimates. The Markov chain Monte Carlo (MCMC) controls are
`bglr_n_iter`, `bglr_burn_in`, `bglr_thin`, and `bglr_seed`.

## Evaluation designs

| Analysis | Fold unit | Prediction target |
|----|----|----|
| Single-record cross-validation | Genotype | One held-out genotype per fold assignment |
| G × E CV1 | Genotype | Untested genotypes in all environments |
| G × E CV2 | Genotype × environment record | Tested genotypes in environments where their record is held out |
| Fixed split | Genotype | A predefined target population |

Folds may be balanced with: (i) a supplied population or stratum column,
(ii) genomic strata derived from GRM principal components, or (iii) no
genotype stratification.

Fixed-split G × E analysis uses CV1 because training and target
membership is defined at genotype level.

## Output files

Every analysis writes files with the prefix `<run_label>_`. Existing
outputs are preserved unless `overwrite_outputs = TRUE`.

### GRM construction

| Suffix | Content |
|----|----|
| `_grm.csv` | Genomic relationship matrix |
| `_sample_qc.csv` | Sample missingness and filter status |
| `_marker_qc.csv` | Marker statistics and filter status |
| `_marker_map.csv` | Retained marker coordinates and alleles |
| `_qc_report.txt` | Quality-control counts and thresholds |
| `_grm_info.txt` | Matrix diagnostics and inversion result |
| `_clusters.csv` | Optional GRM principal components and cluster assignments |
| `_imputed_genotypes.csv` | Optional final marker matrix |

### Training-population optimisation

| Suffix | Content |
|----|----|
| `_selected_training.csv` | Selected training genotypes |
| `_target_ids.csv` | Target genotypes for fixed-split prediction |
| `_ga_replicates.csv` | Objective value and seed for every genetic-algorithm replicate |
| `_ga_convergence.csv` | Objective trajectory for every replicate |
| `_random_benchmark.csv` | Objectives for equal-size random training sets |
| `_pc_variance.csv` | GRM eigenvalues and principal components used |

### GBLUP evaluation

| Suffix | Content |
|----|----|
| `_metrics.csv` | Fold- and iteration-level statistics |
| `_summary.csv` | Mean and standard deviation across iterations |
| `_cv_predictions.csv` | Out-of-sample genotype predictions |
| `_cv_predictions_by_env.csv` | G × E predictions by genotype and environment |
| `_train_predictions.csv` | In-sample genomic predictions |
| `_train_se.csv` | ASReml standard errors or BGLR posterior standard deviations |
| `_split_ids.csv` | Training and target membership for a fixed split |
| `_derived_strata*.csv` | Optional genomic-strata assignments and diagnostics |
| `_run_info.txt` | Analysis settings and package versions |

The principal performance statistics are prediction ability, accuracy,
root mean squared error (RMSE), standardised RMSE, proportion of
variance explained (PVE), and heritability. Environment-specific and
stratum-specific statistics are added when their corresponding analyses
are requested.

## Example data

The complete example set is available as R objects:

``` r
data(gs_example_all)
names(gs_example_all)
dim(gs_example_all$geno)  # 120 × 240
dim(gs_example_all$grm)   # 120 × 120
```

Equivalent numeric, HapMap, and VCF genotype files are installed under
`extdata`. The same directory contains the GRM, single-record BLUEs, G ×
E BLUEs, marker weights, candidate IDs, and target IDs used in the
vignette.

## Editable workflow templates

The exported functions form the package application programming
interface (API). Editable templates are also available for analyses
organised through settings files:

``` r
gspip_script("grm_builder", copy_to = "analysis/")
gspip_script("optimize_training_population", copy_to = "analysis/")
gspip_script("gblup_cv_single", copy_to = "analysis/")
gspip_script("gblup_cv_gxe", copy_to = "analysis/")
```

Each template is a thin wrapper around the corresponding exported
function. Therefore, direct function calls and template-based analyses
use the same implementation.

## Documentation

- [Complete workflow
  tutorial](https://FAkohoue.github.io/GSPip/docs/articles/GSPip-workflow.md)
- [Function reference](https://FAkohoue.github.io/GSPip/reference/)
- [Package website](https://FAkohoue.github.io/GSPip/)
