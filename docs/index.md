# GSPip

An R package for marker quality control, construction of genomic
relationship matrices (GRMs) with optional marker weights, STPGA
training-population optimization, and GBLUP evaluation with or without
genotype-by-environment effects.

The programmatic package API provides
[`build_grm()`](https://FAkohoue.github.io/GSPip/reference/build_grm.md),
[`optimize_training_population()`](https://FAkohoue.github.io/GSPip/reference/optimize_training_population.md),
[`run_gblup_single()`](https://FAkohoue.github.io/GSPip/reference/run_gblup_single.md),
and
[`run_gblup_gxe()`](https://FAkohoue.github.io/GSPip/reference/run_gblup_gxe.md).
Editable thin-wrapper workflows are installed as templates and can be
located or copied with
[`gspip_script()`](https://FAkohoue.github.io/GSPip/reference/gspip_script.md).

The scripts are designed to fail early on invalid settings or inputs.
They check identifiers, matrix dimensions and positive definiteness,
finite trait values and weights, fold settings, output collisions, and
the conditions needed by each model. Existing outputs are never replaced
unless `overwrite_outputs <- TRUE`.

| Order | Script | Purpose |
|----|----|----|
| 1 | `grm_builder.R` | Marker QC, optional imputation, equal- or supplied-marker-weight GRM, optional clusters |
| 2 | `optimize_training_population.R` | Select a fixed-size training set with STPGA and compare it with equal-size random sets |
| 3a | `gblup_cv_single.R` | Repeated CV or fixed training/target evaluation with one BLUE per genotype |
| 3b | `gblup_cv_gxe.R` | Repeated CV or fixed training/target evaluation with one BLUE per genotype x environment |

To run a script, edit the `SETTINGS` section at the top, then run the
whole script. Every output file name starts with the `run_label`
setting.

## Installation and quick start

Install the development version from GitHub:

``` r
install.packages("pak")
pak::pak("FAkohoue/GSPip")
```

Install from the local checkout during development:

``` r
install.packages("D:/OneDrive - CGIAR/GitHub libraries/GSPip",
                 repos = NULL, type = "source")
```

Build a GRM directly from R:

``` r
library(GSPip)

result <- build_grm(
  geno_file = "path/to/genotypes.csv",
  output_dir = "path/to/results",
  input_format = "auto",
  run_label = "trial1",
  imputation_method = "none"
)
print(result)
```

Optimize a training population programmatically with
[`optimize_training_population()`](https://FAkohoue.github.io/GSPip/reference/optimize_training_population.md),
or copy an editable workflow template:

Copy an editable workflow template:

``` r
gspip_script("optimize_training_population", copy_to = "analysis/")
gspip_script("gblup_cv_single", copy_to = "analysis/")
```

ASReml-R is required only when
[`run_gblup_single()`](https://FAkohoue.github.io/GSPip/reference/run_gblup_single.md),
[`run_gblup_gxe()`](https://FAkohoue.github.io/GSPip/reference/run_gblup_gxe.md),
or their editable wrappers are executed. It still needs a separate VSNi
installation and valid license. Use a new `run_label` for each analysis;
the workflows stop before replacing existing results unless explicitly
allowed.

Paths may be absolute or relative to the R working directory. For a
reproducible shared analysis, prefer project-relative paths and record
the commit hash used for the run.

## Complete example data

GSPip ships a complete, mutually aligned example set with 120 genotypes,
240 markers, three populations, and three environments. Load every
R-native object as one list:

``` r
data(gs_example_all)
names(gs_example_all)
dim(gs_example_all$geno) # 120 x 240
dim(gs_example_all$grm)  # 120 x 120
```

Individual datasets are also available through
[`data()`](https://rdrr.io/r/utils/data.html), including
`gs_example_geno`, `gs_example_marker_map`, `gs_example_marker_weights`,
`gs_example_blues_single`, `gs_example_blues_gxe`,
`gs_example_candidate_ids`, `gs_example_target_ids`, `gs_example_grm`,
and `gs_example_truth`.

Equivalent numeric, HapMap, and VCF genotype files and all workflow
inputs are installed under `extdata`:

``` r
numeric_file <- system.file(
  "extdata", "example_genotypes_numeric.csv", package = "GSPip"
)
grm_file <- system.file("extdata", "example_grm.csv", package = "GSPip")
single_blues <- system.file(
  "extdata", "example_BLUEs_single.csv", package = "GSPip"
)
gxe_blues <- system.file(
  "extdata", "example_BLUEs_gxe.csv", package = "GSPip"
)
```

The set also includes external marker weights, candidate and target ID
tables, simulation truth, and `example_manifest.csv`, which maps every
file to its purpose. Regenerate all packaged data from the repository
root with:

``` r
source("data-raw/generate_example_data.R")
```

The generator uses base R, protects existing outputs by default, and
accepts `GS_EXAMPLE_OVERWRITE=true` for deliberate non-interactive
regeneration.

> \[!IMPORTANT\] **Check for a weight or SE column before running
> cross-validation.** The two CV scripts do not behave the same way when
> both columns are absent. Do not create an artificial `Weight = 1`
> column while keeping `fix_stage1_error <- TRUE`; that would force the
> residual variance to 1 rather than estimate it from the data.

### If the data has no Weight or SE column

| Script | Can it run without Weight and SE? | Required action |
|----|----|----|
| `gblup_cv_single.R` | Yes | Set `weight_col <- NULL`, `fix_stage1_error <- FALSE`, and `nonmarker_effects <- FALSE` |
| `gblup_cv_gxe.R` | No | Supply a valid weight column or SE column. The script stops before model fitting if both are absent |

For the single-genotype script, the no-weight configuration is:

``` r
weight_col        <- NULL
fix_stage1_error  <- FALSE
nonmarker_effects <- FALSE
```

This gives every observed BLUE equal weight and lets ASReml estimate the
residual variance. The `se_col` setting is ignored. The non-marker
effect must be disabled because, with one record per genotype and an
estimated residual scale, it cannot be separated from the residual.

For the G × E script, at least one real precision measure is required:

- If the weight column exists, its positive finite values are used.
- If a weight is missing but the SE is available, the script uses
  `1 / SE^2` for that row.
- The entire weight column may be absent when the SE column is present,
  and vice versa.
- A row without a trait value needs neither a weight nor an SE because
  it does not contribute to model fitting or metrics.
- If an observed row has neither a valid weight nor a valid SE, the
  script stops and identifies the problem before fitting any model.

### If the data has no Population or stratum column

Both CV scripts support three explicit choices. Do not leave
`strata_col <- "Population"` when that column is absent—the input check
will stop and identify the missing column.

#### Choice 1: run without stratification

``` r
strata_col             <- NULL
derive_strata_from_grm <- FALSE
```

In `gblup_cv_single.R` and G × E CV1, genotypes are allocated across
folds without population balancing. In G × E CV2, rows remain balanced
within environments, but not within populations. Population-specific
`PA_CV_within_*` metrics are not produced.

#### Choice 2: derive genomic strata from the GRM

``` r
strata_col              <- NULL
derive_strata_from_grm  <- TRUE
grm_strata_n_clusters   <- NULL # Automatically select k
grm_strata_max_clusters <- 10   # Evaluate k = 2, ..., 10 when the data size allows
grm_strata_n_pcs        <- 2    # May be increased, for example to 5 or 10
grm_strata_seed         <- 123
```

The user does not need to know the number of clusters in advance. With
`grm_strata_n_clusters <- NULL`, the script:

1.  Calculates the requested number of leading PCs from the full GRM.
2.  Fits k-means solutions for each candidate from 2 through
    `grm_strata_max_clusters`, subject to the number of genotypes and
    distinct PC positions.
3.  Calculates the average silhouette width for every candidate
    solution.
4.  Selects the candidate with the largest average silhouette width.
5.  Uses the selected clusters to balance the CV folds.

Average silhouette width compares within-cluster cohesion with
separation from neighboring clusters; larger values indicate a
better-separated candidate partition. The method is used here to select
an appropriate candidate count, not to claim a known true number of
biological populations.

`grm_strata_n_pcs <- 2` is only the default, not a restriction. Two PCs
are easy to visualize and reduce the influence of low-variance noise,
but they can miss structure expressed on later axes. Set this to a
larger positive integer when the eigenvalue profile or biological
knowledge supports doing so. K-means uses all requested PCs with their
eigenvalue-based scale, so later, lower-variance PCs contribute less
distance than leading PCs. Adding many weak PCs can dilute real
separation, so compare the diagnostics rather than assuming that more
PCs are always better.

The result is reproducible for fixed inputs, settings, R version, and
package versions. The script writes three diagnostic files:

| File | Content |
|----|----|
| `<run_label>_derived_strata_selection.csv` | Every candidate `K`, its average silhouette, minimum cluster size, total within-cluster sum of squares, and the selected solution |
| `<run_label>_derived_strata.csv` | `GenoID`, every PC used (`PC1` through `PCn`), and the selected `Stratum` for every genotype in the GRM |
| `<run_label>_derived_strata_pc_variance.csv` | Eigenvalue, individual and cumulative variance percentage, and whether each PC was used for clustering |

Inspect all three files before interpreting the CV results. If prior
biological knowledge or the diagnostics justify a fixed value, set
`grm_strata_n_clusters` to that integer; this is an explicit manual
override and only that value is fitted.

After inspecting the automatic diagnostics, the selected count can also
be supplied as `n_clusters` in `grm_builder.R` to generate a reusable
`_clusters.csv`; its `Population` column can then be merged into the
phenotype table before using Choice 3.

Method reference: P. J. Rousseeuw (1987), [“Silhouettes: A graphical aid
to the interpretation and validation of cluster
analysis”](https://doi.org/10.1016/0377-0427(87)90125-7).

> \[!CAUTION\] GRM-derived clusters are unsupervised genomic groups, not
> automatically known biological populations or named ancestry groups.
> Automatic silhouette selection is a data-driven aid, not proof of a
> true population count. Because clustering uses markers but not trait
> values, it does not introduce phenotype leakage; however, balancing
> these clusters across folds evaluates prediction within the
> represented genomic structure. It does not evaluate prediction into a
> completely unseen population.

#### Choice 3: use a supplied column

``` r
strata_col             <- "Population"
derive_strata_from_grm <- FALSE
```

The named column must exist, contain no missing or empty values, and
assign each genotype to only one stratum. A supplied, scientifically
justified population assignment is preferable to GRM-derived clusters
when available.

`strata_col` and GRM derivation are mutually exclusive. If both are
enabled, the script stops before reading the data and asks the user to
select one source.

## Training-population optimization with STPGA

`optimize_training_population.R` selects a requested number of candidate
genotypes before phenotyping or model fitting. It uses genomic
principal-component predictors derived from the GRM; no trait values are
read by the optimizer. This separates training-set design from
downstream assessment and avoids selecting a set because it happened to
predict the same observed phenotypes well.

The optimizer has two modes:

| Mode | Required ID files | Objective population |
|----|----|----|
| `optimization_mode <- "targeted"` | Target file; candidate file optional | The supplied target genotypes. This is the recommended mode when the future prediction population is already known |
| `optimization_mode <- "no_test"` | Candidate file optional | For each proposed training set, its complement in the candidate population is the test set |

Candidate and target files are CSV files with one unique, non-empty ID
per row in the column named by `id_col` (default `GenoID`). Every ID
must be in the GRM. The two sets must be disjoint. With
`candidate_ids_file <- NULL`, targeted mode uses every GRM genotype not
in the target file as a candidate; no-test mode uses every GRM genotype
as a candidate.

Minimal targeted configuration:

``` r
grm_file           <- "results/example_grm.csv"
candidate_ids_file <- "data/candidate_ids.csv" # Or NULL
target_ids_file    <- "data/target_ids.csv"
optimization_mode  <- "targeted"
training_size      <- 200
criterion          <- "PEVMEAN"
```

The available criteria are `PEVMEAN`, `PEVMEAN2`, `PEVMAX`, `CDMEAN`,
`CDMEAN2`, and `CDMAX`. They are STPGA design criteria and are all
minimized by the algorithm. PEV criteria summarize prediction-error
variance, while CD criteria summarize coefficient-of-determination
design loss. They are not the `PVE` variance-component metric reported
by the GBLUP scripts.

### PCs and genetic-algorithm replication

The default `stpga_n_pcs <- NULL` automatically uses enough
positive-eigenvalue PCs to reach `stpga_pc_variance_target`, subject to
`stpga_max_pcs`. Set `stpga_n_pcs` to a positive integer for a manual
override. The optimizer writes the complete eigenvalue profile and
indicates exactly which PCs were used.

STPGA’s search is stochastic and does not guarantee the global optimum.
The script therefore runs `ga_replicates` independent searches, records
each seed and solution, and selects the replicate with the smallest
recomputed objective. It also evaluates `random_replicates`
independently sampled training sets of exactly the same size. A selected
set should not be treated as useful merely because the genetic algorithm
returned it: inspect convergence, consistency among replicates, and its
improvement over the random objective distribution.

### Optimizer outputs

| File | Content |
|----|----|
| `_selected_training.csv` | Best training set; directly accepted by both GBLUP scripts |
| `_target_ids.csv` | Supplied target set in targeted mode, or the best training set’s candidate complement in no-test mode; directly accepted by both GBLUP scripts |
| `_candidate_ids.csv` | Candidate population actually used |
| `_ga_replicates.csv` | Seed, objective, recorded iterations, and best-replicate flag |
| `_ga_training_sets.csv` | Every genotype selected in every GA replicate |
| `_ga_convergence.csv` | Objective trajectory returned by each GA replicate |
| `_random_benchmark.csv`, `_random_training_sets.csv` | Equal-size random objective values and membership |
| `_pc_variance.csv` | GRM eigenvalues, variance percentages, and PCs used |
| `_run_info.txt` | Design settings, best and random results, paths, and package versions |

### Evaluate the selected set in GBLUP

Set the following in `gblup_cv_single.R` or `gblup_cv_gxe.R`:

``` r
evaluation_mode   <- "fixed_split"
training_ids_file <- "results/STPGA_training_selected_training.csv"
target_ids_file   <- "results/STPGA_training_target_ids.csv"
```

Fixed-split mode performs one fit: every selected training genotype is
assigned to training and every target genotype is held out together.
`k_folds`, `n_iterations`, and `base_seed` do not define this split. The
phenotype file must contain all listed IDs; training genotypes must have
observed trait values. Target trait values may be present for
retrospective assessment or missing for prediction only. The scripts
write `_split_ids.csv` so the exact roles are preserved with the
results.

For `gblup_cv_gxe.R`, fixed-split evaluation is genotype-level CV1 only:
all environment records of a target genotype are held out. Set
`cv_scheme <- "CV1"`. The script deliberately stops if fixed-split mode
is combined with CV2, because a genotype-level STPGA selection does not
define which individual genotype-by-environment cells should be
observed.

The same GRM must be used by the optimizer and the GBLUP analysis. For
an honest performance comparison, define target membership without using
its trait values and compare the optimized set with multiple equal-size
random training sets under the same GBLUP model. Optimization improves a
genomic design criterion; it does not guarantee better prediction
ability for every trait or environment.

Package reference: [`STPGA` 5.2.1
documentation](https://cran.r-universe.dev/STPGA/doc/manual.html). The
package defines the targeted and no-test genetic algorithms, their
minimized criteria, and recommends replicated searches because the
genetic algorithm is not guaranteed to find a global optimum.

| Script | Packages |
|----|----|
| `grm_builder.R` | `data.table`, `AGHmatrix`, `ASRgenomics`; `rrBLUP` only for `rrblup_em` imputation |
| `optimize_training_population.R` | `data.table`, `STPGA` |
| CV scripts | `data.table`, `ASRgenomics`, `cluster`, and `asreml` with a valid VSNi license |

The repository’s GitHub Actions workflow parses every `.R` file and runs
synthetic numeric, HapMap, and VCF reader tests on each push and pull
request. This catches syntax errors and genotype-format regressions
without requiring ASReml or private data. A full numerical run still
requires the remaining packages, license, and representative input data.

## Reproducibility and safe sharing

- `base_seed`, `qc_seed`, `cluster_seed`, `grm_strata_seed`, `ga_seed`,
  and `random_seed` make fold assignment, imputation, clustering,
  optimization, and random benchmarks reproducible for fixed inputs and
  package versions.
- Each run-info file records the settings, normalized input paths, input
  counts, and
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).
- `data/` and `results/` are ignored by Git to reduce the risk of
  publishing private inputs or large generated files. Review
  `git status` before every push.
- Cross-validation metrics can legitimately be `NA` when a validation
  subset has too few usable observations, no variation, a zero trait
  mean, or undefined Cullis heritability. The summary’s `N_iter` reports
  how many finite iteration-level values contributed.
- Non-converged ASReml fits are retained and flagged by `Converged = 0`;
  review them before interpreting the summary.

Before making the repository public, add the license and citation
information appropriate for your organization and confirm that your
ASReml license permits the intended use. ASReml itself is not
redistributed here.

### Pre-release smoke test

Run this once on a small, representative copy of your real data in the
same R and ASReml environment your colleagues will use:

1.  Set `n_iterations <- 2`, `k_folds <- 2`, and use a temporary
    `run_label` for repeated CV.
2.  Run the GRM script and confirm that `_grm_info.txt` reports an
    invertible matrix, a positive smallest eigenvalue, and plausible
    diagonal/off-diagonal means.
3.  If using optimization, first use a small `training_size`,
    `ga_iterations`, `ga_replicates`, and `random_replicates`. Confirm
    the selected IDs are candidates, are unique, and are disjoint from
    the targets. Then restore the production settings.
4.  Run the relevant GBLUP script and confirm that every expected output
    is present and `Converged` is 1 for every fit. For repeated CV,
    `N_iter` should be 2 for the principal metrics; fixed-split output
    has one iteration.
5.  Check that prediction row counts and genotype/environment IDs match
    the inputs. In fixed-split mode, verify `_split_ids.csv` against the
    optimizer files.
6.  Run the same script again without changing `run_label`; it should
    stop before overwriting any result.
7.  Restore the intended fold, iteration, optimizer, and run-label
    settings before the production run.

## GRM script

### Accepted genotype formats

Format handling is implemented by this project’s internal [`read_geno()`
module](https://FAkohoue.github.io/GSPip/R/read_geno.R). It is not
imported from HapBlockR and is not exported as a public package
function.
[`build_grm()`](https://FAkohoue.github.io/GSPip/reference/build_grm.md)
and the editable GRM template therefore use one implementation rather
than maintaining separate parsers.

Set `input_format` to `"numeric"`, `"hapmap"`, or `"vcf"`. With
`"auto"`, `.vcf` and `.vcf.gz` are recognized as VCF, `.hmp`,
`.hmp.txt`, `.hapmap`, and `.hapmap.txt` as HapMap; every other
extension is treated as a numeric table. Set the format explicitly when
a file has a nonstandard extension.

| Format | Required layout | Internal dosage |
|----|----|----|
| Numeric table | Markers in rows; `marker_col`; at least two genotype columns containing 0, 1, 2, or missing | Already 0/1/2 |
| HapMap | The 11 standard HapMap columns through `QCcode`, followed by samples | Count of the second allele in `alleles`; common letter, two-letter, IUPAC, and 0/1/2 calls are accepted |
| VCF or VCF.GZ | Nine fixed VCF columns through `FORMAT`, followed by samples; a diploid `GT` field | Count of ALT: `0/0 = 0`, heterozygote = 1, `1/1 = 2` |

Only biallelic markers are used. `multiallelic_action <- "remove"`
reports and removes unsupported records; `"error"` stops instead.
Haploid and polyploid VCF genotypes are not silently recoded. Missing
numeric calls may be `NA`, blank, or `.`; missing HapMap/VCF calls use
their normal format-specific notation.

Marker IDs must be unique. If a VCF record has no `ID`, the script
creates `CHROM:POS:REF:ALT`; use that generated value in an external
marker-weight table.

For numeric input, every column not named by `marker_col`, the
coordinate/allele settings, `meta_cols`, or `marker_weight_col` is
treated as a genotype. List every non-genotype column. The following is
a valid numeric-table example; `CHR`, `POS`, `REF`, and `ALT` are
optional unless Beagle is selected.

| SNP      | CHR | POS   | REF | ALT | G001 | G002 | G003 | G004 | G005 | G006 | …   |
|----------|-----|-------|-----|-----|------|------|------|------|------|------|-----|
| S1_10234 | 1   | 10234 | T   | G   | 0    | 2    | 2    | 2    | 0    | 0    | …   |
| S1_18945 | 1   | 18945 | G   | T   | 0    | 0    | 0    | 0    | 2    | 0    | …   |
| S1_27656 | 1   | 27656 | C   | T   | 0    | NA   | 2    | 2    | 0    | 2    | …   |
| S1_36367 | 1   | 36367 | A   | C   | 0    | 0    | 0    | 2    | 0    | 0    | …   |
| S1_45078 | 1   | 45078 | G   | T   | 2    | 2    | 2    | 2    | 2    | 2    | …   |
| S1_53789 | 1   | 53789 | A   | T   | 0    | 0    | 2    | 2    | 2    | 2    | …   |
| S1_62500 | 1   | 62500 | G   | T   | 0    | 0    | 0    | 0    | NA   | 0    | …   |
| S1_71211 | 1   | 71211 | T   | A   | 0    | 0    | 2    | 0    | 0    | 2    | …   |
| S1_79922 | 1   | 79922 | A   | G   | 2    | 2    | 2    | 2    | 2    | 2    | …   |
| S1_88633 | 1   | 88633 | A   | C   | 2    | 2    | 2    | 2    | 2    | 2    | …   |

### Marker weights

Leave both weight settings `NULL` for AGHmatrix’s equal marker weights.
Numeric input may name an inline column:

``` r
marker_weight_col <- "MarkerWeight"
```

All three formats may instead use a separate table:

``` r
marker_weights_file      <- "marker_weights.csv"
marker_weights_id_col    <- "MarkerID"
marker_weights_value_col <- "Weight"
```

The external table must contain one finite, non-negative weight for
every biallelic marker accepted from the genotype file before QC; extra
marker IDs are ignored. The script aligns the weights by marker ID after
QC and passes the retained vector unchanged to
`AGHmatrix::Gmatrix(weights = ...)`. At least one retained weight must
be positive. Use weights only with `grm_method <- "VanRaden"`; the
script stops on weights with `"Su"`. Weight scale is not normalized. See
the official [`AGHmatrix::Gmatrix()`
documentation](https://cran.r-universe.dev/AGHmatrix/doc/manual.html).

### QC and imputation

QC is calculated internally from observed calls, before imputation:

1.  Remove genotypes whose missing fraction exceeds `qc_sweep_sample`.
2.  For the remaining genotypes, calculate every marker’s call rate,
    ALT/second-allele frequency, MAF, and observed heterozygosity.
3.  Apply `qc_call_rate`, `qc_maf`, and `qc_heterozygosity`, and remove
    markers with no calls or no allele variation.
4.  Apply the selected `imputation_method` to the retained matrix.

The two CSV QC reports retain both passing and failing rows, their
statistics, and a failure reason.

| `imputation_method` | What happens |
|----|----|
| `"none"` | Missing calls remain `NA`. **The script does not stop because calls are missing.** It calls `AGHmatrix::Gmatrix(..., impute.method = "none")`, so AGHmatrix receives the missing calls directly. |
| `"wright"` | Random 0/1/2 draws use each marker’s allele frequency and a Wright inbreeding estimate for each genotype. Reproducible with `qc_seed`. |
| `"random_hwe"` | Random 0/1/2 draws follow Hardy-Weinberg frequencies at each marker. Reproducible with `qc_seed`. |
| `"rrblup_em"` | [`rrBLUP::A.mat()`](https://rdrr.io/pkg/rrBLUP/man/A.mat.html) performs multivariate EM imputation. It can return fractional dosages and is therefore limited here to the additive VanRaden GRM. |
| `"beagle"` | The retained data are written to a temporary VCF, Beagle is run through Java, and the imputed calls at the original markers are read back. |

`none` means no imputation by this script or by the `Gmatrix()`
argument; it does not mean that the script rejects `NA`. Check
`_grm_info.txt` to see exactly how many calls were passed to AGHmatrix.

### Beagle requirements and safeguards

Beagle is an external Java program, not an R package. Download its JAR
separately, install Java, and set `beagle_jar`. Beagle imputation is LD-
and position-aware, so every retained marker must have a real
chromosome, positive integer position, REF allele, and one ALT allele:

- VCF already supplies these fields.
- Standard HapMap supplies chromosome, position, and the two alleles.
- Numeric input must provide the columns named by `chrom_col`,
  `position_col`, `ref_col`, and `alt_col`.

The script refuses to invent coordinates or alleles. It accepts numeric,
HapMap, and VCF input because each is converted to the same temporary
VCF after QC; users do not need to convert the file manually. Configure
memory and threads with `beagle_memory` and `beagle_threads`. Optional
`beagle_ref`, `beagle_map`, and `beagle_ne` are passed to Beagle. With
no map, Beagle assumes 1 cM per Mb. Its default effective population
size was designed for a large outbred population and may be
inappropriate for a breeding population, so set `beagle_ne` from
defensible population knowledge. If a reference panel introduces
additional variants, this script deliberately keeps only the original
post-QC markers for the GRM. The command and Beagle output are saved in
`_beagle.log`. See the official [Beagle
page](https://faculty.washington.edu/browning/beagle/beagle.html) and
[Beagle 5.5
documentation](https://faculty.washington.edu/browning/beagle/beagle_5.5_17Dec24.pdf).

### Remaining GRM steps

1.  [`AGHmatrix::Gmatrix()`](https://rdrr.io/pkg/AGHmatrix/man/Gmatrix.html)
    computes the `"VanRaden"` additive or `"Su"` dominance GRM. Its
    imputation argument is always explicitly `"none"`, because any
    requested imputation has already occurred.
2.  `G.inverse()` tests the GRM. If necessary, the script tries
    `G.tuneup(blend = TRUE)`, then `G.tuneup(bend = TRUE)` on the
    originally computed GRM.
3.  If `n_clusters` is set, k-means uses the first `cluster_n_pcs` GRM
    PCs. `cluster_n_pcs` may be greater than two. Because the number of
    clusters is unknown until the clustering decision is made, cluster
    generation is optional and downstream scripts can instead derive
    populations from the GRM when configured to do so.

### Outputs

| File | Content |
|----|----|
| `_grm.csv` | The GRM: first column `GenoID`, then one column per genotype |
| `_clusters.csv` | Only if `n_clusters` is set: `GenoID`, every PC used (`PC1` through `PCn`), and `Population` |
| `_marker_weights_used.csv` | Only when weights are supplied: post-QC marker IDs and exact weights passed to `Gmatrix()` |
| `_sample_qc.csv` | Missingness, keep flag, and failure reason for every genotype |
| `_marker_qc.csv` | Call rate, allele frequency, MAF, heterozygosity, polymorphism, keep flag, and failure reason; unsupported/multiallelic IDs are also listed |
| `_marker_map.csv` | Post-QC marker IDs, chromosome, position, REF, and ALT |
| `_qc_report.txt` | Concise counts and thresholds |
| `_imputed_genotypes.csv` | Optional final marker-by-genotype matrix when `write_imputed_genotypes <- TRUE` |
| `_beagle.log` | Beagle command and log when Beagle actually runs |
| `_grm_info.txt` | Input format, imputation and NA counts, `G.inverse()` checks, tune-up, diagnostics, and [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) |

Example of `_grm.csv` for 10 inbred lines, with a VanRaden GRM after
blending. With inbred lines, the diagonal is close to 2.

| GenoID | G001   | G002   | G003   | G004   | G005   | G006   | G007   | G008   | G009   | G010   |
|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| G001   | 1.973  | -0.219 | -0.030 | -0.242 | -0.156 | -0.202 | -0.189 | -0.304 | -0.382 | -0.230 |
| G002   | -0.219 | 1.913  | -0.216 | -0.103 | -0.245 | -0.277 | -0.160 | -0.172 | -0.319 | -0.180 |
| G003   | -0.030 | -0.216 | 1.751  | -0.205 | -0.353 | -0.158 | -0.290 | -0.163 | -0.075 | -0.241 |
| G004   | -0.242 | -0.103 | -0.205 | 1.929  | -0.288 | -0.224 | -0.225 | -0.147 | -0.238 | -0.238 |
| G005   | -0.156 | -0.245 | -0.353 | -0.288 | 2.047  | -0.207 | -0.194 | -0.365 | -0.034 | -0.186 |
| G006   | -0.202 | -0.277 | -0.158 | -0.224 | -0.207 | 1.920  | -0.185 | -0.266 | -0.205 | -0.177 |
| G007   | -0.189 | -0.160 | -0.290 | -0.225 | -0.194 | -0.185 | 1.890  | -0.315 | -0.102 | -0.212 |
| G008   | -0.304 | -0.172 | -0.163 | -0.147 | -0.365 | -0.266 | -0.315 | 1.978  | -0.113 | -0.113 |
| G009   | -0.382 | -0.319 | -0.075 | -0.238 | -0.034 | -0.205 | -0.102 | -0.113 | 1.885  | -0.398 |
| G010   | -0.230 | -0.180 | -0.241 | -0.238 | -0.186 | -0.177 | -0.212 | -0.113 | -0.398 | 1.996  |

This file is also the GRM input of the CV scripts. A GRM from another
source can be used if it has the same layout: a square, symmetric matrix
with genotype IDs in the first column and as column names.

A marker set chosen using the trait, for example from GWAS results on
the same genotypes, carries phenotype information into the GRM, even
without marker weights.

## Cross-validation scripts

### Inputs

- The GRM from `grm_builder.R`, or a compatible GRM from another source.
- A table of stage-1 BLUEs, with BLUE weights passed to
  `asreml(weights = )`:
  - The weight column (`weight_col`), for example Smith weights, is used
    where it has a value.
  - Where it has none, the weight is 1/SE^2, from the SE column
    (`se_col`).
  - A row with a trait value but neither a weight nor an SE stops the
    run.
- Rows with no trait value are kept and predicted. They do not
  contribute to the fit or to the metrics.
- `gblup_cv_single.R` only: `weight_col <- NULL` means no BLUE weights,
  and the SE column is ignored. This also needs
  `fix_stage1_error <- FALSE` and `nonmarker_effects <- FALSE`. The
  residual variance is then estimated, so the analysis no longer treats
  the stage-1 errors as known.
- `gblup_cv_gxe.R` does not support a no-weight analysis. It requires a
  weight or SE for every row with an observed trait value.
- Every genotype in the BLUE table must be in the GRM. Genotypes in the
  GRM with no row in the BLUE table are allowed, and get NA in the
  prediction tables.

Example for `gblup_cv_single.R`, one row per genotype:

| GenoID | NBL   | Population | SE    | Weight |
|--------|-------|------------|-------|--------|
| G001   | 10.66 | 1          | 0.566 | 3.122  |
| G002   | 10.92 | 1          | 0.385 | 6.747  |
| G003   | 10.38 | 2          | 0.457 | 4.788  |
| G004   | 8.96  | 1          | 0.461 | NA     |
| G005   | 9.49  | 2          | 0.449 | 4.960  |
| G006   | 8.18  | 2          | 0.483 | 4.287  |
| G007   | 11.18 | 1          | 0.499 | 4.016  |
| G008   | NA    | 2          | NA    | NA     |
| G009   | 8.41  | 1          | 0.417 | 5.751  |
| G010   | 9.90  | 2          | 0.528 | 3.587  |

- G004 has no weight, so its weight is 1/0.461² = 4.705.
- G008 has no trait value, so it needs no weight or SE. It is predicted,
  but not used in the fit or the metrics.
- `Population` is only needed if `strata_col` is set. A file with SE
  only, or Weight only, also works.

Example for `gblup_cv_gxe.R`, one row per genotype x environment:

| GenoID | Env | NBL   | Population | SE    | Weight |
|--------|-----|-------|------------|-------|--------|
| G001   | E1  | 10.02 | 1          | 0.400 | 6.250  |
| G001   | E2  | 11.68 | 1          | 0.459 | 4.747  |
| G002   | E1  | 11.01 | 1          | 0.448 | 4.982  |
| G002   | E2  | 10.27 | 1          | 0.363 | 7.589  |
| G003   | E1  | 9.19  | 2          | 0.584 | 2.932  |
| G003   | E2  | NA    | 2          | NA    | NA     |
| G004   | E1  | 8.37  | 1          | 0.422 | 5.615  |
| G004   | E2  | 11.97 | 1          | 0.568 | 3.100  |
| G005   | E1  | 10.12 | 2          | 0.346 | NA     |
| G005   | E2  | 10.81 | 2          | 0.304 | 10.821 |

- Each genotype keeps the same `Population` in every environment.
- G003 in E2 has no trait value: it is predicted, but not used in the
  fit or the metrics.
- G005 in E1 has no weight, so its weight is 1/0.346² = 8.353.

### Models

| Script | Model |
|----|----|
| `gblup_cv_single.R` | y = mu + g + r + e |
| `gblup_cv_gxe.R` | y = mu + E + g + r + gE + rE + e, with environment E fixed |

| Term | Distribution | Included |
|----|----|----|
| g, genomic | N(0, G Vg) | Always |
| gE, genomic G x E | N(0, I_E ⊗ G Vge) | G x E script |
| r, non-marker genotype effect | N(0, I Vr) | If `nonmarker_effects <- TRUE` |
| rE, non-marker G x E | N(0, I Vre) | G x E script, if `nonmarker_effects <- TRUE` |
| e, stage-1 error | N(0, diag(1/w) s) | Always |

With `fix_stage1_error <- TRUE`, s = 1, so the stage-1 error variances
are fixed at 1/w. The terms r and rE hold genetic variation the markers
do not capture.

### Evaluation modes

`evaluation_mode <- "repeated_cv"` uses the fold design below.
`evaluation_mode <- "fixed_split"` uses the supplied training and target
files once, as described in the STPGA section. Fixed-split target
predictions remain out-of-sample; training predictions are in-sample and
must not be used to claim predictive performance.

| Scheme | Folds assigned to | What is predicted | Prediction |
|----|----|----|----|
| Single script | Genotypes, within stratum | Each genotype, once per iteration | mu + g |
| CV1 | Genotypes, within stratum | Untested genotypes | mu + E + g + gE |
| CV2 | BLUEs, within environment x stratum | Tested genotypes in environments where they were held out | mu + E + g + r + gE |

Under CV2, a genotype with no other BLUE in training is predicted as in
CV1, and is counted in `Untested_val_rows`. Folds have no population
stratification only when `strata_col <- NULL` and
`derive_strata_from_grm <- FALSE`; when GRM derivation is enabled, the
derived genomic clusters are used as the strata.

### Outputs

| File | Content |
|----|----|
| `_metrics.csv` | Metrics per iteration and fold. Fold 0 holds metrics on predictions pooled within an iteration |
| `_summary.csv` | Mean and SD of each metric across iterations |
| `_cv_predictions.csv` | Out-of-sample GEBV per genotype, per iteration, and their mean. Single: mu + g. G x E: mu + mean of E + g. Under CV2, averaged over the folds where a BLUE of the genotype was held out |
| `_cv_predictions_by_env.csv` | G x E only: observed BLUE and out-of-sample prediction for each row |
| `_train_predictions.csv`, `_train_se.csv` | In-sample predictions of the `vm` term and their SE, averaged over the folds where the genotype was in training. Not for assessing predictive performance |
| `_intercepts.csv` | Intercept of each fitted model |
| `_run_info.txt` | Settings, input checks and package versions |
| `_split_ids.csv` | Only in fixed-split mode: every included genotype and its training or target role |
| `_derived_strata.csv` | Only when `derive_strata_from_grm <- TRUE`: genomic PC coordinates and the cluster used for CV stratification |
| `_derived_strata_selection.csv` | Only when GRM strata are derived: candidate cluster counts, silhouette diagnostics, and the selected solution |
| `_derived_strata_pc_variance.csv` | Only when GRM strata are derived: the GRM eigenvalue profile and PCs used for clustering |

### Metrics

| Metric | Definition |
|----|----|
| `Vg`, `Vge`, `Vr`, `Vre` | REML estimates from the training set of the fold |
| `Residual_scale` | s\. Equals 1 when `fix_stage1_error <- TRUE` |
| `Ve` | s × mean(1/w): the mean stage-1 error variance. With `weight_col <- NULL`, the residual variance |
| `PVE` | Proportion of variance from g. Single: Vg Dg / (Vg Dg + Vr + Ve). G x E: Vg Dg / ((Vg + Vge) Dg + Vr + Vre + Ve). Dg = mean(diag(G)) − mean(G) |
| `Heritability` | Cullis heritability: 1 − avsed² / (2 varG), with avsed the mean SED from [`predict()`](https://rdrr.io/r/stats/predict.html) on the `vm` term, and varG = Vg × (mean diagonal − mean off-diagonal of G). G x E: on a genotype-mean basis |
| `Heritability_min` | Same, with the largest SED |
| `Heritability_max` | Same, with the smallest SED |
| `Heritability_SE`, `_min_SE`, `_max_SE` | Delta-method SE from `vpredict()`. Covers uncertainty in Vg only; the SED is treated as a constant |
| `Heritability_unscaled` | 1 − avsed² / (2 Vg), for reference. Not used for `Accuracy` |
| `GOF` | In-sample correlation of observed BLUEs with fitted values, including r. G x E: within environment, then averaged |
| `PredictionAbility` | Correlation of observed BLUEs with out-of-sample predictions. G x E: within environment, then averaged |
| `Accuracy` | PredictionAbility / √Heritability. G x E with `env_h2`: mean over environments of r_j / √h2_j |
| `RMSE` | Root mean squared error of out-of-sample predictions |
| `Relative_RMSE` | RMSE as a percentage of the absolute trait mean; `NA` when the mean is zero |
| `Standardized_RMSE` | RMSE divided by the SD of observed BLUEs. G x E: SD after removing environment means |
| `*_CV_allfolds` | The metric on predictions pooled within an iteration. `Accuracy_CV_allfolds` uses the mean heritability of the folds |
| `PA_CV_env_*`, `Accuracy_CV_env_*` | G x E only: pooled, within one environment |
| `PA_CV_within_*` | Pooled, within one supplied or GRM-derived stratum, when at least `min_stratum_n` genotypes have usable observed/out-of-sample prediction pairs |
| `Untested_val_rows` | G x E only: held-out BLUEs whose genotype has no BLUE in training |
| `Converged` | 1 if ASReml reported convergence, otherwise 0 |

Correlations and RMSE count each BLUE equally: BLUE weights are not used
in them.
