# Optimize a genomic training population

Selects a fixed-size training set with replicated STPGA searches and
compares the selected design with equal-size random training sets.

## Usage

``` r
optimize_training_population(
  grm_file,
  output_dir,
  candidate_ids_file = NULL,
  target_ids_file = NULL,
  run_label = "STPGA_training",
  id_col = "GenoID",
  optimization_mode = c("targeted", "no_test"),
  training_size = 200L,
  criterion = c("PEVMEAN", "PEVMEAN2", "PEVMAX", "CDMEAN", "CDMEAN2", "CDMAX"),
  lambda = 1e-06,
  stpga_n_pcs = NULL,
  stpga_pc_variance_target = 0.95,
  stpga_max_pcs = 100L,
  ga_population = 100L,
  ga_elite = 5L,
  ga_iterations = 500L,
  ga_min_before_stop = 200L,
  ga_mutation_prob = 0.8,
  ga_mutation_intensity = 1,
  ga_replicates = 20L,
  ga_seed = 123L,
  ga_cores = 1L,
  random_replicates = 100L,
  random_seed = 456L,
  overwrite_outputs = FALSE
)
```

## Arguments

- grm_file:

  Path to a square GRM CSV whose first column contains genotype IDs.

- output_dir:

  Directory receiving optimization outputs.

- candidate_ids_file:

  Optional candidate-ID CSV; \`NULL\` derives candidates from the GRM.

- target_ids_file:

  Target-ID CSV required for targeted optimization.

- run_label:

  Prefix used for output files.

- id_col:

  Genotype-ID column in candidate and target files.

- optimization_mode:

  Either \`"targeted"\` or \`"no_test"\`.

- training_size:

  Number of candidate genotypes selected.

- criterion:

  STPGA design criterion.

- lambda:

  Positive ridge parameter passed to STPGA.

- stpga_n_pcs:

  Optional fixed PC count; \`NULL\` selects by cumulative variance.

- stpga_pc_variance_target, stpga_max_pcs:

  Automatic PC-selection controls.

- ga_population, ga_elite, ga_iterations, ga_min_before_stop:

  Genetic-algorithm controls.

- ga_mutation_prob, ga_mutation_intensity:

  Mutation controls.

- ga_replicates, ga_seed, ga_cores:

  Replication, seed, and parallel controls.

- random_replicates, random_seed:

  Random benchmark controls.

- overwrite_outputs:

  Permit replacement of outputs with the same run label.

## Value

Invisibly, a \`gspip_training_optimization\` list containing selected
IDs, diagnostics, benchmark results, and output paths.
