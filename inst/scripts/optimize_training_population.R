################################################################################
# GSPip training-population optimization workflow template
#
# Edit SETTINGS, then run the whole script. The analysis itself is implemented
# by GSPip::optimize_training_population().
################################################################################

library(GSPip)

# Files
grm_file           <- "path/to/GRM.csv"
candidate_ids_file <- NULL
target_ids_file    <- "path/to/target_ids.csv"
output_dir         <- "path/to/results"
run_label          <- "STPGA_training"
id_col             <- "GenoID"

# Design
optimization_mode <- "targeted" # "targeted" or "no_test"
training_size     <- 200L
criterion         <- "PEVMEAN"
lambda            <- 1e-6
stpga_n_pcs              <- NULL
stpga_pc_variance_target <- 0.95
stpga_max_pcs            <- 100L

# Genetic algorithm
ga_population         <- 100L
ga_elite              <- 5L
ga_iterations         <- 500L
ga_min_before_stop    <- 200L
ga_mutation_prob      <- 0.8
ga_mutation_intensity <- 1
ga_replicates         <- 20L
ga_seed               <- 123L
ga_cores              <- 1L

# Random benchmark and output protection
random_replicates <- 100L
random_seed       <- 456L
overwrite_outputs <- FALSE

result <- optimize_training_population(
  grm_file = grm_file,
  output_dir = output_dir,
  candidate_ids_file = candidate_ids_file,
  target_ids_file = target_ids_file,
  run_label = run_label,
  id_col = id_col,
  optimization_mode = optimization_mode,
  training_size = training_size,
  criterion = criterion,
  lambda = lambda,
  stpga_n_pcs = stpga_n_pcs,
  stpga_pc_variance_target = stpga_pc_variance_target,
  stpga_max_pcs = stpga_max_pcs,
  ga_population = ga_population,
  ga_elite = ga_elite,
  ga_iterations = ga_iterations,
  ga_min_before_stop = ga_min_before_stop,
  ga_mutation_prob = ga_mutation_prob,
  ga_mutation_intensity = ga_mutation_intensity,
  ga_replicates = ga_replicates,
  ga_seed = ga_seed,
  ga_cores = ga_cores,
  random_replicates = random_replicates,
  random_seed = random_seed,
  overwrite_outputs = overwrite_outputs
)

print(result$replicates)
