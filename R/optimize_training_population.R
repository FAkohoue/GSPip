#' Optimize a genomic training population
#'
#' Selects a fixed-size training set with replicated STPGA searches and
#' compares the selected design with equal-size random training sets.
#'
#' @param grm_file Path to a square GRM CSV whose first column contains genotype IDs.
#' @param output_dir Directory receiving optimization outputs.
#' @param candidate_ids_file Optional candidate-ID CSV; `NULL` derives candidates from the GRM.
#' @param target_ids_file Target-ID CSV required for targeted optimization.
#' @param run_label Prefix used for output files.
#' @param id_col Genotype-ID column in candidate and target files.
#' @param optimization_mode Either `"targeted"` or `"no_test"`.
#' @param training_size Number of candidate genotypes selected.
#' @param criterion STPGA design criterion.
#' @param lambda Positive ridge parameter passed to STPGA.
#' @param stpga_n_pcs Optional fixed PC count; `NULL` selects by cumulative variance.
#' @param stpga_pc_variance_target,stpga_max_pcs Automatic PC-selection controls.
#' @param ga_population,ga_elite,ga_iterations,ga_min_before_stop Genetic-algorithm controls.
#' @param ga_mutation_prob,ga_mutation_intensity Mutation controls.
#' @param ga_replicates,ga_seed,ga_cores Replication, seed, and parallel controls.
#' @param random_replicates,random_seed Random benchmark controls.
#' @param overwrite_outputs Permit replacement of outputs with the same run label.
#'
#' @return Invisibly, a `gspip_training_optimization` list containing selected
#'   IDs, diagnostics, benchmark results, and output paths.
#' @export
optimize_training_population <- function(
    grm_file, output_dir, candidate_ids_file = NULL, target_ids_file = NULL,
    run_label = "STPGA_training", id_col = "GenoID",
    optimization_mode = c("targeted", "no_test"), training_size = 200L,
    criterion = c("PEVMEAN", "PEVMEAN2", "PEVMAX", "CDMEAN",
                  "CDMEAN2", "CDMAX"), lambda = 1e-6,
    stpga_n_pcs = NULL, stpga_pc_variance_target = 0.95,
    stpga_max_pcs = 100L, ga_population = 100L, ga_elite = 5L,
    ga_iterations = 500L, ga_min_before_stop = 200L,
    ga_mutation_prob = 0.8, ga_mutation_intensity = 1,
    ga_replicates = 20L, ga_seed = 123L, ga_cores = 1L,
    random_replicates = 100L, random_seed = 456L,
    overwrite_outputs = FALSE) {
  optimization_mode <- match.arg(optimization_mode)
  criterion <- match.arg(criterion)
  #===============================================================================
  # PACKAGES
  #===============================================================================
  required_pkgs <- c("data.table", "STPGA")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Please install: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
  }
  stpga_api <- c("GenAlgForSubsetSelection", "GenAlgForSubsetSelectionNoTest",
                 "PEVMEAN", "PEVMEAN2", "PEVMAX", "CDMEAN", "CDMEAN2", "CDMAX")
  missing_api <- setdiff(stpga_api, getNamespaceExports("STPGA"))
  if (length(missing_api)) {
    stop("The installed STPGA version does not expose the API required by this script: ",
         paste(missing_api, collapse = ", "),
         ". Install the CRAN STPGA 5.2.1 release or a compatible version.", call. = FALSE)
  }
  
  #===============================================================================
  # FUNCTIONS
  #===============================================================================
  log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
  output_path <- function(suffix) file.path(output_dir, paste0(run_label, "_", suffix))
  show_ids <- function(ids, n = 5) {
    paste0(paste(utils::head(ids, n), collapse = ", "), if (length(ids) > n) ", ..." else "")
  }
  
  check_settings <- function() {
    scalar_logical <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
    whole_number <- function(x, lower = 0) {
      is.numeric(x) && length(x) == 1L && is.finite(x) && x >= lower && x == floor(x)
    }
    nonempty_string <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    nullable_string <- function(x) is.null(x) || nonempty_string(x)
  
    if (!nonempty_string(grm_file) || !nonempty_string(output_dir) || !nonempty_string(run_label) ||
        !nonempty_string(id_col) || !nullable_string(candidate_ids_file) ||
        !nullable_string(target_ids_file)) {
      stop("File, run-label and ID-column settings must be non-empty strings where required.", call. = FALSE)
    }
    if (grepl('[<>:"/\\\\|?*]', run_label)) stop("run_label contains reserved file-name characters.", call. = FALSE)
    if (!nonempty_string(optimization_mode) || !(optimization_mode %in% c("targeted", "no_test"))) {
      stop('optimization_mode must be "targeted" or "no_test".', call. = FALSE)
    }
    if (optimization_mode == "targeted" && !nonempty_string(target_ids_file)) {
      stop("target_ids_file is required for targeted optimization.", call. = FALSE)
    }
    allowed_criteria <- c("PEVMEAN", "PEVMEAN2", "PEVMAX", "CDMEAN", "CDMEAN2", "CDMAX")
    if (!is.character(criterion) || length(criterion) != 1L || !(criterion %in% allowed_criteria)) {
      stop("criterion must be one of: ", paste(allowed_criteria, collapse = ", "), call. = FALSE)
    }
    if (!whole_number(training_size, 2)) stop("training_size must be a whole number of at least 2.", call. = FALSE)
    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
      stop("lambda must be one positive finite number.", call. = FALSE)
    }
    if (!is.null(stpga_n_pcs) && !whole_number(stpga_n_pcs, 1)) {
      stop("stpga_n_pcs must be NULL or a positive whole number.", call. = FALSE)
    }
    if (!is.numeric(stpga_pc_variance_target) || length(stpga_pc_variance_target) != 1L ||
        !is.finite(stpga_pc_variance_target) || stpga_pc_variance_target <= 0 ||
        stpga_pc_variance_target > 1) {
      stop("stpga_pc_variance_target must be above 0 and at most 1.", call. = FALSE)
    }
    integer_settings <- list(
      stpga_max_pcs = c(stpga_max_pcs, 1), ga_population = c(ga_population, 2),
      ga_elite = c(ga_elite, 1), ga_iterations = c(ga_iterations, 1),
      ga_min_before_stop = c(ga_min_before_stop, 1), ga_replicates = c(ga_replicates, 1),
      ga_cores = c(ga_cores, 1), random_replicates = c(random_replicates, 0),
      ga_seed = c(ga_seed, 0), random_seed = c(random_seed, 0)
    )
    bad <- names(integer_settings)[!vapply(integer_settings, function(z) whole_number(z[1], z[2]), logical(1))]
    if (length(bad)) stop("Invalid whole-number setting(s): ", paste(bad, collapse = ", "), call. = FALSE)
    if (ga_elite >= ga_population) stop("ga_elite must be smaller than ga_population.", call. = FALSE)
    if (ga_min_before_stop > ga_iterations) stop("ga_min_before_stop cannot exceed ga_iterations.", call. = FALSE)
    if (!is.numeric(ga_mutation_prob) || length(ga_mutation_prob) != 1L ||
        !is.finite(ga_mutation_prob) || ga_mutation_prob < 0 || ga_mutation_prob > 1) {
      stop("ga_mutation_prob must be from 0 to 1.", call. = FALSE)
    }
    if (!is.numeric(ga_mutation_intensity) || length(ga_mutation_intensity) != 1L ||
        !is.finite(ga_mutation_intensity) || ga_mutation_intensity < 0) {
      stop("ga_mutation_intensity must be a non-negative finite number.", call. = FALSE)
    }
    if (ga_seed + ga_replicates > .Machine$integer.max ||
        random_seed + random_replicates > .Machine$integer.max) {
      stop("Seeds plus replicate counts must remain within R's integer range.", call. = FALSE)
    }
    if (.Platform$OS.type == "windows" && ga_cores != 1L) {
      stop("STPGA parallel GA execution is not portable on Windows; set ga_cores <- 1.", call. = FALSE)
    }
    if (!scalar_logical(overwrite_outputs)) stop("overwrite_outputs must be TRUE or FALSE.", call. = FALSE)
  }
  
  check_output_free <- function() {
    if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create output_dir: ", output_dir, call. = FALSE)
    }
    if (file.access(output_dir, 2) != 0) stop("output_dir is not writable: ", output_dir, call. = FALSE)
    suffixes <- c("selected_training.csv", "target_ids.csv", "candidate_ids.csv",
                  "ga_replicates.csv", "ga_training_sets.csv", "ga_convergence.csv",
                  "random_benchmark.csv", "random_training_sets.csv", "pc_variance.csv", "run_info.txt")
    existing <- suffixes[file.exists(vapply(suffixes, output_path, character(1)))]
    if (!overwrite_outputs && length(existing)) {
      stop("Output file(s) already exist: ", paste(existing, collapse = ", "),
           ". Change run_label or set overwrite_outputs <- TRUE.", call. = FALSE)
    }
  }
  
  read_grm <- function(path) {
    if (!file.exists(path)) stop("GRM file not found: ", path, call. = FALSE)
    x <- as.data.frame(data.table::fread(path, header = TRUE, check.names = FALSE,
                             colClasses = list(character = 1L)))
    if (ncol(x) < 3L) stop("The GRM must contain an ID column and at least two matrix columns.", call. = FALSE)
    ids <- as.character(x[[1]])
    G <- as.matrix(x[, -1, drop = FALSE])
    if (!is.numeric(G) || any(!is.finite(G))) stop("The GRM must contain only finite numeric values.", call. = FALSE)
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("GRM row IDs must be non-empty and unique.", call. = FALSE)
    if (nrow(G) != ncol(G) || !setequal(ids, colnames(G))) stop("The GRM must be square with matching row/column IDs.", call. = FALSE)
    rownames(G) <- ids
    G <- G[ids, ids, drop = FALSE]
    if (!isSymmetric(unname(G), tol = 1e-6)) stop("The GRM is not symmetric.", call. = FALSE)
    G
  }
  
  read_ids <- function(path, label) {
    if (is.null(path)) return(NULL)
    if (!file.exists(path)) stop(label, " file not found: ", path, call. = FALSE)
    x <- as.data.frame(data.table::fread(path, header = TRUE, check.names = FALSE,
                             colClasses = "character"))
    if (!id_col %in% names(x)) stop(label, " file is missing column '", id_col, "'.", call. = FALSE)
    ids <- as.character(x[[id_col]])
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
      stop(label, " IDs must be non-empty and unique.", call. = FALSE)
    }
    ids
  }
  
  make_pc_predictors <- function(G, n_pcs, variance_target, max_pcs) {
    eig <- eigen(G, symmetric = TRUE)
    eigen_tol <- max(1, max(abs(eig$values))) * nrow(G) * .Machine$double.eps * 100
    if (min(eig$values) < -eigen_tol) {
      stop(sprintf("The design GRM is not positive semidefinite (smallest eigenvalue %.6g).",
                   min(eig$values)), call. = FALSE)
    }
    positive <- ifelse(eig$values > eigen_tol, eig$values, 0)
    if (sum(positive) <= 0) stop("The design GRM has no positive eigenvalues.", call. = FALSE)
    percent <- positive / sum(positive)
    cumulative <- cumsum(percent)
    available <- min(nrow(G), sum(positive > 0))
    if (is.null(n_pcs)) {
      chosen <- which(cumulative >= variance_target)[1]
      chosen <- min(chosen, max_pcs, available)
    } else {
      if (n_pcs > available) stop("stpga_n_pcs exceeds the positive PCs available in the design GRM.", call. = FALSE)
      chosen <- n_pcs
    }
    chosen <- max(1L, as.integer(chosen))
    used <- seq_len(chosen)
    P <- sweep(eig$vectors[, used, drop = FALSE], 2L, sqrt(positive[used]), `*`)
    rownames(P) <- rownames(G)
    colnames(P) <- paste0("PC", used)
    info <- data.frame(
      PC = paste0("PC", seq_along(eig$values)), Eigenvalue = eig$values,
      Percent = 100 * percent, Cumulative_percent = 100 * cumulative,
      Used = seq_along(eig$values) <= chosen, stringsAsFactors = FALSE
    )
    list(P = P, info = info, n_pcs = chosen, cumulative = 100 * cumulative[chosen])
  }
  
  criterion_value <- function(train, test, P) {
    fun <- getExportedValue("STPGA", criterion)
    value <- fun(Train = train, Test = test, P = P, lambda = lambda, C = NULL)
    value <- as.numeric(value)
    if (length(value) != 1L || !is.finite(value)) stop("The STPGA criterion returned a non-finite value.", call. = FALSE)
    value
  }
  
  #===============================================================================
  # READ AND PREPARE THE DESIGN
  #===============================================================================
  check_settings()
  check_output_free()
  G_all <- read_grm(grm_file)
  grm_ids <- rownames(G_all)
  target_ids <- if (optimization_mode == "targeted") read_ids(target_ids_file, "Target") else NULL
  candidate_ids <- read_ids(candidate_ids_file, "Candidate")
  
  if (!is.null(target_ids)) {
    absent <- setdiff(target_ids, grm_ids)
    if (length(absent)) stop("Target ID(s) absent from the GRM: ", show_ids(absent), call. = FALSE)
  }
  if (is.null(candidate_ids)) candidate_ids <- setdiff(grm_ids, target_ids)
  absent <- setdiff(candidate_ids, grm_ids)
  if (length(absent)) stop("Candidate ID(s) absent from the GRM: ", show_ids(absent), call. = FALSE)
  if (length(intersect(candidate_ids, target_ids))) stop("Candidate and target IDs must be disjoint.", call. = FALSE)
  if (training_size >= length(candidate_ids)) {
    stop("training_size must be smaller than the candidate population (", length(candidate_ids), ").", call. = FALSE)
  }
  if (optimization_mode == "targeted" && length(target_ids) < 1L) stop("The targeted design has no target IDs.", call. = FALSE)
  
  design_ids <- if (optimization_mode == "targeted") c(candidate_ids, target_ids) else candidate_ids
  G_design <- G_all[design_ids, design_ids, drop = FALSE]
  pc <- make_pc_predictors(G_design, stpga_n_pcs, stpga_pc_variance_target, stpga_max_pcs)
  P <- pc$P
  P_ga <- if (optimization_mode == "targeted") P else P[candidate_ids, , drop = FALSE]
  
  log_msg(sprintf("%s optimization: %d candidates, %s, select %d; %s; %d PCs (%.1f%% cumulative variance)",
                  optimization_mode, length(candidate_ids),
                  if (optimization_mode == "targeted") paste(length(target_ids), "targets") else "test = candidate complement",
                  training_size, criterion, pc$n_pcs, pc$cumulative))
  
  #===============================================================================
  # GENETIC-ALGORITHM REPLICATES
  #===============================================================================
  replicate_rows <- vector("list", ga_replicates)
  training_rows <- vector("list", ga_replicates)
  convergence_rows <- vector("list", ga_replicates)
  
  for (replicate in seq_len(ga_replicates)) {
    seed <- ga_seed + replicate
    set.seed(seed)
    log_msg("GA replicate ", replicate, " of ", ga_replicates, " (seed ", seed, ")")
    args <- list(
      P = P_ga, ntoselect = training_size, npop = ga_population, nelite = ga_elite,
      keepbest = TRUE, tabu = TRUE, tabumemsize = 1, mutprob = ga_mutation_prob,
      mutintensity = ga_mutation_intensity, niterations = ga_iterations,
      minitbefstop = ga_min_before_stop, niterreg = 5, lambda = lambda,
      plotiters = FALSE, errorstat = criterion, C = NULL, mc.cores = ga_cores,
      InitPop = NULL, tolconv = 1e-7, Vg = NULL, Ve = NULL, Fedorov = FALSE
    )
    if (optimization_mode == "targeted") {
      args$Candidates <- candidate_ids
      args$Test <- target_ids
      fit <- do.call(STPGA::GenAlgForSubsetSelection, args)
      objective_test <- target_ids
    } else {
      fit <- do.call(STPGA::GenAlgForSubsetSelectionNoTest, args)
    }
    selected <- as.character(fit[[1]])
    if (length(selected) != training_size || anyDuplicated(selected) ||
        length(setdiff(selected, candidate_ids))) {
      stop("STPGA returned an invalid training set in replicate ", replicate, ".", call. = FALSE)
    }
    if (optimization_mode == "no_test") objective_test <- setdiff(candidate_ids, selected)
    objective <- criterion_value(selected, objective_test, P_ga)
    trajectory <- suppressWarnings(as.numeric(fit[[length(fit)]]))
    trajectory <- trajectory[is.finite(trajectory)]
    replicate_rows[[replicate]] <- data.frame(
      Replicate = replicate, Seed = seed, Objective = objective,
      GA_iterations_recorded = length(trajectory), stringsAsFactors = FALSE
    )
    training_rows[[replicate]] <- data.frame(
      Replicate = replicate, Seed = seed, Objective = objective, GenoID = selected,
      stringsAsFactors = FALSE
    )
    convergence_rows[[replicate]] <- data.frame(
      Replicate = rep(replicate, length(trajectory)),
      Iteration = seq_along(trajectory), Objective = trajectory,
      stringsAsFactors = FALSE
    )
    data.table::fwrite(do.call(rbind, replicate_rows[seq_len(replicate)]), output_path("ga_replicates.csv"))
  }
  
  replicates <- do.call(rbind, replicate_rows)
  all_training <- do.call(rbind, training_rows)
  all_convergence <- do.call(rbind, convergence_rows)
  best_rep <- replicates$Replicate[which.min(replicates$Objective)]
  selected_ids <- all_training$GenoID[all_training$Replicate == best_rep]
  evaluation_target_ids <- if (optimization_mode == "targeted") {
    target_ids
  } else {
    setdiff(candidate_ids, selected_ids)
  }
  replicates$Selected_best <- replicates$Replicate == best_rep
  
  #===============================================================================
  # RANDOM EQUAL-SIZE BENCHMARKS
  #===============================================================================
  random_summary <- list()
  random_sets <- list()
  if (random_replicates > 0L) {
    for (replicate in seq_len(random_replicates)) {
      seed <- random_seed + replicate
      set.seed(seed)
      train <- sample(candidate_ids, training_size, replace = FALSE)
      test <- if (optimization_mode == "targeted") target_ids else setdiff(candidate_ids, train)
      objective <- criterion_value(train, test, P_ga)
      random_summary[[replicate]] <- data.frame(
        Replicate = replicate, Seed = seed, Objective = objective, stringsAsFactors = FALSE
      )
      random_sets[[replicate]] <- data.frame(
        Replicate = replicate, Seed = seed, Objective = objective, GenoID = train,
        stringsAsFactors = FALSE
      )
    }
    random_summary <- do.call(rbind, random_summary)
    random_sets <- do.call(rbind, random_sets)
  } else {
    random_summary <- data.frame(Replicate = integer(), Seed = integer(), Objective = numeric())
    random_sets <- data.frame(Replicate = integer(), Seed = integer(), Objective = numeric(), GenoID = character())
  }
  
  #===============================================================================
  # EXPORT
  #===============================================================================
  data.table::fwrite(data.frame(GenoID = selected_ids, stringsAsFactors = FALSE), output_path("selected_training.csv"))
  data.table::fwrite(data.frame(GenoID = candidate_ids, stringsAsFactors = FALSE), output_path("candidate_ids.csv"))
  data.table::fwrite(data.frame(GenoID = evaluation_target_ids, stringsAsFactors = FALSE), output_path("target_ids.csv"))
  data.table::fwrite(replicates, output_path("ga_replicates.csv"))
  data.table::fwrite(all_training, output_path("ga_training_sets.csv"))
  data.table::fwrite(all_convergence, output_path("ga_convergence.csv"))
  data.table::fwrite(random_summary, output_path("random_benchmark.csv"))
  data.table::fwrite(random_sets, output_path("random_training_sets.csv"))
  data.table::fwrite(pc$info, output_path("pc_variance.csv"))
  
  random_median <- if (nrow(random_summary)) stats::median(random_summary$Objective) else NA_real_
  random_better <- if (nrow(random_summary)) mean(random_summary$Objective <= min(replicates$Objective)) * 100 else NA_real_
  writeLines(c(
    paste("Run:", run_label), paste("Date:", format(Sys.time())),
    paste("GRM:", normalizePath(grm_file)), paste("Mode:", optimization_mode),
    paste("STPGA version:", as.character(utils::packageVersion("STPGA"))),
    paste("Candidate ID file:", if (is.null(candidate_ids_file)) "none; derived from GRM"
          else normalizePath(candidate_ids_file)),
    if (optimization_mode == "targeted") paste("Target ID file:", normalizePath(target_ids_file)),
    paste("Candidate IDs:", length(candidate_ids)),
    paste("Target IDs written for fixed-split evaluation:", length(evaluation_target_ids),
          if (optimization_mode == "no_test") "(complement of the best training set)" else ""),
    paste("Training size:", training_size), paste("Criterion:", criterion), paste("Lambda:", lambda),
    paste("PCs used:", pc$n_pcs, "| cumulative variance (%):", signif(pc$cumulative, 6)),
    paste("PC selection:", if (is.null(stpga_n_pcs)) "automatic" else "manual",
          "| target:", stpga_pc_variance_target, "| maximum:", stpga_max_pcs),
    paste("GA population:", ga_population, "| elite:", ga_elite, "| iterations:", ga_iterations,
          "| replicates:", ga_replicates, "| base seed:", ga_seed),
    paste("GA minimum before stop:", ga_min_before_stop, "| mutation probability:", ga_mutation_prob,
          "| mutation intensity:", ga_mutation_intensity, "| cores:", ga_cores),
    paste("Best GA replicate:", best_rep, "| objective:", signif(min(replicates$Objective), 8)),
    paste("Random replicates:", random_replicates, "| base seed:", random_seed,
          "| median objective:", signif(random_median, 8)),
    paste("Random sets with objective <= best GA (%):", signif(random_better, 6)),
    "", capture.output(sessionInfo())
  ), output_path("run_info.txt"))
  
  log_msg(sprintf("Finished. Best replicate %d, objective %.6g; results in %s",
                  best_rep, min(replicates$Objective), normalizePath(output_dir)))

  output_suffixes <- c("selected_training.csv", "target_ids.csv",
                       "candidate_ids.csv", "ga_replicates.csv",
                       "ga_training_sets.csv", "ga_convergence.csv",
                       "random_benchmark.csv", "random_training_sets.csv",
                       "pc_variance.csv", "run_info.txt")
  output_files <- stats::setNames(vapply(output_suffixes, output_path, character(1)),
                                  output_suffixes)
  output_files <- output_files[file.exists(output_files)]
  result <- list(selected_ids = selected_ids, target_ids = evaluation_target_ids,
                 candidate_ids = candidate_ids, replicates = replicates,
                 training_sets = all_training, convergence = all_convergence,
                 random_benchmark = random_summary, random_training_sets = random_sets,
                 pc_variance = pc$info, output_files = output_files)
  class(result) <- c("gspip_training_optimization", "list")
  invisible(result)
}

#' @export
print.gspip_training_optimization <- function(x, ...) {
  best <- x$replicates[which.min(x$replicates$Objective), , drop = FALSE]
  cat("GSPip training-population optimization\n")
  cat("  Selected genotypes:", length(x$selected_ids), "\n")
  cat("  Target genotypes  :", length(x$target_ids), "\n")
  cat("  Best replicate    :", best$Replicate, "\n")
  cat("  Best objective    :", signif(best$Objective, 8), "\n")
  invisible(x)
}
