#' Run single-environment GBLUP evaluation
#'
#' Fits genomic prediction models from one stage-1 BLUE per genotype and
#' evaluates them by repeated cross-validation or a fixed training/target split.
#' The ASReml backend uses REML; the open-source BGLR backend uses Bayesian
#' RKHS regression with the supplied GRM as its kernel.
#'
#' @param grm_file Path to a square GRM CSV whose first column contains genotype IDs.
#' @param pheno_file Path to the stage-1 BLUE table.
#' @param output_dir Directory receiving analysis outputs.
#' @param engine Explicit model engine: \code{"asreml"} or \code{"bglr"}.
#' @param run_label Prefix used for output files.
#' @param id_col,trait_col,strata_col,weight_col,se_col Input column names.
#'   `strata_col` and `weight_col` may be `NULL` where documented.
#' @param derive_strata_from_grm Derive CV strata from GRM PCs.
#' @param grm_strata_n_clusters Optional fixed cluster count; `NULL` selects automatically.
#' @param grm_strata_max_clusters,grm_strata_n_pcs,grm_strata_seed GRM clustering controls.
#' @param fix_stage1_error Fix the stage-1 residual scale to one.
#' @param nonmarker_effects Include the non-marker genotype effect.
#' @param k_folds,n_iterations,base_seed Cross-validation controls.
#' @param evaluation_mode Either `"repeated_cv"` or `"fixed_split"`.
#' @param training_ids_file,target_ids_file ID files required for a fixed split.
#' @param within_stratum_pa,min_stratum_n Controls within-stratum prediction ability.
#' @param overwrite_outputs Permit replacement of outputs with the same run label.
#' @param asreml_workspace Workspace passed to ASReml-R.
#' @param bglr_n_iter,bglr_burn_in,bglr_thin BGLR MCMC controls.
#' @param bglr_seed Seed used for BGLR sampling.
#' @param bglr_R2 Proportion of response variance assigned to BGLR linear
#'   predictors when its default priors are initialized.
#' @param bglr_verbose Print BGLR iteration progress.
#'
#' @return Invisibly, a `gspip_gblup_single` list containing metrics, predictions,
#'   input data used by the model, and output paths.
#' @export
run_gblup_single <- function(
    grm_file, pheno_file, output_dir, engine,
    run_label = "GBLUP_NBL",
    id_col = "GenoID", trait_col = "NBL", strata_col = "Population",
    weight_col = "Weight", se_col = "SE",
    derive_strata_from_grm = FALSE, grm_strata_n_clusters = NULL,
    grm_strata_max_clusters = 10L, grm_strata_n_pcs = 2L,
    grm_strata_seed = 123L, fix_stage1_error = TRUE,
    nonmarker_effects = TRUE, k_folds = 5L, n_iterations = 50L,
    base_seed = 0L, evaluation_mode = c("repeated_cv", "fixed_split"),
    training_ids_file = NULL, target_ids_file = NULL,
    within_stratum_pa = TRUE, min_stratum_n = 10L,
    overwrite_outputs = FALSE, asreml_workspace = 32e07,
    bglr_n_iter = 12000L, bglr_burn_in = 2000L, bglr_thin = 10L,
    bglr_seed = 123L, bglr_R2 = 0.5, bglr_verbose = FALSE) {
  engine <- match.arg(engine, c("asreml", "bglr"))
  evaluation_mode <- match.arg(evaluation_mode)
  #===============================================================================
  # PACKAGES
  #===============================================================================
  required_pkgs <- c(
    "data.table", "cluster",
    if (engine == "asreml") c("asreml", "ASRgenomics") else "BGLR"
  )
  missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    engine_note <- if (engine == "asreml") {
      " ASReml-R is licensed software available from VSNi, not from CRAN."
    } else {
      " Install BGLR with install.packages(\"BGLR\")."
    }
    stop("Please install: ", paste(missing_pkgs, collapse = ", "), ".", engine_note,
         call. = FALSE)
  }
  # Bind ASReml model-term constructors in the function environment so formulas
  # can resolve them without attaching ASReml while the package is loaded.
  if (engine == "asreml") {
    asreml_namespace <- asNamespace("asreml")
    vm  <- get("vm",  envir = asreml_namespace, inherits = FALSE)
    ide <- get("ide", envir = asreml_namespace, inherits = FALSE)
    rm(asreml_namespace)
  }
  
  #===============================================================================
  # FUNCTIONS
  #===============================================================================
  # General ----------------------------------------------------------------------
  log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
  
  show_ids <- function(ids, n = 5) {
    paste0(paste(utils::head(ids, n), collapse = ", "), if (length(ids) > n) ", ..." else "")
  }
  
  output_path <- function(suffix) file.path(output_dir, paste0(run_label, "_", suffix))
  
  check_output_free <- function() {
    if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create output_dir: ", output_dir, call. = FALSE)
    }
    if (file.access(output_dir, 2) != 0) stop("output_dir is not writable: ", output_dir, call. = FALSE)
    suffixes <- c("metrics.csv", "summary.csv", "cv_predictions.csv", "train_predictions.csv",
                  "train_se.csv", "intercepts.csv", "derived_strata.csv",
                  "derived_strata_selection.csv", "derived_strata_pc_variance.csv", "split_ids.csv",
                  "run_info.txt")
    existing <- suffixes[file.exists(vapply(suffixes, output_path, character(1)))]
    if (!overwrite_outputs && length(existing) > 0) {
      stop("Output file(s) already exist for run_label '", run_label, "': ",
           paste(existing, collapse = ", "),
           ". Change run_label or set overwrite_outputs <- TRUE.", call. = FALSE)
    }
  }
  
  # Inputs -------------------------------------------------------------------------
  check_settings <- function() {
    scalar_logical <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
    whole_number <- function(x, lower = 0) {
      is.numeric(x) && length(x) == 1L && is.finite(x) && x >= lower && x == floor(x)
    }
    nonempty_string <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    nullable_string <- function(x) is.null(x) || nonempty_string(x)
  
    if (!nonempty_string(grm_file) || !nonempty_string(pheno_file) || !nonempty_string(output_dir)) {
      stop("grm_file, pheno_file and output_dir must be non-empty paths.", call. = FALSE)
    }
    if (!nonempty_string(run_label) || grepl('[<>:"/\\\\|?*]', run_label)) {
      stop("run_label must be non-empty and cannot contain file-name separators or reserved characters.",
           call. = FALSE)
    }
    if (!all(vapply(list(id_col, trait_col, se_col), nonempty_string, logical(1))) ||
        !all(vapply(list(strata_col, weight_col), nullable_string, logical(1)))) {
      stop("Column settings must be non-empty names (or NULL where documented).", call. = FALSE)
    }
    active_columns <- c(id_col, trait_col, strata_col, weight_col, if (!is.null(weight_col)) se_col)
    if (anyDuplicated(active_columns)) stop("Column settings must refer to distinct columns.", call. = FALSE)
    if (!scalar_logical(fix_stage1_error)) {
      stop("fix_stage1_error must be TRUE or FALSE.", call. = FALSE)
    }
    if (!scalar_logical(nonmarker_effects)) {
      stop("nonmarker_effects must be TRUE or FALSE.", call. = FALSE)
    }
    if (engine == "bglr" && fix_stage1_error) {
      stop("The BGLR backend estimates its residual variance scale. Set ",
           "fix_stage1_error <- FALSE.", call. = FALSE)
    }
    if (engine == "bglr" && nonmarker_effects) {
      stop("The BGLR backend does not add a BLUE-level identity effect because ",
           "one record per genotype cannot identify it separately from residual error. ",
           "Set nonmarker_effects <- FALSE.", call. = FALSE)
    }
    if (!nonempty_string(evaluation_mode) || !evaluation_mode %in% c("repeated_cv", "fixed_split")) {
      stop('evaluation_mode must be "repeated_cv" or "fixed_split".', call. = FALSE)
    }
    if (evaluation_mode == "fixed_split" &&
        (!nonempty_string(training_ids_file) || !nonempty_string(target_ids_file))) {
      stop("training_ids_file and target_ids_file are required for evaluation_mode <- \"fixed_split\".",
           call. = FALSE)
    }
    if (!scalar_logical(derive_strata_from_grm)) {
      stop("derive_strata_from_grm must be TRUE or FALSE.", call. = FALSE)
    }
    if (derive_strata_from_grm && !is.null(strata_col)) {
      stop("Choose one stratum source: set strata_col to a column name, or set strata_col <- NULL ",
           "with derive_strata_from_grm <- TRUE.", call. = FALSE)
    }
    if (!is.null(grm_strata_n_clusters) && !whole_number(grm_strata_n_clusters, 2)) {
      stop("grm_strata_n_clusters must be NULL (automatic) or a whole number of at least 2.",
           call. = FALSE)
    }
    if (!whole_number(grm_strata_max_clusters, 2)) {
      stop("grm_strata_max_clusters must be a whole number of at least 2.", call. = FALSE)
    }
    if (!whole_number(grm_strata_n_pcs, 1)) {
      stop("grm_strata_n_pcs must be a positive whole number.", call. = FALSE)
    }
    seed_offset <- max(grm_strata_max_clusters,
                       if (is.null(grm_strata_n_clusters)) 0 else grm_strata_n_clusters)
    if (!whole_number(grm_strata_seed, 0) || grm_strata_seed + seed_offset > .Machine$integer.max) {
      stop("grm_strata_seed must be non-negative and remain within R's integer range after adding ",
           "grm_strata_max_clusters.", call. = FALSE)
    }
    if (is.null(weight_col) && fix_stage1_error) {
      stop("No-weight analysis requested with weight_col <- NULL. Set fix_stage1_error <- FALSE ",
           "and nonmarker_effects <- FALSE so all BLUEs receive equal weight and ASReml estimates ",
           "the residual variance.", call. = FALSE)
    }
    if (nonmarker_effects && !fix_stage1_error) {
      stop("nonmarker_effects needs fix_stage1_error <- TRUE: with an estimated residual ",
           "scale, the BLUE-level non-marker term and the residual cannot be separated.", call. = FALSE)
    }
    if (!whole_number(k_folds, 2)) stop("k_folds must be a whole number of at least 2.", call. = FALSE)
    if (!whole_number(n_iterations, 1)) stop("n_iterations must be a positive whole number.", call. = FALSE)
    if (!whole_number(base_seed, 0) || base_seed + n_iterations > .Machine$integer.max) {
      stop("base_seed must be a non-negative whole number with base_seed + n_iterations within R's integer range.",
           call. = FALSE)
    }
    if (!scalar_logical(within_stratum_pa) || !whole_number(min_stratum_n, 3)) {
      stop("within_stratum_pa must be TRUE/FALSE and min_stratum_n a whole number of at least 3.",
           call. = FALSE)
    }
    if (!scalar_logical(overwrite_outputs)) stop("overwrite_outputs must be TRUE or FALSE.", call. = FALSE)
    if (engine == "asreml" &&
        (!is.numeric(asreml_workspace) || length(asreml_workspace) != 1L ||
         !is.finite(asreml_workspace) || asreml_workspace <= 0)) {
      stop("asreml_workspace must be one positive finite number.", call. = FALSE)
    }
    if (engine == "bglr") {
      .gspip_validate_bglr_controls(
        bglr_n_iter, bglr_burn_in, bglr_thin, bglr_seed, bglr_R2, bglr_verbose
      )
      max_bglr_seed <- bglr_seed + n_iterations * k_folds
      if (!is.finite(max_bglr_seed) || max_bglr_seed > .Machine$integer.max) {
        stop("bglr_seed plus the maximum iteration/fold offset exceeds R's integer range.",
             call. = FALSE)
      }
    }
  }
  
  read_grm <- function(path) {
    if (!file.exists(path)) stop("GRM file not found: ", path, call. = FALSE)
    grm <- as.data.frame(data.table::fread(path, header = TRUE, colClasses = list(character = 1L)))
    ids <- as.character(grm[[1]])
    grm <- as.matrix(grm[, -1, drop = FALSE])
    if (!is.numeric(grm)) stop("The GRM contains non-numeric values.", call. = FALSE)
    rownames(grm) <- ids
    if (nrow(grm) != ncol(grm)) {
      stop(sprintf("The GRM is not square: %d rows, %d columns.", nrow(grm), ncol(grm)), call. = FALSE)
    }
    if (nrow(grm) < 2L) stop("The GRM must contain at least two genotypes.", call. = FALSE)
    if (anyNA(ids) || any(!nzchar(ids))) stop("The GRM contains missing or empty genotype IDs.", call. = FALSE)
    if (anyDuplicated(ids)) stop("Duplicated IDs in the GRM: ", show_ids(ids[duplicated(ids)]), call. = FALSE)
    if (anyNA(colnames(grm)) || any(!nzchar(colnames(grm))) || anyDuplicated(colnames(grm))) {
      stop("GRM column IDs must be non-empty and unique.", call. = FALSE)
    }
    if (!setequal(ids, colnames(grm))) stop("GRM row and column IDs differ.", call. = FALSE)
    grm <- grm[ids, ids, drop = FALSE]
    if (any(!is.finite(grm))) stop("The GRM contains missing or non-finite values.", call. = FALSE)
    if (!isSymmetric(unname(grm), tol = 1e-6)) stop("The GRM is not symmetric.", call. = FALSE)
    eigenvalues <- eigen(grm, symmetric = TRUE, only.values = TRUE)$values
    eigen_tol <- max(1, max(abs(eigenvalues))) * nrow(grm) * .Machine$double.eps * 100
    if (min(eigenvalues) <= eigen_tol) {
      stop(sprintf("The GRM must be positive definite; smallest eigenvalue is %.4g. Tune or bend the GRM first.",
                   min(eigenvalues)), call. = FALSE)
    }
    grm
  }
  
  # Weight where present, 1/SE^2 where missing
  resolve_weights <- function(pheno) {
    if (is.null(weight_col)) {
      return(list(w = rep(1, nrow(pheno)), source = "none (weight_col is NULL)"))
    }
    observed <- !is.na(pheno[[trait_col]])
    has_w  <- weight_col %in% names(pheno)
    has_se <- se_col %in% names(pheno)
    if (!has_w && !has_se) {
      stop("The phenotype file has neither the weight column ('", weight_col, "') nor the SE column ('",
           se_col, "'). To run without either column, set weight_col <- NULL, ",
           "fix_stage1_error <- FALSE, and nonmarker_effects <- FALSE.", call. = FALSE)
    }
    w <- rep(NA_real_, nrow(pheno))
    if (has_w) {
      if (!is.numeric(pheno[[weight_col]]) && !all(is.na(pheno[[weight_col]]))) {
        stop("Weight column is not numeric: ", weight_col, call. = FALSE)
      }
      w <- as.numeric(pheno[[weight_col]])
    }
    from_se <- observed & is.na(w)
    if (any(from_se) && has_se) {
      if (!is.numeric(pheno[[se_col]]) && !all(is.na(pheno[[se_col]]))) {
        stop("SE column is not numeric: ", se_col, call. = FALSE)
      }
      se <- as.numeric(pheno[[se_col]])
      if (any(se[from_se] <= 0, na.rm = TRUE)) stop("SE values must be positive.", call. = FALSE)
      w[from_se] <- 1 / se[from_se]^2
    }
    if (anyNA(w[observed])) {
      stop(sum(is.na(w[observed])), " row(s) with a trait value have neither a weight nor an SE. First: ",
           show_ids(as.character(pheno[[id_col]][observed & is.na(w)])), call. = FALSE)
    }
    if (any(!is.finite(w[observed]) | w[observed] <= 0)) {
      stop("Weights must be positive and finite.", call. = FALSE)
    }
    w[!observed] <- median(w[observed])  # No trait value: no effect on the fit
    list(w = w, source = sprintf("%d row(s) from '%s', %d row(s) as 1/SE^2 from '%s'",
                                 sum(observed & !from_se), weight_col, sum(observed & from_se), se_col))
  }
  
  derive_grm_strata <- function(G, k, max_k, n_pcs, seed) {
    eig <- eigen(G, symmetric = TRUE)
    if (n_pcs > nrow(G)) {
      stop("grm_strata_n_pcs exceeds the number of PCs available from the GRM.", call. = FALSE)
    }
    used <- seq_len(n_pcs)
    pcs <- sweep(eig$vectors[, used, drop = FALSE], 2L, sqrt(pmax(eig$values[used], 0)), `*`)
    colnames(pcs) <- paste0("PC", used)
    positive <- eig$values > 0
    positive_total <- sum(eig$values[positive])
    variance_percent <- 100 * pmax(eig$values, 0) / positive_total
    pc_variance <- data.frame(
      PC = paste0("PC", seq_along(eig$values)), Eigenvalue = eig$values,
      Percent = variance_percent, Cumulative_percent = cumsum(variance_percent),
      Used_for_clustering = seq_along(eig$values) <= n_pcs, stringsAsFactors = FALSE
    )
    n_distinct <- nrow(unique(as.data.frame(pcs)))
    if (n_distinct < 2L || nrow(G) < 3L) {
      stop("At least three genotypes and two distinct PC positions are needed to derive GRM strata.",
           call. = FALSE)
    }
  
    if (is.null(k)) {
      upper_k <- min(max_k, nrow(G) - 1L, n_distinct)
      candidates <- seq.int(2L, upper_k)
      selection_mode <- "automatic: maximum average silhouette width"
    } else {
      if (k > min(nrow(G) - 1L, n_distinct)) {
        stop("grm_strata_n_clusters is too large for the number of genotypes or distinct PC positions.",
             call. = FALSE)
      }
      candidates <- as.integer(k)
      selection_mode <- "manual override"
    }
  
    distances <- stats::dist(pcs)
    fits <- vector("list", length(candidates))
    scores <- data.frame(
      K = candidates, Average_silhouette = NA_real_, Min_cluster_size = NA_integer_,
      Total_withinss = NA_real_, N_PCs = n_pcs,
      Cumulative_PC_variance = sum(variance_percent[used]), stringsAsFactors = FALSE
    )
    for (i in seq_along(candidates)) {
      set.seed(seed + candidates[i])
      fit <- stats::kmeans(pcs, centers = candidates[i], nstart = 25)
      fits[[i]] <- fit
      scores$Average_silhouette[i] <- mean(cluster::silhouette(fit$cluster, distances)[, "sil_width"])
      scores$Min_cluster_size[i] <- min(tabulate(fit$cluster, nbins = candidates[i]))
      scores$Total_withinss[i] <- fit$tot.withinss
    }
    if (any(!is.finite(scores$Average_silhouette))) {
      stop("Silhouette scoring returned a non-finite value; GRM strata could not be selected.", call. = FALSE)
    }
    selected_i <- which.max(scores$Average_silhouette)
    selected_k <- candidates[selected_i]
    scores$Selected <- seq_len(nrow(scores)) == selected_i
    scores$Selection_mode <- selection_mode
    cluster_id <- fits[[selected_i]]$cluster
  
    table <- data.frame(GenoID = rownames(G), as.data.frame(pcs, check.names = FALSE),
                        Stratum = paste0("GRM_cluster_", cluster_id),
                        check.names = FALSE, stringsAsFactors = FALSE)
    list(
      table = table, strata = setNames(table$Stratum, table$GenoID), pc_variance = pc_variance,
      n_pcs = n_pcs, cumulative_pc_variance = sum(variance_percent[used]),
      selection = scores, selected_k = selected_k,
      selected_silhouette = scores$Average_silhouette[selected_i], selection_mode = selection_mode
    )
  }
  
  stratum_values <- function(pheno, ids, grm_strata = NULL) {
    if (!is.null(grm_strata)) return(unname(grm_strata[ids]))
    if (is.null(strata_col)) rep("All", nrow(pheno)) else as.character(pheno[[strata_col]])
  }
  
  read_pheno <- function(path, grm_ids, k, grm_strata = NULL) {
    if (!file.exists(path)) stop("Phenotype file not found: ", path, call. = FALSE)
    pheno <- as.data.frame(data.table::fread(path, header = TRUE,
                                 colClasses = list(character = id_col)))
  
    absent <- setdiff(c(id_col, trait_col, strata_col), names(pheno))
    if (length(absent) > 0) {
      stop("Column(s) missing from the phenotype file: ", paste(absent, collapse = ", "), call. = FALSE)
    }
    ids <- as.character(pheno[[id_col]])
    if (anyNA(ids) || any(!nzchar(ids))) stop("Missing or empty genotype IDs in the phenotype file.", call. = FALSE)
    not_in_grm <- setdiff(ids, grm_ids)
    if (length(not_in_grm) > 0) {
      stop(length(not_in_grm), " phenotyped genotype(s) not in the GRM: ", show_ids(not_in_grm), call. = FALSE)
    }
    if (!is.numeric(pheno[[trait_col]])) stop("Trait column is not numeric: ", trait_col, call. = FALSE)
    if (any(!is.finite(pheno[[trait_col]][!is.na(pheno[[trait_col]])]))) {
      stop("Trait values must be finite where present: ", trait_col, call. = FALSE)
    }
    if (!is.null(strata_col) && anyNA(pheno[[strata_col]])) {
      stop("Missing values in stratum column: ", strata_col, call. = FALSE)
    }
    if (!is.null(strata_col) && any(!nzchar(as.character(pheno[[strata_col]])))) {
      stop("Empty values in stratum column: ", strata_col, call. = FALSE)
    }
    n_no_trait <- sum(is.na(pheno[[trait_col]]))
    if (n_no_trait == nrow(pheno)) stop("The trait column has no values: ", trait_col, call. = FALSE)
    if (anyDuplicated(ids)) {
      stop("Each genotype must have one row. For BLUEs per environment, use gblup_cv_gxe.R. Duplicated: ",
           show_ids(unique(ids[duplicated(ids)])), call. = FALSE)
    }
    if (nrow(pheno) < k) stop("k_folds exceeds the number of genotypes in the phenotype file.", call. = FALSE)
    if (sum(!is.na(pheno[[trait_col]])) < k) {
      stop("Fewer observed trait values than k_folds; each fold needs validation data.", call. = FALSE)
    }
    weights <- resolve_weights(pheno)
  
    out <- data.frame(GenoID = ids, y_trait = pheno[[trait_col]], w = weights$w,
                      Stratum = stratum_values(pheno, ids, grm_strata), stringsAsFactors = FALSE)
    small <- names(which(table(out$Stratum) < k))
    if (length(small) > 0) {
      warning("Stratum smaller than k_folds, so absent from some folds: ",
              paste(small, collapse = ", "), call. = FALSE)
    }
    list(data = out, weight_source = weights$source, n_no_trait = n_no_trait)
  }
  
  read_fixed_split <- function(training_path, target_path, phenotype_ids) {
    read_one <- function(path, label) {
      if (!file.exists(path)) stop(label, " ID file not found: ", path, call. = FALSE)
      x <- as.data.frame(data.table::fread(path, header = TRUE, check.names = FALSE,
                               colClasses = "character"))
      if (!id_col %in% names(x)) stop(label, " ID file is missing column '", id_col, "'.", call. = FALSE)
      ids <- as.character(x[[id_col]])
      if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
        stop(label, " IDs must be non-empty and unique.", call. = FALSE)
      }
      ids
    }
    training <- read_one(training_path, "Training")
    target <- read_one(target_path, "Target")
    overlap <- intersect(training, target)
    if (length(overlap)) stop("Training and target ID files overlap: ", show_ids(overlap), call. = FALSE)
    absent_training <- setdiff(training, phenotype_ids)
    absent_target <- setdiff(target, phenotype_ids)
    if (length(absent_training)) stop("Training ID(s) absent from the phenotype file: ",
                                      show_ids(absent_training), call. = FALSE)
    if (length(absent_target)) stop("Target ID(s) absent from the phenotype file: ",
                                    show_ids(absent_target), call. = FALSE)
    if (length(training) < 2L) stop("The fixed split needs at least two training genotypes.", call. = FALSE)
    if (length(target) < 1L) stop("The fixed split has no target genotypes.", call. = FALSE)
    list(training = training, target = target)
  }
  
  # Cross-validation and metrics -------------------------------------------------
  make_folds <- function(strata, k) {
    if (length(strata) < k) stop("Cannot create more folds than records.", call. = FALSE)
    folds <- rep(NA_integer_, length(strata))
    fold_load <- integer(k)
    for (s in sample(unique(strata))) {
      idx <- sample(which(strata == s))
      base_n <- length(idx) %/% k
      labels <- rep(seq_len(k), each = base_n)
      fold_load <- fold_load + base_n
      remainder <- length(idx) %% k
      if (remainder > 0L) {
        tie_break <- stats::runif(k)
        extra <- order(fold_load, tie_break)[seq_len(remainder)]
        labels <- c(labels, extra)
        fold_load[extra] <- fold_load[extra] + 1L
      }
      folds[idx] <- sample(labels)
    }
    if (any(tabulate(folds, nbins = k) == 0L)) {
      stop("Could not create non-empty folds; reduce k_folds.", call. = FALSE)
    }
    folds
  }
  
  metric_rows <- function(iter, fold, ...) {
    v <- c(...)
    data.frame(Iter = iter, Fold = fold, Metric = names(v), Value = unname(v), stringsAsFactors = FALSE)
  }
  
  safe_cor <- function(x, y) {
    ok <- complete.cases(x, y) & is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
    cor(x[ok], y[ok])
  }
  
  prediction_metrics <- function(obs, pred) {
    ok <- complete.cases(obs, pred) & is.finite(obs) & is.finite(pred)
    if (!any(ok)) return(c(PA = NA_real_, RMSE = NA_real_, Relative_RMSE = NA_real_,
                           Standardized_RMSE = NA_real_))
    rmse <- sqrt(mean((obs[ok] - pred[ok])^2))
    obs_mean <- mean(obs[ok])
    obs_sd <- stats::sd(obs[ok])
    c(PA = safe_cor(obs, pred),
      RMSE = rmse,
      Relative_RMSE = if (is.finite(obs_mean) && obs_mean != 0) rmse / abs(obs_mean) * 100 else NA_real_,
      Standardized_RMSE = if (is.finite(obs_sd) && obs_sd > 0) rmse / obs_sd else NA_real_)
  }
  
  accuracy_metric <- function(pa, cullis_h2) {
    if (is.finite(pa) && is.finite(cullis_h2) && cullis_h2 > 0) pa / sqrt(cullis_h2) else NA_real_
  }
  
  safe_ratio <- function(numerator, denominator) {
    if (is.finite(numerator) && is.finite(denominator) && denominator > 0) {
      numerator / denominator
    } else {
      NA_real_
    }
  }
  
  mean_finite <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  
  # Fold metrics are averaged within iteration first: SD is between CV replicates
  summarise_metrics <- function(metrics) {
    m <- data.table::as.data.table(metrics)
    m[, Level := data.table::fifelse(Fold == 0L, "pooled", "fold")]
    by_iter <- m[, .(Value = if (any(is.finite(Value))) mean(Value[is.finite(Value)]) else NA_real_),
                 by = .(Level, Metric, Iter)]
    by_iter[, {
      valid <- Value[is.finite(Value)]
      .(Mean = if (length(valid)) mean(valid) else NA_real_,
        SD = if (length(valid) > 1L) stats::sd(valid) else NA_real_,
        N_iter = length(valid))
    }, by = .(Level, Metric)]
  }
  
  # ASReml outputs ----------------------------------------------------------------
  # asreml::vpredict() addresses Vg by position (V1, V2, ...); a mismatch stops the run
  cullis_vpredict <- function(fit, sed, Vg_pos, G_scale, expected) {
    if (!all(is.finite(c(sed, G_scale, expected))) || G_scale <= 0) return(NA_real_)
    h2_formula <- as.formula(sprintf("h2 ~ 1 - ((%.15g^2) / (2 * %.15g * V%d))", sed, G_scale, Vg_pos),
                             env = globalenv())
    out <- as.data.frame(asreml::vpredict(fit, h2_formula))
    if (abs(out[1, "Estimate"] - expected) > 1e-6 * max(1, abs(expected))) {
      stop(sprintf("asreml::vpredict() gave h2 = %.6f but 1 - SED^2 / (2 * varG) = %.6f, so V%d is not Vg.",
                   out[1, "Estimate"], expected, Vg_pos), call. = FALSE)
    }
    out[1, "SE"]
  }
  
  # Cullis heritability with rescaled Vg: varG = Vg * (mean diagonal - mean off-diagonal)
  cullis_heritability <- function(fit, pred, varcomp, G) {
    avsed <- as.data.frame(pred$avsed)[, 1]  # min, mean, max
    if (length(avsed) < 3L || any(!is.finite(avsed[1:3]))) {
      stop("predict() did not return finite minimum, mean and maximum SED values.", call. = FALSE)
    }
    n <- nrow(G)
    diag_mean    <- mean(diag(G))
    offdiag_mean <- (sum(G) - sum(diag(G))) / (n * (n - 1))
    G_scale <- diag_mean - offdiag_mean
    varG    <- varcomp$Vg * G_scale
  
    if (!is.finite(varG) || varG <= 0 || !is.finite(G_scale) || G_scale <= 0) {
      warning("Cullis heritability is undefined because Vg or the GRM scale is not positive.", call. = FALSE)
      return(list(
        h2 = NA_real_, Dg = diag_mean - mean(G),
        metrics = c(Heritability = NA_real_, Heritability_SE = NA_real_,
                    Heritability_min = NA_real_, Heritability_min_SE = NA_real_,
                    Heritability_max = NA_real_, Heritability_max_SE = NA_real_,
                    Heritability_unscaled = NA_real_)
      ))
    }
  
    h2_min <- 1 - avsed[3]^2 / (2 * varG)
    h2     <- 1 - avsed[2]^2 / (2 * varG)
    h2_max <- 1 - avsed[1]^2 / (2 * varG)
  
    list(
      h2 = h2,
      Dg = diag_mean - mean(G),
      metrics = c(
        Heritability = h2,
        Heritability_SE = cullis_vpredict(fit, avsed[2], varcomp$Vg_pos, G_scale, h2),
        Heritability_min = h2_min,
        Heritability_min_SE = cullis_vpredict(fit, avsed[3], varcomp$Vg_pos, G_scale, h2_min),
        Heritability_max = h2_max,
        Heritability_max_SE = cullis_vpredict(fit, avsed[1], varcomp$Vg_pos, G_scale, h2_max),
        Heritability_unscaled = 1 - avsed[2]^2 / (2 * varcomp$Vg)
      )
    )
  }
  
  read_train_predictions <- function(pred) {
    data.frame(
      GenoID = sub("vm(GenoID, training_grm_inv.sparse)_", "", as.character(pred$pvals[, 1]), fixed = TRUE),
      GEBV   = pred$pvals[, 2],
      SE     = pred$pvals[, 3],
      stringsAsFactors = FALSE
    )
  }
  
  read_main_term <- function(re, tag, ids) {
    nms <- rownames(re)
    is_term <- startsWith(nms, tag)
    v <- re[is_term, 1][match(ids, substring(nms[is_term], nchar(tag) + 1))]
    if (anyNA(v)) stop("Could not match '", tag, "' effects to training genotypes.", call. = FALSE)
    setNames(as.numeric(v), ids)
  }
  
  extract_varcomp <- function(fit) {
    vc <- as.data.frame(summary(fit)$varcomp)
    terms <- rownames(vc)
    has_vm  <- grepl("vm(GenoID", terms, fixed = TRUE)
    has_ide <- grepl("ide(GenoID", terms, fixed = TRUE)
    get_one <- function(keep, label) {
      v <- vc[keep, 1]
      if (length(v) != 1) stop("Could not find ", label, " in: ", paste(terms, collapse = ", "), call. = FALSE)
      v
    }
    list(
      Vg = get_one(has_vm, "Vg"),
      Vg_pos = which(has_vm),
      Vr = if (nonmarker_effects) get_one(has_ide, "Vr") else 0,
      s  = if (fix_stage1_error) 1 else get_one(!has_vm & !has_ide, "the residual scale")
    )
  }
  
  extract_intercept <- function(fit) {
    fx <- as.matrix(fit$coefficients$fixed)
    mu <- fx[rownames(fx) == "(Intercept)", 1]
    if (length(mu) != 1) stop("Could not find the intercept in the fixed effects.", call. = FALSE)
    unname(mu)
  }
  
  # Run record ---------------------------------------------------------------------
  write_run_info <- function() {
    writeLines(c(
      paste("Run:", run_label),
      paste("Date:", format(Sys.time())),
      paste("GRM:", normalizePath(grm_file)),
      paste("Stage-1 BLUEs:", normalizePath(pheno_file)),
      paste("Engine:", engine),
      if (engine == "bglr") {
        paste("BGLR MCMC: nIter =", bglr_n_iter, "| burnIn =", bglr_burn_in,
              "| thin =", bglr_thin, "| retained =",
              floor(bglr_n_iter / bglr_thin) - floor(bglr_burn_in / bglr_thin),
              "| R2 =", bglr_R2, "| seed =", bglr_seed)
      },
      paste("Evaluation mode:", evaluation_mode),
      if (!is.null(fixed_split)) paste("Training ID file:", normalizePath(training_ids_file),
                                       "| genotypes:", length(fixed_split$training)),
      if (!is.null(fixed_split)) paste("Target ID file:", normalizePath(target_ids_file),
                                       "| genotypes:", length(fixed_split$target)),
      paste("Trait:", trait_col, "| Stratum:", if (is.null(strata_col)) "none" else strata_col),
      paste("Stratum source:", stratum_source),
      if (!is.null(grm_strata_result)) {
        sprintf("GRM clustering: %d clusters (%s; average silhouette %.3f) using %d PCs (cumulative variance %.1f%%), seed %s",
                grm_strata_result$selected_k, grm_strata_result$selection_mode,
                grm_strata_result$selected_silhouette, grm_strata_result$n_pcs,
                grm_strata_result$cumulative_pc_variance, grm_strata_seed)
      },
      paste("Stage-1 error:", if (fix_stage1_error) "fixed, residual = 1/w" else
            paste0(if (engine == "bglr") "posterior " else "",
                   "scale estimated, residual = s/w")),
      paste("Non-marker effect r:", if (nonmarker_effects) "yes" else "no"),
      paste("Weights:", pheno_input$weight_source),
      if (!is.null(weight_col)) paste("Weight range:", paste(signif(range(pheno$w[!is.na(pheno$y_trait)]), 4), collapse = " to ")),
      paste("Folds:", run_folds, "| Iterations:", run_iterations,
            if (evaluation_mode == "repeated_cv") paste("| Base seed:", base_seed) else ""),
      paste("Genotypes in GRM:", nrow(grm)),
      paste("Genotypes in the phenotype file:", nrow(pheno)),
      paste("Rows without a trait value (predicted, not in metrics):", sum(is.na(pheno$y_trait))),
      paste("In GRM without a phenotype row (NA in prediction tables):", length(setdiff(all_genoids, pheno$GenoID))),
      "", "Stratum sizes (genotypes):",
      capture.output(print(table(pheno$Stratum))),
      "", capture.output(sessionInfo())
    ), output_path("run_info.txt"))
  }
  
  #===============================================================================
  # STEP 1. READ AND CHECK INPUTS
  #===============================================================================
  check_settings()
  check_output_free()
  
  log_msg("Reading inputs")
  grm         <- read_grm(grm_file)
  all_genoids <- rownames(grm)
  grm_strata_result <- if (derive_strata_from_grm) {
    derive_grm_strata(grm, grm_strata_n_clusters, grm_strata_max_clusters,
                      grm_strata_n_pcs, grm_strata_seed)
  } else {
    NULL
  }
  if (!is.null(grm_strata_result)) {
    log_msg(sprintf("GRM strata: selected k = %d by %s using %d PCs (average silhouette %.3f)",
                    grm_strata_result$selected_k, grm_strata_result$selection_mode,
                    grm_strata_result$n_pcs, grm_strata_result$selected_silhouette))
  }
  stratum_source <- if (derive_strata_from_grm) {
    "GRM PC k-means clusters"
  } else if (is.null(strata_col)) {
    "none"
  } else {
    paste0("phenotype column '", strata_col, "'")
  }
  input_k <- if (evaluation_mode == "fixed_split") 2L else k_folds
  pheno_input <- read_pheno(pheno_file, all_genoids, input_k,
                            if (is.null(grm_strata_result)) NULL else grm_strata_result$strata)
  pheno       <- pheno_input$data
  fixed_split <- NULL
  if (evaluation_mode == "fixed_split") {
    fixed_split <- read_fixed_split(training_ids_file, target_ids_file, pheno$GenoID)
    keep_ids <- c(fixed_split$training, fixed_split$target)
    pheno <- pheno[pheno$GenoID %in% keep_ids, , drop = FALSE]
    observed_training_ids <- pheno$GenoID[pheno$GenoID %in% fixed_split$training &
                                           !is.na(pheno$y_trait)]
    no_training_values <- setdiff(fixed_split$training, observed_training_ids)
    if (length(no_training_values)) {
      stop("Selected training genotype(s) have no observed trait: ",
           show_ids(no_training_values), call. = FALSE)
    }
    if (sum(!is.na(pheno$y_trait[pheno$GenoID %in% fixed_split$training])) < 2L) {
      stop("The fixed training set has fewer than two observed trait values.", call. = FALSE)
    }
    split_table <- data.frame(
      GenoID = keep_ids,
      Role = c(rep("Training", length(fixed_split$training)), rep("Target", length(fixed_split$target))),
      stringsAsFactors = FALSE
    )
    data.table::fwrite(split_table, output_path("split_ids.csv"))
  }
  run_iterations <- if (evaluation_mode == "fixed_split") 1L else n_iterations
  run_folds <- if (evaluation_mode == "fixed_split") 1L else k_folds
  
  if (!is.null(grm_strata_result)) {
    data.table::fwrite(grm_strata_result$table, output_path("derived_strata.csv"))
    data.table::fwrite(grm_strata_result$selection, output_path("derived_strata_selection.csv"))
    data.table::fwrite(grm_strata_result$pc_variance, output_path("derived_strata_pc_variance.csv"))
  }
  write_run_info()
  if (evaluation_mode == "fixed_split") {
    log_msg(sprintf("Fixed split: %d training and %d target genotypes; %s; strata: %s",
                    length(fixed_split$training), length(fixed_split$target),
                    pheno_input$weight_source, stratum_source))
  } else {
    log_msg(sprintf("%d genotypes, %d to predict only; %s; strata: %s; %d iterations of %d-fold CV",
                    nrow(pheno), sum(is.na(pheno$y_trait)), pheno_input$weight_source, stratum_source,
                    run_iterations, run_folds))
  }
  
  #===============================================================================
  # STEP 2. DEFINE THE MODEL
  #===============================================================================
  # Top-level objects: predict() re-evaluates the model call
  fixed_formula <- if (engine == "asreml") y_trait ~ 1 else NULL
  random_formula <- if (engine == "asreml") {
    if (nonmarker_effects) {
      ~ vm(GenoID, training_grm_inv.sparse) + ide(GenoID)
    } else {
      ~ vm(GenoID, training_grm_inv.sparse)
    }
  } else {
    NULL
  }
  residual_dispersion <- if (engine == "asreml" && fix_stage1_error) 1 else NA
  
  #===============================================================================
  # STEP 3. PREPARE RESULT TABLES
  #===============================================================================
  empty_iter_table <- function(row_names = all_genoids) {
    matrix(NA_real_, nrow = length(row_names), ncol = run_iterations,
           dimnames = list(row_names, paste0(trait_col, "_iter_", seq_len(run_iterations))))
  }
  cv_predictions    <- empty_iter_table()
  train_predictions <- empty_iter_table()
  train_se          <- empty_iter_table()
  intercepts        <- vector("list", run_iterations)
  metrics_list      <- list()
  
  #===============================================================================
  # STEP 4. RUN THE REQUESTED EVALUATION
  #===============================================================================
  # asreml::asreml() and predict() stay at top level: predict() needs the same objects
  for (iteration in seq_len(run_iterations)) {
    log_msg(sprintf("Iteration %d of %d", iteration, run_iterations))
  
    # Assign folds
    if (evaluation_mode == "fixed_split") {
      pheno$Fold <- ifelse(pheno$GenoID %in% fixed_split$target, 1L, 0L)
    } else {
      set.seed(base_seed + iteration)
      pheno$Fold <- make_folds(pheno$Stratum, run_folds)
    }
  
    observed_by_fold <- tabulate(pheno$Fold[!is.na(pheno$y_trait)], nbins = run_folds)
    if (any(observed_by_fold == 0L)) {
      warning("Iteration ", iteration, ": fold(s) without observed validation traits: ",
              paste(which(observed_by_fold == 0L), collapse = ", "), call. = FALSE)
    }
  
    fold_train <- vector("list", run_folds)
    fold_h2    <- numeric(run_folds)
    fold_blue  <- vector("list", run_folds)
  
    for (fold in seq_len(run_folds)) {
      # Split training and validation sets
      training_data   <- pheno[pheno$Fold != fold, ]
      validation_data <- pheno[pheno$Fold == fold, ]
      training_ids    <- training_data$GenoID
      validation_ids  <- validation_data$GenoID
  
      if (sum(!is.na(training_data$y_trait)) < 2L) {
        stop(sprintf("Iteration %d, fold %d has fewer than two observed training values. Check the split or add data.",
                     iteration, fold), call. = FALSE)
      }
  
      training_grm <- grm[training_ids, training_ids, drop = FALSE]
      observed_y <- training_data$y_trait
      if (engine == "asreml") {
        training_grm_inv.sparse <- ASRgenomics::G.inverse(
          G = training_grm, sparseform = TRUE
        )$Ginv
        training_data$GenoID <- factor(training_data$GenoID, levels = training_ids)

        GBLUPcv <- asreml::asreml(
          fixed      = fixed_formula,
          random     = random_formula,
          weights    = "w",
          family     = asreml::asr_gaussian(dispersion = residual_dispersion),
          na.action  = asreml::na.method(x = "include", y = "include"),
          pworkspace = asreml_workspace,
          workspace  = asreml_workspace,
          data       = training_data
        )
        converged <- isTRUE(GBLUPcv$converge)
        if (!converged) {
          warning(sprintf("ASReml did not report convergence: iteration %d, fold %d",
                          iteration, fold), call. = FALSE)
        }

        varcomp <- extract_varcomp(GBLUPcv)
        pred_train <- predict(
          GBLUPcv, classify = "vm(GenoID, training_grm_inv.sparse)", sed = TRUE,
          pworkspace = asreml_workspace, workspace = asreml_workspace
        )
        heritability <- cullis_heritability(
          GBLUPcv, pred_train, varcomp, training_grm
        )
        fold_h2[fold] <- heritability$h2

        intercept   <- extract_intercept(GBLUPcv)
        random_coef <- as.matrix(GBLUPcv$coefficients$random)
        u_train <- read_main_term(
          random_coef, "vm(GenoID, training_grm_inv.sparse)_", training_ids
        )
        r_train <- if (nonmarker_effects) {
          read_main_term(random_coef, "ide(GenoID)_", training_ids)
        } else {
          0
        }

        u_val <- grm[validation_ids, training_ids, drop = FALSE] %*%
          solve(training_grm, u_train)
        val_pred <- intercept + as.vector(u_val)
        cv_predictions[validation_ids, iteration] <- val_pred

        Vg <- varcomp$Vg
        Vr <- varcomp$Vr
        Ve <- varcomp$s * mean(1 / training_data$w[!is.na(observed_y)])
        val <- prediction_metrics(validation_data$y_trait, val_pred)

        metrics_list[[length(metrics_list) + 1L]] <- metric_rows(
          iteration, fold,
          Vg = Vg,
          Vr = Vr,
          Residual_scale = varcomp$s,
          Ve = Ve,
          PVE = safe_ratio(Vg * heritability$Dg, Vg * heritability$Dg + Vr + Ve),
          heritability$metrics,
          GOF = safe_cor(observed_y, intercept + u_train + r_train),
          PredictionAbility = unname(val["PA"]),
          Accuracy = accuracy_metric(unname(val["PA"]), heritability$h2),
          RMSE = unname(val["RMSE"]),
          Relative_RMSE = unname(val["Relative_RMSE"]),
          Standardized_RMSE = unname(val["Standardized_RMSE"]),
          Converged = as.numeric(converged)
        )

        fold_train[[fold]] <- read_train_predictions(pred_train)
      } else {
        y_bglr <- pheno$y_trait
        validation_rows <- pheno$Fold == fold
        y_bglr[validation_rows] <- NA_real_
        bglr_fit <- .gspip_fit_bglr_single(
          y = y_bglr,
          ids = pheno$GenoID,
          grm = grm,
          weights = pheno$w,
          n_iter = bglr_n_iter,
          burn_in = bglr_burn_in,
          thin = bglr_thin,
          seed = bglr_seed + (iteration - 1L) * run_folds + fold,
          R2 = bglr_R2,
          verbose = bglr_verbose
        )

        val_pred <- bglr_fit$prediction[validation_rows]
        cv_predictions[validation_ids, iteration] <- val_pred
        train_rows <- match(training_ids, pheno$GenoID)
        Dg <- mean(diag(training_grm)) - mean(training_grm)
        Ve <- bglr_fit$Ve_scale *
          mean(1 / training_data$w[!is.na(observed_y)])
        h2 <- safe_ratio(bglr_fit$Vg * Dg, bglr_fit$Vg * Dg + Ve)
        fold_h2[fold] <- h2
        val <- prediction_metrics(validation_data$y_trait, val_pred)

        metrics_list[[length(metrics_list) + 1L]] <- metric_rows(
          iteration, fold,
          Vg_posterior = bglr_fit$Vg,
          Residual_scale_posterior = bglr_fit$Ve_scale,
          Ve_posterior = Ve,
          PVE_posterior = h2,
          Heritability_posterior = h2,
          DIC = bglr_fit$DIC,
          MCMC_retained = bglr_fit$retained_samples,
          GOF = safe_cor(
            observed_y,
            bglr_fit$prediction[train_rows]
          ),
          PredictionAbility = unname(val["PA"]),
          Accuracy = accuracy_metric(unname(val["PA"]), h2),
          RMSE = unname(val["RMSE"]),
          Relative_RMSE = unname(val["Relative_RMSE"]),
          Standardized_RMSE = unname(val["Standardized_RMSE"])
        )

        fold_train[[fold]] <- data.frame(
          GenoID = training_ids,
          GEBV = bglr_fit$genomic[train_rows],
          SE = bglr_fit$prediction_sd[train_rows],
          stringsAsFactors = FALSE
        )
        intercept <- bglr_fit$intercept
      }
      fold_blue[[fold]] <- data.frame(
        Iter = iteration, Fold = fold, Intercept = intercept,
        Engine = engine, stringsAsFactors = FALSE
      )
    }
  
    # Pool the out-of-sample predictions of this iteration
    pooled_pred <- cv_predictions[pheno$GenoID, iteration]
    all_folds   <- prediction_metrics(pheno$y_trait, pooled_pred)
    metrics_list[[length(metrics_list) + 1L]] <- metric_rows(
      iteration, 0L,
      PA_CV_allfolds = unname(all_folds["PA"]),
      Accuracy_CV_allfolds = accuracy_metric(unname(all_folds["PA"]), mean_finite(fold_h2)),
      RMSE_CV_allfolds = unname(all_folds["RMSE"]),
      Standardized_RMSE_CV_allfolds = unname(all_folds["Standardized_RMSE"])
    )
  
    # Prediction ability within strata
    if (within_stratum_pa && (derive_strata_from_grm || !is.null(strata_col))) {
      for (s in unique(pheno$Stratum)) {
        in_s <- pheno$Stratum == s
        usable_s <- in_s & is.finite(pheno$y_trait) & is.finite(pooled_pred)
        if (sum(usable_s) >= min_stratum_n) {
          pa_s <- safe_cor(pheno$y_trait[in_s], pooled_pred[in_s])
          names(pa_s) <- paste0("PA_CV_within_", s)
          metrics_list[[length(metrics_list) + 1L]] <- metric_rows(iteration, 0L, pa_s)
        }
      }
    }
  
    # Average in-sample predictions over folds, and save progress
    train_all  <- do.call(rbind, fold_train)
    train_mean <- tapply(train_all$GEBV, train_all$GenoID, mean, na.rm = TRUE)
    se_mean    <- tapply(train_all$SE,   train_all$GenoID, mean, na.rm = TRUE)
    train_predictions[names(train_mean), iteration] <- train_mean
    train_se[names(se_mean), iteration] <- se_mean
    intercepts[[iteration]] <- do.call(rbind, fold_blue)
  
    data.table::fwrite(do.call(rbind, metrics_list), output_path("metrics.csv"))
  }
  
  #===============================================================================
  # STEP 5. EXPORT RESULTS
  #===============================================================================
  row_means <- function(mat) {
    m <- rowMeans(mat, na.rm = TRUE)
    m[is.nan(m)] <- NA_real_
    m
  }
  export_iter_table <- function(mat, suffix) {
    out <- data.frame(GenoID = rownames(mat), mat, check.names = FALSE, stringsAsFactors = FALSE)
    out$Mean <- row_means(mat)
    data.table::fwrite(out, output_path(suffix))
  }
  
  export_iter_table(cv_predictions, "cv_predictions.csv")
  export_iter_table(train_predictions, "train_predictions.csv")
  export_iter_table(train_se, "train_se.csv")
  data.table::fwrite(do.call(rbind, intercepts), output_path("intercepts.csv"))
  
  metrics <- do.call(rbind, metrics_list)
  metrics_summary <- summarise_metrics(metrics)
  data.table::fwrite(metrics_summary, output_path("summary.csv"))
  
  n_not_converged <- if (engine == "asreml") {
    sum(metrics$Metric == "Converged" & metrics$Value == 0)
  } else {
    0L
  }
  if (n_not_converged > 0L) {
    warning(n_not_converged, " fold model(s) did not report convergence. See 'Converged' in the metrics file.",
            call. = FALSE)
  }
  
  log_msg("Finished. Results are in: ", normalizePath(output_dir))
  print(metrics_summary)

  output_suffixes <- c("metrics.csv", "summary.csv", "cv_predictions.csv",
                       "train_predictions.csv", "train_se.csv", "intercepts.csv",
                       "derived_strata.csv", "derived_strata_selection.csv",
                       "derived_strata_pc_variance.csv", "split_ids.csv",
                       "run_info.txt")
  output_files <- stats::setNames(vapply(output_suffixes, output_path, character(1)),
                                  output_suffixes)
  output_files <- output_files[file.exists(output_files)]
  result <- list(engine = engine, metrics = metrics, summary = metrics_summary,
                 cv_predictions = cv_predictions,
                 train_predictions = train_predictions, train_se = train_se,
                 intercepts = intercepts, phenotype = pheno, fixed_split = fixed_split,
                 derived_strata = grm_strata_result, output_files = output_files)
  class(result) <- c("gspip_gblup_single", "list")
  invisible(result)
}

#' @export
print.gspip_gblup_single <- function(x, ...) {
  cat("GSPip single-environment GBLUP evaluation\n")
  cat("  Engine         :", if (is.null(x$engine)) "asreml" else x$engine, "\n")
  cat("  Prediction rows:", nrow(x$cv_predictions), "\n")
  cat("  Output files   :", length(x$output_files), "\n")
  print(x$summary)
  invisible(x)
}
