#' Run genotype-by-environment GBLUP evaluation
#'
#' Fits genomic prediction models from stage-1 BLUEs observed across
#' environments and evaluates CV1, CV2, or a fixed genotype-level split.
#' The ASReml backend uses REML; the open-source BGLR backend uses genomic-main
#' and genomic-by-environment RKHS kernels.
#'
#' @param grm_file Path to a square GRM CSV whose first column contains genotype IDs.
#' @param pheno_file Path to the genotype-by-environment BLUE table.
#' @param output_dir Directory receiving analysis outputs.
#' @param engine Explicit model engine: \code{"asreml"} or \code{"bglr"}.
#' @param run_label Prefix used for output files.
#' @param id_col,trait_col,env_col,strata_col,weight_col,se_col Input column names.
#' @param derive_strata_from_grm Derive CV strata from GRM PCs.
#' @param grm_strata_n_clusters Optional fixed cluster count; `NULL` selects automatically.
#' @param grm_strata_max_clusters,grm_strata_n_pcs,grm_strata_seed GRM clustering controls.
#' @param fix_stage1_error Fix the stage-1 residual scale to one.
#' @param nonmarker_effects Include non-marker genotype and G by E effects.
#' @param env_h2 Optional named heritability vector by environment.
#' @param cv_scheme Either `"CV1"` or `"CV2"`.
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
#' @return Invisibly, a `gspip_gblup_gxe` list containing metrics, predictions,
#'   input data used by the model, and output paths.
#' @export
run_gblup_gxe <- function(
    grm_file, pheno_file, output_dir, engine,
    run_label = "GBLUP_GxE_NBL",
    id_col = "GenoID", trait_col = "NBL", env_col = "Env",
    strata_col = "Population", weight_col = "Weight", se_col = "SE",
    derive_strata_from_grm = FALSE, grm_strata_n_clusters = NULL,
    grm_strata_max_clusters = 10L, grm_strata_n_pcs = 2L,
    grm_strata_seed = 123L, fix_stage1_error = TRUE,
    nonmarker_effects = TRUE, env_h2 = NULL, cv_scheme = c("CV1", "CV2"),
    k_folds = 5L, n_iterations = 50L, base_seed = 0L,
    evaluation_mode = c("repeated_cv", "fixed_split"),
    training_ids_file = NULL, target_ids_file = NULL,
    within_stratum_pa = TRUE, min_stratum_n = 10L,
    overwrite_outputs = FALSE, asreml_workspace = 32e07,
    bglr_n_iter = 12000L, bglr_burn_in = 2000L, bglr_thin = 10L,
    bglr_seed = 123L, bglr_R2 = 0.5, bglr_verbose = FALSE) {
  engine <- match.arg(engine, c("asreml", "bglr"))
  cv_scheme <- match.arg(cv_scheme)
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
    idv <- get("idv", envir = asreml_namespace, inherits = FALSE)
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
    suffixes <- c("metrics.csv", "summary.csv", "cv_predictions.csv", "cv_predictions_by_env.csv",
                  "train_predictions.csv", "train_se.csv", "intercepts.csv", "derived_strata.csv",
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
    if (is.null(weight_col)) {
      stop("gblup_cv_gxe.R does not support weight_col <- NULL. Supply a real weight column or ",
           "an SE column; see 'If the data has no Weight or SE column' in README.md.", call. = FALSE)
    }
    if (!all(vapply(list(id_col, trait_col, env_col, weight_col, se_col),
                    nonempty_string, logical(1))) || !nullable_string(strata_col)) {
      stop("Column settings must be non-empty names (strata_col may be NULL).", call. = FALSE)
    }
    required_columns <- c(id_col, trait_col, env_col, weight_col, se_col)
    optional_columns <- c(strata_col)
    if (anyDuplicated(c(required_columns, optional_columns))) {
      stop("Column settings must refer to distinct columns.", call. = FALSE)
    }
    if (!nonempty_string(cv_scheme) || !(cv_scheme %in% c("CV1", "CV2"))) {
      stop('cv_scheme must be "CV1" or "CV2".', call. = FALSE)
    }
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
      stop("The BGLR backend currently fits genomic main and genomic G x E ",
           "kernels only. Set nonmarker_effects <- FALSE.", call. = FALSE)
    }
    if (!nonempty_string(evaluation_mode) || !evaluation_mode %in% c("repeated_cv", "fixed_split")) {
      stop('evaluation_mode must be "repeated_cv" or "fixed_split".', call. = FALSE)
    }
    if (evaluation_mode == "fixed_split" &&
        (!nonempty_string(training_ids_file) || !nonempty_string(target_ids_file))) {
      stop("training_ids_file and target_ids_file are required for evaluation_mode <- \"fixed_split\".",
           call. = FALSE)
    }
    if (evaluation_mode == "fixed_split" && cv_scheme != "CV1") {
      stop("Fixed STPGA evaluation is genotype-level and requires cv_scheme <- \"CV1\".", call. = FALSE)
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
    observed <- !is.na(pheno[[trait_col]])
    has_w  <- weight_col %in% names(pheno)
    has_se <- se_col %in% names(pheno)
    if (!has_w && !has_se) {
      stop("The G x E phenotype file has neither the weight column ('", weight_col,
           "') nor the SE column ('", se_col, "'). This script requires at least one real ",
           "precision measure and will not create artificial unit weights. See README.md.", call. = FALSE)
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
  
    absent <- setdiff(c(id_col, trait_col, strata_col, env_col), names(pheno))
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
  
    if (anyNA(pheno[[env_col]])) stop("Missing values in environment column: ", env_col, call. = FALSE)
    env <- as.character(pheno[[env_col]])
    if (any(!nzchar(env))) stop("Empty values in environment column: ", env_col, call. = FALSE)
    if (length(unique(env)) < 2) {
      stop("Only one environment found. For one BLUE per genotype, use gblup_cv_single.R.", call. = FALSE)
    }
    no_data <- setdiff(env, env[!is.na(pheno[[trait_col]])])
    if (length(no_data) > 0) {
      stop("Environment(s) with no trait values: ", paste(no_data, collapse = ", "), call. = FALSE)
    }
    if (anyDuplicated(paste(ids, env, sep = "\r"))) {
      stop("Each genotype x environment must have one row. Duplicated genotype(s): ",
           show_ids(unique(ids[duplicated(paste(ids, env, sep = "\r"))])), call. = FALSE)
    }
  
    strata <- stratum_values(pheno, ids, grm_strata)
    n_strata <- tapply(strata, ids, function(x) length(unique(x)))
    if (any(n_strata > 1)) {
      stop("Each genotype must belong to one stratum. Inconsistent: ",
           show_ids(names(n_strata)[n_strata > 1]), call. = FALSE)
    }
    weights <- resolve_weights(pheno)
  
    n_genotypes <- length(unique(ids))
    if (n_genotypes < k) stop("k_folds exceeds the number of genotypes in the phenotype file.", call. = FALSE)
    if (sum(!is.na(pheno[[trait_col]])) < k) {
      stop("Fewer observed trait values than k_folds; each fold needs validation data.", call. = FALSE)
    }
  
    out <- data.frame(GenoID = ids, Env = env, y_trait = pheno[[trait_col]],
                      w = weights$w, Stratum = strata, stringsAsFactors = FALSE)
    geno_strata <- out$Stratum[!duplicated(out$GenoID)]
    small <- names(which(table(geno_strata) < k))
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
  
  check_env_h2 <- function(h2, envs) {
    if (is.null(h2)) return(invisible(NULL))
    if (!is.numeric(h2) || is.null(names(h2)) || anyNA(h2) || any(h2 <= 0 | h2 > 1) ||
        anyNA(names(h2)) || any(!nzchar(names(h2))) || anyDuplicated(names(h2))) {
      stop("env_h2 must be a named numeric vector with values above 0 and up to 1.", call. = FALSE)
    }
    absent <- setdiff(envs, names(h2))
    if (length(absent) > 0) stop("env_h2 has no value for: ", paste(absent, collapse = ", "), call. = FALSE)
    unknown <- setdiff(names(h2), envs)
    if (length(unknown) > 0) stop("env_h2 names not found in the data: ", paste(unknown, collapse = ", "), call. = FALSE)
    invisible(NULL)
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
  
  env_cors <- function(obs, pred, env) {
    vapply(split(seq_along(obs), env), function(i) safe_cor(obs[i], pred[i]), numeric(1))
  }
  mean_or_na <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  mean_env_cor <- function(obs, pred, env) mean_or_na(env_cors(obs, pred, env))
  
  # With env_h2: mean of r_j / sqrt(h2_j). Otherwise: mean r / sqrt(Cullis h2)
  accuracy_metric <- function(obs, pred, env, cullis_h2) {
    r <- env_cors(obs, pred, env)
    if (!is.null(env_h2)) return(mean_or_na(r / sqrt(env_h2[names(r)])))
    if (is.finite(cullis_h2) && cullis_h2 > 0) mean_or_na(r) / sqrt(cullis_h2) else NA_real_
  }
  
  prediction_metrics <- function(obs, pred, env) {
    ok <- complete.cases(obs, pred, env) & is.finite(obs) & is.finite(pred)
    if (!any(ok)) return(c(PA = NA_real_, RMSE = NA_real_, Relative_RMSE = NA_real_,
                           Standardized_RMSE = NA_real_))
    rmse <- sqrt(mean((obs[ok] - pred[ok])^2))
    obs_mean <- mean(obs[ok])
    centered_obs <- obs[ok] - ave(obs[ok], env[ok], FUN = mean)
    obs_sd <- stats::sd(centered_obs)
    c(PA = mean_env_cor(obs, pred, env),
      RMSE = rmse,
      Relative_RMSE = if (is.finite(obs_mean) && obs_mean != 0) rmse / abs(obs_mean) * 100 else NA_real_,
      Standardized_RMSE = if (is.finite(obs_sd) && obs_sd > 0) rmse / obs_sd else NA_real_)
  }
  
  safe_ratio <- function(numerator, denominator) {
    if (is.finite(numerator) && is.finite(denominator) && denominator > 0) {
      numerator / denominator
    } else {
      NA_real_
    }
  }
  
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
  
  # Coefficient names: <...>_<env>:<term>_<genotype>. Absent non-marker cells are 0
  read_env_term <- function(re, tag, ids, env_levels, absent_as_zero = FALSE) {
    nms <- rownames(re)
    full_tag <- paste0(":", tag)
    pos <- regexpr(full_tag, nms, fixed = TRUE)
    idx <- which(pos > 0)
    prefix <- substr(nms[idx], 1, pos[idx] - 1)
    geno <- substring(nms[idx], pos[idx] + nchar(full_tag))
    out <- matrix(NA_real_, length(ids), length(env_levels), dimnames = list(ids, env_levels))
    free <- rep(TRUE, length(idx))
    for (lev in env_levels[order(-nchar(env_levels))]) {  # Longest first: "E_1" before "1"
      hit <- free & endsWith(prefix, paste0("_", lev)) & geno %in% ids
      out[cbind(geno[hit], rep(lev, sum(hit)))] <- re[idx[hit], 1]
      free[hit] <- FALSE
    }
    if (absent_as_zero && any(!is.na(out))) out[is.na(out)] <- 0
    if (anyNA(out)) {
      stop("Could not read the '", tag, "' environment effects. Random coefficient names look like: ",
           show_ids(nms, 3), call. = FALSE)
    }
    out
  }
  
  extract_varcomp <- function(fit) {
    vc <- as.data.frame(summary(fit)$varcomp)
    terms <- rownames(vc)
    has_vm  <- grepl("vm(GenoID", terms, fixed = TRUE)
    has_ide <- grepl("ide(GenoID", terms, fixed = TRUE)
    has_env <- grepl("Env", terms, fixed = TRUE)
    get_one <- function(keep, label) {
      v <- vc[keep, 1]
      if (length(v) != 1) stop("Could not find ", label, " in: ", paste(terms, collapse = ", "), call. = FALSE)
      v
    }
    list(
      Vg  = get_one(has_vm & !has_env, "Vg"),
      Vg_pos = which(has_vm & !has_env),
      Vge = get_one(has_vm & has_env, "Vge"),
      Vr  = if (nonmarker_effects) get_one(has_ide & !has_env, "Vr") else 0,
      Vre = if (nonmarker_effects) get_one(has_ide & has_env, "Vre") else 0,
      s   = if (fix_stage1_error) 1 else get_one(!has_vm & !has_ide, "the residual scale")
    )
  }
  
  extract_fixed <- function(fit, env_levels) {
    fx <- as.matrix(fit$coefficients$fixed)
    mu <- fx[rownames(fx) == "(Intercept)", 1]
    if (length(mu) != 1) stop("Could not find the intercept in the fixed effects.", call. = FALSE)
    env_eff <- setNames(fx[match(paste0("Env_", env_levels), rownames(fx)), 1], env_levels)
    if (anyNA(env_eff)) stop("Could not match environment effects. Found: ", show_ids(rownames(fx)), call. = FALSE)
    list(mu = unname(mu), env = env_eff)
  }
  
  extract_random <- function(fit, ids, env_levels) {
    re <- as.matrix(fit$coefficients$random)
    vm_tag  <- "vm(GenoID, training_grm_inv.sparse)_"
    ide_tag <- "ide(GenoID)_"
    list(
      u  = read_main_term(re, vm_tag, ids),
      ge = read_env_term(re, vm_tag, ids, env_levels),
      r  = if (nonmarker_effects) read_main_term(re, ide_tag, ids) else setNames(rep(0, length(ids)), ids),
      re = if (nonmarker_effects) {
        read_env_term(re, ide_tag, ids, env_levels, absent_as_zero = TRUE)
      } else {
        matrix(0, length(ids), length(env_levels), dimnames = list(ids, env_levels))
      }
    )
  }
  
  # Genotypes absent from training: g and gE from G_new,train %*% G_train^-1; r and rE = 0
  extend_to_new_genotypes <- function(effects, new_ids, train_ids) {
    if (length(new_ids) == 0) return(effects)
    projected <- grm[new_ids, train_ids, drop = FALSE] %*%
      solve(grm[train_ids, train_ids, drop = FALSE], cbind(effects$u, effects$ge))
    envs <- colnames(effects$ge)
    zero_re <- matrix(0, length(new_ids), length(envs), dimnames = list(new_ids, envs))
    list(
      u  = c(effects$u, setNames(projected[, 1], new_ids)),
      ge = rbind(effects$ge, matrix(projected[, -1], nrow = length(new_ids), dimnames = list(new_ids, envs))),
      r  = c(effects$r, setNames(rep(0, length(new_ids)), new_ids)),
      re = rbind(effects$re, zero_re)
    )
  }
  
  # mu + E_j + g_i + r_i + gE_ij + rE_ij. An environment with no training data gives NA
  predict_rows <- function(data, fixed, effects) {
    out <- fixed$mu + unname(fixed$env[data$Env]) + unname(effects$u[data$GenoID]) + unname(effects$r[data$GenoID])
    known <- data$Env %in% colnames(effects$ge)
    cells <- cbind(data$GenoID[known], data$Env[known])
    out[known] <- out[known] + effects$ge[cells] + effects$re[cells]
    out
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
      paste("Trait:", trait_col, "| Stratum:", if (is.null(strata_col)) "none" else strata_col,
            "| Environment:", env_col),
      paste("Stratum source:", stratum_source),
      if (!is.null(grm_strata_result)) {
        sprintf("GRM clustering: %d clusters (%s; average silhouette %.3f) using %d PCs (cumulative variance %.1f%%), seed %s",
                grm_strata_result$selected_k, grm_strata_result$selection_mode,
                grm_strata_result$selected_silhouette, grm_strata_result$n_pcs,
                grm_strata_result$cumulative_pc_variance, grm_strata_seed)
      },
      paste("CV scheme:", cv_scheme),
      paste("Stage-1 error:", if (fix_stage1_error) "fixed, residual = 1/w" else
            paste0(if (engine == "bglr") "posterior " else "",
                   "scale estimated, residual = s/w")),
      paste("Non-marker effects r + rE:", if (nonmarker_effects) "yes" else "no"),
      paste("Accuracy uses:", if (is.null(env_h2)) {
        if (engine == "asreml") "Cullis genotype-mean heritability, per fold" else
          "BGLR posterior genotype-mean heritability approximation, per fold"
      } else {
        paste0("env_h2: ", paste(names(env_h2), env_h2, sep = "=", collapse = ", "))
      }),
      paste("Weights:", pheno_input$weight_source),
      paste("Weight range:", paste(signif(range(pheno$w[!is.na(pheno$y_trait)]), 4), collapse = " to ")),
      paste("Folds:", run_folds, "| Iterations:", run_iterations,
            if (evaluation_mode == "repeated_cv") paste("| Base seed:", base_seed) else ""),
      paste("Genotypes in GRM:", nrow(grm)),
      paste("Genotypes in the phenotype file:", nrow(geno_table)),
      paste("Rows with a trait value:", sum(!is.na(pheno$y_trait))),
      paste("Rows without a trait value (predicted, not in metrics):", sum(is.na(pheno$y_trait))),
      paste("In GRM without a phenotype row (NA in prediction tables):", length(setdiff(all_genoids, pheno$GenoID))),
      "", "Stratum sizes (genotypes):",
      capture.output(print(table(geno_table$Stratum))),
      "", "Rows per environment:", capture.output(print(table(pheno$Env))),
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
    fixed_split <- read_fixed_split(training_ids_file, target_ids_file, unique(pheno$GenoID))
    keep_ids <- c(fixed_split$training, fixed_split$target)
    pheno <- pheno[pheno$GenoID %in% keep_ids, , drop = FALSE]
    observed_training_ids <- unique(pheno$GenoID[pheno$GenoID %in% fixed_split$training &
                                                  !is.na(pheno$y_trait)])
    no_training_values <- setdiff(fixed_split$training, observed_training_ids)
    if (length(no_training_values)) {
      stop("Selected training genotype(s) have no observed trait in any environment: ",
           show_ids(no_training_values), call. = FALSE)
    }
    data.table::fwrite(data.frame(
      GenoID = c(fixed_split$training, fixed_split$target),
      Role = c(rep("Training", length(fixed_split$training)), rep("Target", length(fixed_split$target))),
      stringsAsFactors = FALSE
    ), output_path("split_ids.csv"))
  }
  run_iterations <- if (evaluation_mode == "fixed_split") 1L else n_iterations
  run_folds <- if (evaluation_mode == "fixed_split") 1L else k_folds
  geno_table  <- pheno[!duplicated(pheno$GenoID), c("GenoID", "Stratum")]
  check_env_h2(env_h2, unique(pheno$Env))
  
  if (!is.null(grm_strata_result)) {
    data.table::fwrite(grm_strata_result$table, output_path("derived_strata.csv"))
    data.table::fwrite(grm_strata_result$selection, output_path("derived_strata_selection.csv"))
    data.table::fwrite(grm_strata_result$pc_variance, output_path("derived_strata_pc_variance.csv"))
  }
  write_run_info()
  if (evaluation_mode == "fixed_split") {
    log_msg(sprintf("Fixed genotype split: %d training and %d target genotypes; %d observed rows; %s; strata: %s",
                    length(fixed_split$training), length(fixed_split$target), sum(!is.na(pheno$y_trait)),
                    pheno_input$weight_source, stratum_source))
  } else {
    log_msg(sprintf("%s: %d genotypes, %d rows with a trait value, %d to predict only; %s; strata: %s; %d iterations of %d-fold CV",
                    cv_scheme, nrow(geno_table), sum(!is.na(pheno$y_trait)), sum(is.na(pheno$y_trait)),
                    pheno_input$weight_source, stratum_source, run_iterations, run_folds))
  }
  
  #===============================================================================
  # STEP 2. DEFINE THE MODEL
  #===============================================================================
  # Top-level objects: predict() re-evaluates the model call
  fixed_formula <- if (engine == "asreml") y_trait ~ Env else NULL
  random_formula <- if (engine == "asreml") {
    if (nonmarker_effects) {
      ~ vm(GenoID, training_grm_inv.sparse) +
        idv(Env):vm(GenoID, training_grm_inv.sparse) +
        ide(GenoID) + idv(Env):ide(GenoID)
    } else {
      ~ vm(GenoID, training_grm_inv.sparse) +
        idv(Env):vm(GenoID, training_grm_inv.sparse)
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
  cv_row_predictions <- empty_iter_table(seq_len(nrow(pheno)))
  
  #===============================================================================
  # STEP 4. RUN THE REQUESTED EVALUATION
  #===============================================================================
  # asreml::asreml() and predict() stay at top level: predict() needs the same objects
  for (iteration in seq_len(run_iterations)) {
    log_msg(sprintf("Iteration %d of %d", iteration, run_iterations))
  
    # Assign folds: CV1 to genotypes, CV2 to BLUEs within environment x stratum
    if (evaluation_mode == "fixed_split") {
      pheno$Fold <- ifelse(pheno$GenoID %in% fixed_split$target, 1L, 0L)
    } else {
      set.seed(base_seed + iteration)
      if (cv_scheme == "CV1") {
        geno_table$Fold <- make_folds(geno_table$Stratum, run_folds)
        pheno$Fold <- geno_table$Fold[match(pheno$GenoID, geno_table$GenoID)]
      } else {
        pheno$Fold <- make_folds(paste(pheno$Env, pheno$Stratum, sep = "\r"), run_folds)
      }
    }
  
    observed_by_fold <- tabulate(pheno$Fold[!is.na(pheno$y_trait)], nbins = run_folds)
    if (any(observed_by_fold == 0L)) {
      warning("Iteration ", iteration, ": fold(s) without observed validation traits: ",
              paste(which(observed_by_fold == 0L), collapse = ", "), call. = FALSE)
    }
  
    fold_train <- vector("list", run_folds)
    fold_h2    <- numeric(run_folds)
    fold_blue  <- vector("list", run_folds)
    fold_gebv  <- vector("list", run_folds)
    n_untested <- 0
  
    for (fold in seq_len(run_folds)) {
      # Split training and validation sets
      val_rows        <- which(pheno$Fold == fold)
      training_data   <- pheno[pheno$Fold != fold, ]
      validation_data <- pheno[val_rows, ]
      env_levels      <- sort(unique(training_data$Env[!is.na(training_data$y_trait)]))
      training_data   <- training_data[training_data$Env %in% env_levels, ]
      training_ids    <- unique(training_data$GenoID)
      validation_ids  <- unique(validation_data$GenoID)
      new_ids         <- setdiff(validation_ids, training_ids)
  
      if (sum(!is.na(training_data$y_trait)) < 2L) {
        stop(sprintf("Iteration %d, fold %d has fewer than two observed training values. Check the split or add data.",
                     iteration, fold), call. = FALSE)
      }
      if (length(env_levels) < 2L) {
        stop(sprintf("Iteration %d, fold %d has observed training data in fewer than two environments. Check the split or add data.",
                     iteration, fold), call. = FALSE)
      }
  
      unseen_env <- setdiff(validation_data$Env, env_levels)
      if (length(unseen_env) > 0) {
        warning(sprintf("Iteration %d, fold %d: environment(s) with no training data, not predicted: %s",
                        iteration, fold, paste(unseen_env, collapse = ", ")), call. = FALSE)
      }
  
      training_grm <- grm[training_ids, training_ids, drop = FALSE]
      observed_train <- training_data
      if (engine == "asreml") {
        training_grm_inv.sparse <- ASRgenomics::G.inverse(
          G = training_grm, sparseform = TRUE
        )$Ginv
        training_data$GenoID <- factor(training_data$GenoID, levels = training_ids)
        training_data$Env    <- factor(training_data$Env, levels = env_levels)

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

        fixed_coef  <- extract_fixed(GBLUPcv, env_levels)
        effects     <- extract_random(GBLUPcv, training_ids, env_levels)
        effects_all <- extend_to_new_genotypes(effects, new_ids, training_ids)

        val_row_pred <- predict_rows(validation_data, fixed_coef, effects_all)
        cv_row_predictions[val_rows, iteration] <- val_row_pred
        fold_gebv[[fold]] <- data.frame(
          GenoID = validation_ids,
          GEBV = fixed_coef$mu + mean(fixed_coef$env) +
            unname(effects_all$u[validation_ids]),
          stringsAsFactors = FALSE
        )

        tested_ids <- unique(
          observed_train$GenoID[!is.na(observed_train$y_trait)]
        )
        fold_untested <- sum(
          !is.na(validation_data$y_trait) &
            !validation_data$GenoID %in% tested_ids
        )
        n_untested <- n_untested + fold_untested

        Vg  <- varcomp$Vg
        Vge <- varcomp$Vge
        Vr  <- varcomp$Vr
        Vre <- varcomp$Vre
        Ve  <- varcomp$s *
          mean(1 / training_data$w[!is.na(training_data$y_trait)])
        Dg  <- heritability$Dg
        train_row_pred <- predict_rows(observed_train, fixed_coef, effects)
        val <- prediction_metrics(
          validation_data$y_trait, val_row_pred, validation_data$Env
        )

        metrics_list[[length(metrics_list) + 1L]] <- metric_rows(
          iteration, fold,
          Vg = Vg,
          Vge = Vge,
          Vr = Vr,
          Vre = Vre,
          Residual_scale = varcomp$s,
          Ve = Ve,
          PVE = safe_ratio(Vg * Dg, (Vg + Vge) * Dg + Vr + Vre + Ve),
          heritability$metrics,
          GOF = mean_env_cor(
            observed_train$y_trait, train_row_pred, observed_train$Env
          ),
          PredictionAbility = unname(val["PA"]),
          Accuracy = accuracy_metric(
            validation_data$y_trait, val_row_pred, validation_data$Env,
            heritability$h2
          ),
          RMSE = unname(val["RMSE"]),
          Relative_RMSE = unname(val["Relative_RMSE"]),
          Standardized_RMSE = unname(val["Standardized_RMSE"]),
          Untested_val_rows = fold_untested,
          Converged = as.numeric(converged)
        )

        fold_train[[fold]] <- read_train_predictions(pred_train)
        intercept <- fixed_coef$mu
      } else {
        model_rows <- which(pheno$Env %in% env_levels)
        validation_model <- pheno$Fold[model_rows] == fold
        y_bglr <- pheno$y_trait[model_rows]
        y_bglr[validation_model] <- NA_real_
        bglr_fit <- .gspip_fit_bglr_gxe(
          y = y_bglr,
          ids = pheno$GenoID[model_rows],
          environments = pheno$Env[model_rows],
          grm = grm,
          weights = pheno$w[model_rows],
          n_iter = bglr_n_iter,
          burn_in = bglr_burn_in,
          thin = bglr_thin,
          seed = bglr_seed + (iteration - 1L) * run_folds + fold,
          R2 = bglr_R2,
          verbose = bglr_verbose
        )

        val_model_position <- match(val_rows, model_rows)
        val_row_pred <- rep(NA_real_, length(val_rows))
        predictable <- !is.na(val_model_position)
        val_row_pred[predictable] <-
          bglr_fit$prediction[val_model_position[predictable]]
        cv_row_predictions[val_rows, iteration] <- val_row_pred

        main_by_id <- tapply(
          bglr_fit$genomic_main, pheno$GenoID[model_rows], mean
        )
        fixed_by_env <- tapply(
          bglr_fit$environment_fixed, pheno$Env[model_rows], mean
        )
        fold_gebv[[fold]] <- data.frame(
          GenoID = validation_ids,
          GEBV = bglr_fit$intercept + mean(fixed_by_env) +
            unname(main_by_id[validation_ids]),
          stringsAsFactors = FALSE
        )

        tested_ids <- unique(
          observed_train$GenoID[!is.na(observed_train$y_trait)]
        )
        fold_untested <- sum(
          !is.na(validation_data$y_trait) &
            !validation_data$GenoID %in% tested_ids
        )
        n_untested <- n_untested + fold_untested

        Dg <- mean(diag(training_grm)) - mean(training_grm)
        Ve <- bglr_fit$Ve_scale *
          mean(1 / training_data$w[!is.na(training_data$y_trait)])
        genomic_main_variance <- bglr_fit$Vg * Dg
        genomic_gxe_variance <- bglr_fit$Vge * Dg
        posterior_pve <- safe_ratio(
          genomic_main_variance,
          genomic_main_variance + genomic_gxe_variance + Ve
        )
        posterior_h2_mean <- safe_ratio(
          genomic_main_variance,
          genomic_main_variance +
            genomic_gxe_variance / length(env_levels) +
            Ve / length(env_levels)
        )
        fold_h2[fold] <- posterior_h2_mean

        training_global_rows <- which(
          pheno$Fold != fold & pheno$Env %in% env_levels
        )
        training_model_position <- match(training_global_rows, model_rows)
        train_row_pred <- bglr_fit$prediction[training_model_position]
        val <- prediction_metrics(
          validation_data$y_trait, val_row_pred, validation_data$Env
        )

        metrics_list[[length(metrics_list) + 1L]] <- metric_rows(
          iteration, fold,
          Vg_posterior = bglr_fit$Vg,
          Vge_posterior = bglr_fit$Vge,
          Residual_scale_posterior = bglr_fit$Ve_scale,
          Ve_posterior = Ve,
          PVE_posterior = posterior_pve,
          Heritability_genotype_mean_posterior = posterior_h2_mean,
          DIC = bglr_fit$DIC,
          MCMC_retained = bglr_fit$retained_samples,
          GOF = mean_env_cor(
            observed_train$y_trait, train_row_pred, observed_train$Env
          ),
          PredictionAbility = unname(val["PA"]),
          Accuracy = accuracy_metric(
            validation_data$y_trait, val_row_pred, validation_data$Env,
            posterior_h2_mean
          ),
          RMSE = unname(val["RMSE"]),
          Relative_RMSE = unname(val["Relative_RMSE"]),
          Standardized_RMSE = unname(val["Standardized_RMSE"]),
          Untested_val_rows = fold_untested
        )

        train_id_rows <- pheno$GenoID[training_global_rows]
        train_main <- tapply(
          bglr_fit$genomic_main[training_model_position], train_id_rows, mean
        )
        train_sd <- tapply(
          bglr_fit$prediction_sd[training_model_position], train_id_rows, mean,
          na.rm = TRUE
        )
        train_sd[is.nan(train_sd)] <- NA_real_
        fold_train[[fold]] <- data.frame(
          GenoID = names(train_main),
          GEBV = as.numeric(train_main),
          SE = as.numeric(train_sd[names(train_main)]),
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
    gebv_all  <- do.call(rbind, fold_gebv)
    gebv_mean <- tapply(gebv_all$GEBV, gebv_all$GenoID, mean, na.rm = TRUE)
    cv_predictions[names(gebv_mean), iteration] <- gebv_mean
  
    pooled_pred <- cv_row_predictions[, iteration]
    all_folds   <- prediction_metrics(pheno$y_trait, pooled_pred, pheno$Env)
    metrics_list[[length(metrics_list) + 1L]] <- metric_rows(
      iteration, 0L,
      PA_CV_allfolds = unname(all_folds["PA"]),
      Accuracy_CV_allfolds = accuracy_metric(pheno$y_trait, pooled_pred, pheno$Env, mean_or_na(fold_h2)),
      RMSE_CV_allfolds = unname(all_folds["RMSE"]),
      Standardized_RMSE_CV_allfolds = unname(all_folds["Standardized_RMSE"]),
      Untested_val_rows_CV_allfolds = n_untested
    )
  
    # Prediction ability and accuracy within environments
    for (e in unique(pheno$Env)) {
      in_e  <- pheno$Env == e
      pa_e  <- safe_cor(pheno$y_trait[in_e], pooled_pred[in_e])
      acc_e <- accuracy_metric(pheno$y_trait[in_e], pooled_pred[in_e], pheno$Env[in_e], mean_or_na(fold_h2))
      names(pa_e)  <- paste0("PA_CV_env_", e)
      names(acc_e) <- paste0("Accuracy_CV_env_", e)
      metrics_list[[length(metrics_list) + 1L]] <- metric_rows(iteration, 0L, pa_e, acc_e)
    }
  
    # Prediction ability within strata
    if (within_stratum_pa && (derive_strata_from_grm || !is.null(strata_col))) {
      for (s in unique(geno_table$Stratum)) {
        in_s <- pheno$Stratum == s
        usable_s <- in_s & is.finite(pheno$y_trait) & is.finite(pooled_pred)
        if (length(unique(pheno$GenoID[usable_s])) >= min_stratum_n) {
          pa_s <- mean_env_cor(pheno$y_trait[in_s], pooled_pred[in_s], pheno$Env[in_s])
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
  
  by_env <- data.frame(GenoID = pheno$GenoID, Env = pheno$Env, Observed = pheno$y_trait,
                       cv_row_predictions, check.names = FALSE, stringsAsFactors = FALSE)
  by_env$Mean <- row_means(cv_row_predictions)
  data.table::fwrite(by_env, output_path("cv_predictions_by_env.csv"))
  
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
                       "cv_predictions_by_env.csv", "train_predictions.csv",
                       "train_se.csv", "intercepts.csv", "derived_strata.csv",
                       "derived_strata_selection.csv",
                       "derived_strata_pc_variance.csv", "split_ids.csv",
                       "run_info.txt")
  output_files <- stats::setNames(vapply(output_suffixes, output_path, character(1)),
                                  output_suffixes)
  output_files <- output_files[file.exists(output_files)]
  result <- list(engine = engine, metrics = metrics, summary = metrics_summary,
                 cv_predictions = cv_predictions,
                 predictions_by_environment = by_env,
                 train_predictions = train_predictions, train_se = train_se,
                 intercepts = intercepts, phenotype = pheno, fixed_split = fixed_split,
                 derived_strata = grm_strata_result, output_files = output_files)
  class(result) <- c("gspip_gblup_gxe", "list")
  invisible(result)
}

#' @export
print.gspip_gblup_gxe <- function(x, ...) {
  cat("GSPip genotype-by-environment GBLUP evaluation\n")
  cat("  Engine         :", if (is.null(x$engine)) "asreml" else x$engine, "\n")
  cat("  Prediction rows:", nrow(x$predictions_by_environment), "\n")
  cat("  Output files   :", length(x$output_files), "\n")
  print(x$summary)
  invisible(x)
}
