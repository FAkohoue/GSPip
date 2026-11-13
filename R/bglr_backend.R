# Internal BGLR backend -----------------------------------------------------------

.gspip_validate_bglr_controls <- function(n_iter, burn_in, thin, seed, R2, verbose) {
  whole_number <- function(x, lower = 0) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= lower && x == floor(x)
  }
  if (!whole_number(n_iter, 2)) {
    stop("bglr_n_iter must be a whole number of at least 2.", call. = FALSE)
  }
  if (!whole_number(burn_in, 0) || burn_in >= n_iter) {
    stop("bglr_burn_in must be a non-negative whole number smaller than bglr_n_iter.",
         call. = FALSE)
  }
  if (!whole_number(thin, 1)) {
    stop("bglr_thin must be a positive whole number.", call. = FALSE)
  }
  retained <- floor(n_iter / thin) - floor(burn_in / thin)
  if (retained < 100L) {
    stop("BGLR must retain at least 100 posterior samples. Increase bglr_n_iter ",
         "or reduce bglr_burn_in/bglr_thin.", call. = FALSE)
  }
  if (!whole_number(seed, 0) || seed > .Machine$integer.max) {
    stop("bglr_seed must be a non-negative whole number within R's integer range.",
         call. = FALSE)
  }
  if (!is.numeric(R2) || length(R2) != 1L || !is.finite(R2) ||
      R2 <= 0 || R2 >= 1) {
    stop("bglr_R2 must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("bglr_verbose must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(as.integer(retained))
}

.gspip_bglr_component <- function(component, label, bglr_weights) {
  if (is.null(component) || is.null(component$u) || is.null(component$varU)) {
    stop("BGLR did not return the expected ", label,
         " RKHS component. Check the BGLR version and model inputs.", call. = FALSE)
  }
  u <- as.numeric(component$u)
  var_u <- as.numeric(component$varU)[1L]
  if (length(u) != length(bglr_weights) ||
      any(!is.finite(u)) || length(var_u) != 1L || !is.finite(var_u) || var_u < 0) {
    stop("BGLR returned non-finite values for the ", label, " component.",
         call. = FALSE)
  }
  # With BGLR weights, RKHS K is internally transformed to T K T and the
  # returned row effects are T u. Convert them back to the response scale.
  list(u = u / bglr_weights, varU = var_u)
}

.gspip_bglr_fixed_coefficients <- function(component) {
  if (is.null(component)) return(numeric(0))
  value <- component$b
  if (is.null(value)) value <- component$beta
  if (is.null(value)) {
    stop("BGLR did not return coefficients for the fixed environment effects.",
         call. = FALSE)
  }
  as.numeric(value)
}

.gspip_call_bglr <- function(y, ETA, weights, n_iter, burn_in, thin,
                             seed, R2, verbose) {
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    stop("engine = \"bglr\" requires the open-source BGLR package. Install it ",
         "with install.packages(\"BGLR\").", call. = FALSE)
  }
  if (!is.numeric(y) || length(y) < 2L ||
      any(!is.finite(y[!is.na(y)])) || sum(!is.na(y)) < 2L) {
    stop("BGLR needs at least two finite observed response values.", call. = FALSE)
  }
  if (!is.numeric(weights) || length(weights) != length(y) ||
      any(!is.finite(weights) | weights <= 0)) {
    stop("BGLR observation weights must be positive, finite, and aligned with y.",
         call. = FALSE)
  }

  retained <- .gspip_validate_bglr_controls(
    n_iter, burn_in, thin, seed, R2, verbose
  )
  bglr_weights <- sqrt(weights)
  set.seed(as.integer(seed))
  save_prefix <- paste0(tempfile("gspip_bglr_"), "_")
  on.exit(
    unlink(Sys.glob(paste0(save_prefix, "*")), recursive = FALSE, force = TRUE),
    add = TRUE
  )

  fit <- BGLR::BGLR(
    y = y,
    ETA = ETA,
    response_type = "gaussian",
    nIter = as.integer(n_iter),
    burnIn = as.integer(burn_in),
    thin = as.integer(thin),
    R2 = R2,
    # GSPip w is inverse residual variance. BGLR uses inverse squared weights.
    weights = bglr_weights,
    saveAt = save_prefix,
    rmExistingFiles = TRUE,
    verbose = verbose
  )

  if (is.null(fit$yHat) || length(fit$yHat) != length(y) ||
      any(!is.finite(fit$yHat))) {
    stop("BGLR did not return finite posterior mean predictions.", call. = FALSE)
  }
  if (is.null(fit$varE) || length(fit$varE) != 1L ||
      !is.finite(fit$varE) || fit$varE < 0) {
    stop("BGLR did not return a finite non-negative posterior residual variance.",
         call. = FALSE)
  }

  # BGLR stores yHat and SD.yHat on the internally weighted response scale.
  # Keep the fitted object intact for component extraction, and expose the
  # conversion vector explicitly so every row-level result is returned to the
  # original response scale exactly once.
  list(fit = fit, retained = retained, bglr_weights = bglr_weights)
}

.gspip_fit_bglr_single <- function(y, ids, grm, weights, n_iter, burn_in,
                                    thin, seed, R2, verbose) {
  K <- unname(grm[ids, ids, drop = FALSE])
  called <- .gspip_call_bglr(
    y = y,
    ETA = list(genomic = list(K = K, model = "RKHS")),
    weights = weights,
    n_iter = n_iter,
    burn_in = burn_in,
    thin = thin,
    seed = seed,
    R2 = R2,
    verbose = verbose
  )
  fit <- called$fit
  genomic <- .gspip_bglr_component(
    fit$ETA[["genomic"]], "genomic", called$bglr_weights
  )
  prediction_sd <- fit$SD.yHat
  if (is.null(prediction_sd) || length(prediction_sd) != length(y)) {
    prediction_sd <- rep(NA_real_, length(y))
  } else {
    prediction_sd <- as.numeric(prediction_sd) / called$bglr_weights
  }

  list(
    prediction = as.numeric(fit$yHat) / called$bglr_weights,
    prediction_sd = as.numeric(prediction_sd),
    genomic = genomic$u,
    Vg = genomic$varU,
    Ve_scale = as.numeric(fit$varE),
    intercept = as.numeric(fit$mu)[1L],
    retained_samples = called$retained,
    DIC = if (is.null(fit$fit$DIC)) NA_real_ else as.numeric(fit$fit$DIC)[1L]
  )
}

.gspip_fit_bglr_gxe <- function(y, ids, environments, grm, weights,
                                 n_iter, burn_in, thin, seed, R2, verbose) {
  env_factor <- factor(environments)
  if (nlevels(env_factor) < 2L) {
    stop("The BGLR G x E model needs at least two environments.", call. = FALSE)
  }
  X_env <- stats::model.matrix(~ env_factor)
  X_env <- X_env[, colnames(X_env) != "(Intercept)", drop = FALSE]
  K_main <- unname(grm[ids, ids, drop = FALSE])
  same_environment <- outer(as.character(env_factor), as.character(env_factor), "==")
  K_gxe <- K_main * same_environment

  called <- .gspip_call_bglr(
    y = y,
    ETA = list(
      environment = list(X = X_env, model = "FIXED"),
      genomic_main = list(K = K_main, model = "RKHS"),
      genomic_gxe = list(K = K_gxe, model = "RKHS")
    ),
    weights = weights,
    n_iter = n_iter,
    burn_in = burn_in,
    thin = thin,
    seed = seed,
    R2 = R2,
    verbose = verbose
  )
  fit <- called$fit
  main <- .gspip_bglr_component(
    fit$ETA[["genomic_main"]], "genomic main", called$bglr_weights
  )
  gxe <- .gspip_bglr_component(
    fit$ETA[["genomic_gxe"]], "genomic G x E", called$bglr_weights
  )
  env_coef <- .gspip_bglr_fixed_coefficients(fit$ETA[["environment"]])
  if (length(env_coef) != ncol(X_env) || any(!is.finite(env_coef))) {
    stop("BGLR returned invalid fixed environment coefficients.", call. = FALSE)
  }
  env_fixed <- as.vector(X_env %*% env_coef)
  prediction_sd <- fit$SD.yHat
  if (is.null(prediction_sd) || length(prediction_sd) != length(y)) {
    prediction_sd <- rep(NA_real_, length(y))
  } else {
    prediction_sd <- as.numeric(prediction_sd) / called$bglr_weights
  }

  list(
    prediction = as.numeric(fit$yHat) / called$bglr_weights,
    prediction_sd = as.numeric(prediction_sd),
    genomic_main = main$u,
    genomic_gxe = gxe$u,
    environment_fixed = env_fixed,
    Vg = main$varU,
    Vge = gxe$varU,
    Ve_scale = as.numeric(fit$varE),
    intercept = as.numeric(fit$mu)[1L],
    retained_samples = called$retained,
    DIC = if (is.null(fit$fit$DIC)) NA_real_ else as.numeric(fit$fit$DIC)[1L]
  )
}
