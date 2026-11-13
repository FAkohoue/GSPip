#' Locate or copy an installed GSPip workflow script
#'
#' GSPip provides the programmatic APIs [build_grm()],
#' [optimize_training_population()], [run_gblup_single()], and
#' [run_gblup_gxe()]. Editable scripts under `inst/scripts` are thin wrappers
#' around those functions. This helper locates an installed wrapper or copies
#' it to a user-controlled location.
#'
#' @param workflow One of `"grm_builder"`, `"optimize_training_population"`,
#'   `"gblup_cv_single"`, or `"gblup_cv_gxe"`.
#' @param copy_to Optional destination file or directory. When `NULL`, return
#'   the installed template path without copying it.
#' @param overwrite Allow replacement of an existing destination file.
#'
#' @return The normalized template or copied-file path, invisibly when copied.
#' @export
gspip_script <- function(
    workflow = c("grm_builder", "optimize_training_population",
                 "gblup_cv_single", "gblup_cv_gxe"),
    copy_to = NULL,
    overwrite = FALSE
) {
  workflow <- match.arg(workflow)
  script_name <- paste0(workflow, ".R")
  source_path <- system.file("scripts", script_name, package = "GSPip")
  if (!nzchar(source_path) || !file.exists(source_path)) {
    stop("Installed workflow template is missing: ", script_name, call. = FALSE)
  }
  if (is.null(copy_to)) return(normalizePath(source_path))
  if (!is.character(copy_to) || length(copy_to) != 1L || is.na(copy_to) || !nzchar(copy_to)) {
    stop("copy_to must be NULL or one non-empty path.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.", call. = FALSE)
  }
  destination <- if (dir.exists(copy_to)) file.path(copy_to, script_name) else copy_to
  if (file.exists(destination) && !overwrite) {
    stop("Destination exists; choose another path or set overwrite = TRUE: ",
         destination, call. = FALSE)
  }
  destination_dir <- dirname(destination)
  if (!dir.exists(destination_dir) &&
      !dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create destination directory: ", destination_dir, call. = FALSE)
  }
  if (!file.copy(source_path, destination, overwrite = overwrite, copy.mode = TRUE)) {
    stop("Could not copy workflow template to: ", destination, call. = FALSE)
  }
  invisible(normalizePath(destination))
}

