# Locate or copy an installed GSPip workflow script

GSPip provides the programmatic APIs \[build_grm()\],
\[optimize_training_population()\], \[run_gblup_single()\], and
\[run_gblup_gxe()\]. Editable scripts under \`inst/scripts\` are thin
wrappers around those functions. This helper locates an installed
wrapper or copies it to a user-controlled location.

## Usage

``` r
gspip_script(
  workflow = c("grm_builder", "optimize_training_population", "gblup_cv_single",
    "gblup_cv_gxe"),
  copy_to = NULL,
  overwrite = FALSE
)
```

## Arguments

- workflow:

  One of \`"grm_builder"\`, \`"optimize_training_population"\`,
  \`"gblup_cv_single"\`, or \`"gblup_cv_gxe"\`.

- copy_to:

  Optional destination file or directory. When \`NULL\`, return the
  installed template path without copying it.

- overwrite:

  Allow replacement of an existing destination file.

## Value

The normalized template or copied-file path, invisibly when copied.
