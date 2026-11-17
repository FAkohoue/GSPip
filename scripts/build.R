######################################################################################################
#######################################################################################################

# =============================================================================
# GSPip - clean build script
# Run SESSION A top-to-bottom, let it restart R, then run SESSION B.
# Both sessions must start with the working directory set to the package root
# (the folder containing DESCRIPTION).
# =============================================================================

# ------------------------------- SESSION A -----------------------------------
# Goal: wipe every stale artifact so the next session starts from zero.

# 1. Remove the installed package so no old DLL can be loaded by mistake.
remove.packages("GSPip")

unlink(file.path(.libPaths()[1], "00LOCK-GSPip"), recursive = TRUE)

# 2. Delete the compiled DLL from the source tree (src/*.so / src/*.dll).
#    Do this BEFORE the restart so there is nothing to unload conflicts with.
#devtools::clean_dll()

# 3. Restart R.  All loaded DLLs are released, file locks are cleared.
.rs.restartR()
# -- after restart, continue in SESSION B --------------------------------------



# ------------------------------- SESSION B -----------------------------------
# Goal: rebuild everything from source in the correct dependency order.

# 1. Regenerate example data (.rda files in data/ and flat files in
#    inst/extdata/).  Must come first because devtools::document() will
#    try to lazy-load data/ when it parses roxygen @examples.
source("data-raw/generate_example_data.R")

# 2. First document pass: generates NAMESPACE and man/ from roxygen tags in
#    R/*.R.  RcppExports.R does not exist yet so its tags are skipped.
devtools::document()


# 5. Compile C++ and install into the library.
#    upgrade = FALSE  - do not touch other packages.
devtools::install()

roxygen2::roxygenise()

devtools::load_all()

# 6. Run the test suite.  All C++ symbols are now registered in the
#    installed DLL, so load_all() will find them.
devtools::test()


# 7. Full CRAN check (run after tests pass).
devtools::check()


# 9. Build vignettes
#options(pkgdown.internet = FALSE)
library(GSPip)

# Build everything except home, then build home separately
#pkgdown::build_favicons(overwrite = TRUE)

#pkgdown::check_pkgdown()


pkgdown::build_reference()
pkgdown::build_articles()
pkgdown::build_news()

# Build home with network disabled at the curl level
httr2_mock <- function(...) stop("no network", call. = FALSE)
pkgdown::build_home()

#unloadNamespace("OptSLDP")


#pkgdown::clean_site(force = TRUE)
pkgdown::build_site()

# 10. Build package
devtools::build()

pkgdown::clean_site()
pkgdown::build_site()



setwd("D:/OneDrive - CGIAR/GitHub libraries/GSPip")

dir.create("vignettes", showWarnings = FALSE)

file.copy(
  from = "docs/articles/GSPip-workflow.md",
  to   = "vignettes/GSPip-workflow.Rmd",
  overwrite = TRUE
)

file.exists("vignettes/GSPip-workflow.Rmd")
