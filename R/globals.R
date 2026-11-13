# data.table expressions are intentionally evaluated within data-table columns.
# Declaring these symbols prevents R CMD check from treating that non-standard
# evaluation as missing package-level objects.
utils::globalVariables(c(".", ":=", "Fold", "Iter", "Level", "Metric", "Value"))
