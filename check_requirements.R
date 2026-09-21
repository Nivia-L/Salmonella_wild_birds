# Check packages without installing anything.
required <- c(
  "ape", "adegenet", "dplyr", "ggplot2", "ggnewscale",
  "ggtree", "ggtreeExtra", "phangorn", "reshape2", "treeio"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) == 0) {
  message("All required R packages are installed.")
} else {
  message("Missing packages: ", paste(missing, collapse = ", "))
  message("Install CRAN packages with install.packages() and Bioconductor packages with BiocManager::install().")
}
