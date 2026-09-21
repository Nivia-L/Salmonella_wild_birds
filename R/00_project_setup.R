# ============================================================
# PROJECT SETUP – Salmonella_wild_birds
# This file only defines project paths. It does not alter analyses.
# ============================================================

find_project_root <- function() {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidates <- unique(c(wd, dirname(wd)))
  for (p in candidates) {
    if (file.exists(file.path(p, "Salmonella_wild_birds.Rproj"))) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Project root not found. Open Salmonella_wild_birds.Rproj before running the scripts.")
}

project_root <- find_project_root()
data_dir <- file.path(project_root, "Data")
output_dir <- file.path(project_root, "Outputs")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")

for (d in c(output_dir, figure_dir, table_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Helper for reproducible paths
project_file <- function(...) file.path(project_root, ...)
data_file <- function(...) file.path(data_dir, ...)
figure_file <- function(...) file.path(figure_dir, ...)
table_file <- function(...) file.path(table_dir, ...)
