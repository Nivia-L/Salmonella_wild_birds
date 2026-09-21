# ============================================================
# RUN ALL ANALYSES
# This wrapper only controls execution order. It does not alter
# the analytical code in the individual scripts.
# ============================================================

source("R/00_project_setup.R")

scripts <- c(
  "Scripts/PCA2_raw_rw.R",
  "Scripts/Genome_figure.R",
  "Scripts/Phylogenetic_tree_Salmonella_wild_birds_rw.R",
  "Scripts/Dendr_rect_nodos_or.R",
  "Scripts/Dendrog_and_heatmap.R",
  "Scripts/Table_and_SNPdistance_matrix.R"
)

for (s in scripts) {
  message("\n===== Running: ", s, " =====")
  source(project_file(s), echo = TRUE)
}

message("\nAll analyses completed. Check Outputs/figures and Outputs/tables.")
