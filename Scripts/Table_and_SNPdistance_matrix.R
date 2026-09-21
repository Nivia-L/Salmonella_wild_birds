# Project paths only; no analytical changes.
setup_candidates <- c(file.path("R", "00_project_setup.R"), file.path("..", "R", "00_project_setup.R"))
setup_file <- setup_candidates[file.exists(setup_candidates)][1]
if (is.na(setup_file)) stop("Cannot find R/00_project_setup.R. Open Salmonella_wild_birds.Rproj before running.")
source(setup_file)

library(ape)
library(dplyr)
library(ggplot2)
library(reshape2)

tree <- read.tree(data_file("snp_tree_2023-09-04.nwk"))
metadata <- read.csv(data_file("Supplementary_Table_S1_metadata.csv"), stringsAsFactors = FALSE)
metadata$Isolate_ID[metadata$Isolate_ID == "U2106s (b)"] <- "b"

# ============================================================
# TABLA DE DISTANCIA MEDIA ENTRE GRUPOS (Source)
# CORREGIDO: match por rownames directamente
# ============================================================

# Matriz de distancias
dist_matrix <- cophenetic.phylo(tree)

# ⚠️ FIX: alinear sources AL ORDEN de rownames(dist_matrix)
sources <- metadata$Source[match(rownames(dist_matrix), metadata$Isolate_ID)]

# Verificar que el match es correcto
cat("Verificación primeras 6 filas:\n")
print(data.frame(
  tip    = rownames(dist_matrix)[1:6],
  source = sources[1:6]
))

# Grupos únicos sin NA
grupos <- sort(unique(sources[!is.na(sources)]))

# Matriz resultado
result <- matrix(NA, nrow = length(grupos), ncol = length(grupos),
                 dimnames = list(grupos, grupos))

for (i in seq_along(grupos)) {
  for (j in seq_along(grupos)) {
    ids_i <- rownames(dist_matrix)[!is.na(sources) & sources == grupos[i]]
    ids_j <- rownames(dist_matrix)[!is.na(sources) & sources == grupos[j]]
    
    if (length(ids_i) == 0 || length(ids_j) == 0) next
    
    sub <- dist_matrix[ids_i, ids_j, drop = FALSE]
    
    if (i == j) {
      # Intra-grupo: solo triángulo superior sin diagonal
      vals <- sub[upper.tri(sub)]
    } else {
      vals <- as.vector(sub)
    }
    
    if (length(vals) > 0 && !all(is.na(vals))) {
      result[i, j] <- round(mean(vals, na.rm = TRUE), 6)
    }
  }
}

# Mostrar
cat("\n=== DISTANCIA MEDIA ENTRE GRUPOS (unidades de rama) ===\n")
print(as.data.frame(result))

# Exportar
write.csv(as.data.frame(result), table_file("distancia_media_grupos.csv"), row.names = TRUE)

# Versión larga
dist_long <- data.frame(
  Group_1           = character(),
  Group_2           = character(),
  Mean_distance     = numeric(),
  n_comparisons     = integer(),
  stringsAsFactors  = FALSE
)

for (i in seq_along(grupos)) {
  for (j in i:length(grupos)) {
    ids_i <- rownames(dist_matrix)[!is.na(sources) & sources == grupos[i]]
    ids_j <- rownames(dist_matrix)[!is.na(sources) & sources == grupos[j]]
    sub   <- dist_matrix[ids_i, ids_j, drop = FALSE]
    vals  <- if (i == j) sub[upper.tri(sub)] else as.vector(sub)
    vals  <- vals[!is.na(vals)]
    
    dist_long <- rbind(dist_long, data.frame(
      Group_1       = grupos[i],
      Group_2       = grupos[j],
      Mean_distance = round(mean(vals), 6),
      n_comparisons = length(vals)
    ))
  }
}

cat("\n=== VERSIÓN LARGA (ordenada por distancia) ===\n")
print(dist_long[order(dist_long$Mean_distance), ], row.names = FALSE)

write.csv(dist_long, table_file("distancia_media_grupos_long.csv"), row.names = FALSE)
######################
# ============================================================
# CONVERTIR A SNPs Y GENERAR TABLA PUBLICABLE
# ============================================================

# Salmonella enterica: genoma ~4,800,000 bp
genome_size <- 4800000

# Convertir distancias de rama a SNPs
dist_long$Mean_SNPs <- round(dist_long$Mean_distance * genome_size, 1)

# ============================================================
# TABLA RESUMEN PUBLICABLE – solo pares únicos (sin diagonal duplicada)
# ============================================================

tabla_pub <- dist_long %>%
  mutate(
    Mean_SNPs     = round(Mean_distance * genome_size, 1),
    Comparison    = paste(Group_1, "vs", Group_2)
  ) %>%
  select(Group_1, Group_2, Mean_SNPs, n_comparisons) %>%
  rename(
    "Source 1"       = Group_1,
    "Source 2"       = Group_2,
    "Mean SNP dist." = Mean_SNPs,
    "n pairs"        = n_comparisons
  )

cat("=== TABLA PUBLICABLE (SNP distances) ===\n")
print(tabla_pub, row.names = FALSE)

# Exportar tabla publicable
write.csv(tabla_pub, table_file("Table_SNP_distances.csv"), row.names = FALSE)

# ============================================================
# MATRIZ CUADRADA EN SNPs para suplementario
# ============================================================

result_snp <- round(result * genome_size, 1)

cat("\n=== MATRIZ EN SNPs ===\n")
print(as.data.frame(result_snp))

write.csv(as.data.frame(result_snp),
          table_file("Table_SNP_distances_matrix.csv"))

cat("\nArchivos exportados:\n")
cat("  - Table_SNP_distances.csv      (tabla larga)\n")
cat("  - Table_SNP_distances_matrix.csv (matriz)\n")
#############
# Exportar como tabla formateada
library(ggplot2)
library(reshape2)

result_snp_df <- as.data.frame(result_snp)
result_snp_df$Group1 <- rownames(result_snp_df)

dist_melt <- melt(result_snp_df, id.vars = "Group1",
                  variable.name = "Group2",
                  value.name = "SNP_dist")

ggplot(dist_melt, aes(x = Group2, y = Group1, fill = SNP_dist)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(SNP_dist == 0, "—",
                               as.character(SNP_dist))),
            size = 3, color = "white") +
  scale_fill_viridis_c(
    name     = "Mean SNP\ndistance",
    option   = "viridis",   # opciones: "viridis", "magma", "plasma", "cividis"
    na.value = "grey90"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 9)
  )

ggsave(figure_file("Fig_SNP_distances_heatmap.pdf"),
       width = 18, height = 14, units = "cm", device = cairo_pdf)

ggsave(figure_file("Fig_SNP_distances_heatmap.jpg"),
       width  = 18,
       height = 14,
       units  = "cm",
       dpi    = 600,
       bg     = "white")
