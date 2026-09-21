# Project paths only; no analytical changes.
setup_candidates <- c(file.path("R", "00_project_setup.R"), file.path("..", "R", "00_project_setup.R"))
setup_file <- setup_candidates[file.exists(setup_candidates)][1]
if (is.na(setup_file)) stop("Cannot find R/00_project_setup.R. Open Salmonella_wild_birds.Rproj before running.")
source(setup_file)

## ============================================================
# FIGURA 2 – PCA corregido con aislamientos de este estudio
# ============================================================

library(adegenet)
library(ggplot2)

# ------------------------------------------------------------
# 1. CARGAR DATOS
# ------------------------------------------------------------

Data_PCA2 <- read.snp(data_file("Data_PCA2.snp"))

# ------------------------------------------------------------
# 2. PCA
# ------------------------------------------------------------

pca_res <- glPca(Data_PCA2, nf = 3)
var     <- pca_res$eig / sum(pca_res$eig) * 100

# ------------------------------------------------------------
# 3. DATA FRAME con IDs y población del archivo .snp
# ------------------------------------------------------------

pca_df         <- as.data.frame(pca_res$scores)
pca_df$ID      <- indNames(Data_PCA2)
pca_df$Source  <- as.character(pop(Data_PCA2))  # población original del .snp

# ⚠️ FC-Y-07 no existe en el archivo — se usan los 4 IDs confirmados:
#    CN-Y-11, FC-Y-11, U154s, U168s
this_study_ids <- c("CN-Y-11", "FC-Y-11", "U154s", "U168s")

# Columna Study: distingue "Wild Bird (this study)" del resto
pca_df$Study <- ifelse(
  pca_df$ID %in% this_study_ids,
  "Wild Bird (this study)",
  pca_df$Source
)

# Verificar asignación
cat("Aislamientos de este estudio encontrados:\n")
print(pca_df[pca_df$Study == "Wild Bird (this study)", c("ID","Source","Study")])

cat("\nDistribución completa:\n")
print(table(pca_df$Study))

# ------------------------------------------------------------
# 4. PALETA – todos los grupos visibles
# ------------------------------------------------------------

source_palette_pca <- c(
  "Broilers"                = "#1f78b4",
  "Heavy_breeders"          = "#ff69b4",
  "Layers"                  = "#33a02c",
  "Pigeons"                 = "#66c2d4",
  "Pigs"                    = "#ff8c00",
  "Wild_Duck"               = "#b15928",
  "Human"                   = "#e31a1c",
  "Dogs"                    = "#808000",
  "Supermarket_environment" = "#f0e6b6",
  "Wild Bird (this study)"  = "green"
)

# ------------------------------------------------------------
# 5. GRÁFICO
# ------------------------------------------------------------

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Study)) +
  
  # Líneas de referencia
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  
  # Elipses por grupo (excluyendo Wild Bird — solo 4 puntos)
  stat_ellipse(
    data      = subset(pca_df, Study != "Wild Bird (this study)"),
    aes(group = Study),
    linewidth = 0.35,
    linetype  = "dashed",
    alpha     = 0.5
  ) +
  
  # Puntos regulares
  geom_point(
    data  = subset(pca_df, Study != "Wild Bird (this study)"),
    size  = 2.5,
    alpha = 0.85
  ) +
  
  # Wild Bird (este estudio) — triángulo más grande encima
  geom_point(
    data  = subset(pca_df, Study == "Wild Bird (this study)"),
    size  = 4,
    shape = 17
  ) +
  
  scale_color_manual(values = source_palette_pca) +
  coord_equal() +
  theme_classic(base_size = 12) +
  labs(
    x     = paste0("PC1 (", round(var[1], 1), "%)"),
    y     = paste0("PC2 (", round(var[2], 1), "%)"),
    color = "Source"
  ) +
  theme(
    legend.title    = element_text(face = "bold", size = 10),
    legend.text     = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    plot.margin     = margin(10, 10, 10, 10)
  )

p_pca


# ------------------------------------------------------------
# 5. GRÁFICO (LIMPIO + LÍNEAS DE COORDENADAS)
# ------------------------------------------------------------

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Study)) +
  
  # Líneas de coordenadas (gris suave)
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  
  # Puntos regulares
  geom_point(
    data  = subset(pca_df, Study != "Wild Bird (this study)"),
    size  = 2.5,
    alpha = 0.85
  ) +
  
  # Wild Bird (este estudio) — triángulo destacado
  geom_point(
    data  = subset(pca_df, Study == "Wild Bird (this study)"),
    size  = 4,
    shape = 17
  ) +
  
  scale_color_manual(values = source_palette_pca) +
  
  coord_equal() +
  
  theme_classic() +
  
  labs(
    x     = paste0("PC1 (", round(var[1], 1), "%)"),
    y     = paste0("PC2 (", round(var[2], 1), "%)"),
    color = "Source"
  ) +
  
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 12),
    legend.title    = element_text(face = "bold", size = 12),
    legend.text     = element_text(size = 11),
    legend.key.size = unit(0.4, "cm"),
    plot.margin     = margin(10, 10, 10, 10)
  )

p_pca

# ------------------------------------------------------------
# 6. EXPORTAR
# ------------------------------------------------------------
ggsave(figure_file("Figure2_PCA.pdf"), plot = p_pca,
       width = 8, height = 8, device = cairo_pdf)

ggsave(figure_file("Figure2_PCA.png"), plot = p_pca,
       width = 8, height = 6, dpi = 600, bg = "white")
ggsave(figure_file("Figure2_PCA.jpg"), plot = p_pca,
       width = 8, height = 6, dpi = 600, bg = "white")
##############

