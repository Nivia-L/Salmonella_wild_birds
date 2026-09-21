# Project paths only; no analytical changes.
setup_candidates <- c(file.path("R", "00_project_setup.R"), file.path("..", "R", "00_project_setup.R"))
setup_file <- setup_candidates[file.exists(setup_candidates)][1]
if (is.na(setup_file)) stop("Cannot find R/00_project_setup.R. Open Salmonella_wild_birds.Rproj before running.")
source(setup_file)

# ============================================================
# FIGURA 3 – Dendrograma rectangular publicable
# Salmonella – versión final
# ============================================================

library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(phangorn)

# ============================================================
# PASO 1 – CARGAR DATOS
# ============================================================

tree <- read.tree(data_file("snp_tree_2023-09-04.nwk"))

metadata <- read.csv(data_file("Supplementary_Table_S1_metadata.csv"),
                     stringsAsFactors = FALSE)
metadata$Isolate_ID[metadata$Isolate_ID == "U2106s (b)"] <- "b"

# ============================================================
# PASO 2 – ENRIQUECER CON POBLACIONES DEL .snp
# ============================================================

snp_ids <- c(
  "U154s","U2616s","6CT.9","U2617s","1CT0.163","5CT0.006",
  "U793s","1CT.189","U2072s","1CT0.149","1CT0.158","2CT0.109",
  "U2618s","G12A","G15A","2CT0.117","U5290s","U744s","1CT0.157",
  "3CT0.019","6CT.10","U2449s","4CT0.015","U2619s","1CT0.182",
  "3CT0.020","CN-Y-11","6CT0.006","6CTO.7","2CT.121","2CT.123",
  "G13A","1CT0.162","2CT0.081","1CT0.180","U1860s","U845s",
  "U824s","1CT.193","U2620s","1CT.192","1CT.214","U2406s",
  "U2137s","U2106s (b)","U792s","1CT.210","U5291s","U026s",
  "U2491s","3CT0.029","U2615s","U732s","U5289s","U825s",
  "2CT0.086","U2614s","U168s","1CT0.118","U842s","2CT0.084",
  "G3A","U019s","U796s","U5288s","6CT0.002","1CT0.107",
  "1CT0.128","FC-Y-11","1CTO.183","6CT.11","1CT.196","U822s",
  "6CT.12","2CT0.115"
)

snp_pops <- c(
  "Wild_Duck","Pigeons","Broilers","Pigeons","Broilers","Broilers",
  "Heavy_breeders","Broilers","Pigs","Broilers","Broilers","Broilers",
  "Pigeons","Broilers","Broilers","Broilers","Broilers","Heavy_breeders",
  "Broilers","Broilers","Broilers","Human","Broilers","Pigeons","Broilers",
  "Broilers","Wild_Duck","Broilers","Broilers","Broilers","Broilers",
  "Broilers","Broilers","Broilers","Broilers","Pigs","Layers",
  "Layers","Broilers","Pigeons","Broilers","Broilers","Human",
  "Pigs","Pigs","Heavy_breeders","Broilers","Broilers",
  "Supermarket_environment","Human","Broilers","Dogs","Heavy_breeders",
  "Broilers","Layers","Broilers","Dogs","Wild_Duck","Broilers",
  "Layers","Broilers","Broilers","Supermarket_environment","Heavy_breeders",
  "Broilers","Broilers","Broilers","Broilers","Wild_Duck","Broilers",
  "Broilers","Broilers","Layers","Broilers","Broilers"
)

snp_df <- data.frame(
  Isolate_ID = snp_ids,
  Source_snp = snp_pops,
  stringsAsFactors = FALSE
)
snp_df$Isolate_ID[snp_df$Isolate_ID == "U2106s (b)"] <- "b"

metadata <- left_join(metadata, snp_df, by = "Isolate_ID")
metadata$Source <- ifelse(!is.na(metadata$Source_snp),
                          metadata$Source_snp,
                          metadata$Source)
metadata$Source_snp <- NULL

# ============================================================
# PASO 3 – CORREGIR ISOLADOS DE ESTE ESTUDIO
# ============================================================

this_study_ids <- c("CN-Y-11", "FC-Y-11", "U154s", "U168s")

metadata$Study  <- "Other isolates"
metadata$Study[metadata$Isolate_ID %in% this_study_ids] <-
  "Wild bird isolates (this study)"
metadata$Source[metadata$Isolate_ID %in% this_study_ids] <- "Wild Bird"

cat("Conteo por Study:\n")
print(table(metadata$Study))

# ============================================================
# PASO 4 – PALETAS
# ============================================================

source_palette <- c(
  "Broilers"                = "#1f78b4",
  "Heavy_breeders"          = "#ff69b4",
  "Layers"                  = "#FFD700",
  "Pigeons"                 = "#66c2d4",
  "Pigs"                    = "#ff7f00",
  "Human"                   = "#e31a1c",
  "Wild_Duck"               = "#8B4513",
  "Dogs"                    = "#808000",
  "Supermarket_environment" = "#fdbf6f",
  "Wild Bird"               = "green"
)

# ============================================================
# PASO 5 – LADDERIZE + ROTAR NODOS
# ============================================================

# ============================================================
# PASO 5 – LADDERIZE + ROTAR NODOS
# ⚠️ Usar ape::rotate() sobre phylo directamente
# ============================================================

# ============================================================
# RECARGAR árbol original y rotar con ape::rotate()
# ============================================================

# Recargar desde el archivo original — borra el tree dañado
tree <- read.tree(data_file("snp_tree_2023-09-04.nwk"))
tree <- ladderize(tree, right = FALSE)

# Rotar con ape::rotate() — preserva longitudes de rama
tree <- ape::rotate(tree, 76)
tree <- ape::rotate(tree, 81)
tree <- ape::rotate(tree, 83)
tree <- ape::rotate(tree, 84)
tree <- ape::rotate(tree, 85)
tree <- ape::rotate(tree, 86)
tree <- ape::rotate(tree, 90)
tree <- ape::rotate(tree, 91)
tree <- ape::rotate(tree, 92)
tree <- ape::rotate(tree, 93)
tree <- ape::rotate(tree, 116)
tree <- ape::rotate(tree, 119)
tree <- ape::rotate(tree, 127)
tree <- ape::rotate(tree, 136)

# Verificar longitudes correctas
cat("Rango de longitudes:", range(tree$edge.length, na.rm = TRUE), "\n")

# Realinear metadata
meta_aligned <- metadata[match(tree$tip.label, metadata$Isolate_ID), ]
rownames(meta_aligned) <- meta_aligned$Isolate_ID

# Continuar desde PASO 6 hacia abajo

# ============================================================
# PASO 6 – ÁRBOL BASE
# ============================================================

p <- ggtree(tree, layout = "rectangular", linewidth = 0.6) %<+%
  meta_aligned +
  aes(color = Source) +
  scale_color_manual(
    name     = "Source",
    values   = c(source_palette, "Mixed" = "grey75"),
    na.value = "grey85",
    breaks   = names(source_palette),   # ← solo muestra los de la paleta, excluye NA
    guide    = guide_legend(
      order        = 1,
      override.aes = list(linewidth = 2.5, shape = NA)
    )
  ) +
  # Triángulos para aislamientos de este estudio
  geom_tippoint(
    data    = function(d) subset(d, label %in% this_study_ids),
    mapping = aes(x = x, y = y),
    shape   = 17,
    size    = 2.5,
    color   = "green"
  ) +
  # Tip labels
  geom_tiplab(
    aes(label = label),
    size     = 3.0,
    offset   = 0.001,
    hjust    = 0,
    color    = "black"
  )

# ============================================================
# PASO 7 – LEYENDA STUDY (triángulo fantasma)
# ============================================================

p <- p +
  geom_point(
    data        = data.frame(x = -Inf, y = -Inf,
                             Study = "Wild bird isolates (this study)"),
    aes(x = x, y = y, shape = Study),
    color       = "green",
    size        = 0.001,
    inherit.aes = FALSE
  ) +
  scale_shape_manual(
    name   = "Study",
    values = c("Wild bird isolates (this study)" = 17),
    guide  = guide_legend(
      order        = 2,
      override.aes = list(color = "green", size = 3)
    )
  )

# ============================================================
# PASO 8 – TEMA PUBLICABLE
# ============================================================

p <- p +
  theme_tree2() +
  labs(x = "SNP Distances") +
  theme(
    legend.position    = "right",
    legend.title       = element_text(size = 12, face = "bold"),
    legend.text        = element_text(size = 10),
    legend.key.size    = unit(0.5, "cm"),
    legend.spacing.y   = unit(0.2, "cm"),
    legend.box.spacing = unit(0.3, "cm"),
    axis.text.x        = element_text(size = 10, color = "black"),
    axis.title.x       = element_text(size = 12, face = "bold",
                                      margin = margin(t = 6)),
    plot.margin        = margin(15, 10, 10, 10)
  )

p

# ============================================================
# PASO 9 – EXPORTAR
# ============================================================

ggsave(figure_file("Salmonella_fig3rv.pdf"), plot = p,
       width = 22, height = 35, units = "cm",
       device = cairo_pdf, limitsize = FALSE)

ggsave(figure_file("Salmonella_fig3rv.png"), plot = p,
       width = 22, height = 35, units = "cm",
       dpi = 600, bg = "white", limitsize = FALSE)
ggsave(figure_file("Salmonella_fig3rv.jpg"), plot = p,
       width = 22, height = 35, units = "cm",
       dpi = 600, bg = "white", limitsize = FALSE)
