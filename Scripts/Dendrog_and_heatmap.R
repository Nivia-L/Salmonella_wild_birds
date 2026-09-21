# Project paths only; no analytical changes.
setup_candidates <- c(file.path("R", "00_project_setup.R"), file.path("..", "R", "00_project_setup.R"))
setup_file <- setup_candidates[file.exists(setup_candidates)][1]
if (is.na(setup_file)) stop("Cannot find R/00_project_setup.R. Open Salmonella_wild_birds.Rproj before running.")
source(setup_file)


load(data_file("Data-RSalmonella.RData"))
# ============================================================
# FIGURA 3b – Dendrograma rectangular + heatmap de metadatos
# ggtree + gheatmap – versión publicable
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

# Align metadata to tree tips for plotting and diagnostics.
meta_aligned <- metadata[match(tree$tip.label, metadata$Isolate_ID), , drop = FALSE]
rownames(meta_aligned) <- meta_aligned$Isolate_ID

# ============================================================
# PASO 2b – ORDENAR Y ROTAR EL ÁRBOL
# The rotation sequence below is unchanged from the analysis.
# It only controls the visual arrangement of the final tree.
# ============================================================

tree <- ladderize(tree, right = FALSE)
p_temp <- ggtree(tree) %<+% meta_aligned

p_temp <- rotate(p_temp, 76)
p_temp <- rotate(p_temp, 81)
p_temp <- rotate(p_temp, 83)
p_temp <- rotate(p_temp, 85)
p_temp <- rotate(p_temp, 90)
p_temp <- rotate(p_temp, 116)
p_temp <- rotate(p_temp, 119)
p_temp <- rotate(p_temp, 127)
p_temp <- rotate(p_temp, 136)
p_temp <- rotate(p_temp, 83)
p_temp <- rotate(p_temp, 84)
p_temp <- rotate(p_temp, 86)
p_temp <- rotate(p_temp, 82)
p_temp <- rotate(p_temp, 93)
p_temp <- rotate(p_temp, 92)
p_temp <- rotate(p_temp, 91)

tree <- as.phylo(p_temp)

# ============================================================
# PASO 3 – CORREGIR ISOLADOS DE ESTE ESTUDIO
# ============================================================

this_study_ids <- c("CN-Y-11", "FC-Y-11", "U154s", "U168s")

# ⚠️ Resetear TODOS a "Other isolates" primero
metadata$Study <- "Other isolates"

# Luego asignar solo los 4 confirmados
metadata$Study[metadata$Isolate_ID %in% this_study_ids] <-
  "Wild bird isolates (this study)"
metadata$Source[metadata$Isolate_ID %in% this_study_ids] <- "Wild Bird"

# Verificar — debe dar exactamente 4
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

study_palette <- c(
  "Wild bird isolates (this study)" = "green")


# ============================================================
# PASO 5 – PREPARAR MATRICES PARA GHEATMAP
# Gheatmap necesita un data frame con rownames = tip labels
# ============================================================

# Reconstruir mat_study con el Study corregido
meta_aligned <- metadata[match(tree$tip.label, metadata$Isolate_ID), ]
rownames(meta_aligned) <- meta_aligned$Isolate_ID

mat_study <- data.frame(
  Study = meta_aligned$Study,
  row.names = rownames(meta_aligned)
)

mat_source <- data.frame(
  Source = meta_aligned$Source,
  row.names = rownames(meta_aligned)
)

# Verificar
cat("Celdas Wild bird en mat_study:", 
    sum(mat_study$Study == "Wild bird isolates (this study)"), "\n")
# ============================================================
# PASO 6 – ÁRBOL BASE rectangular
# ============================================================

p <- ggtree(tree, layout = "rectangular", linewidth = 0.6) %<+% meta_aligned +
  aes(color = Source) +
  scale_color_manual(
    values   = c(source_palette, "Mixed" = "grey75"),
    na.value = "grey85",
    guide    = "none"    # la leyenda va en el heatmap
  ) +
  # Triángulos para this_study
  geom_tippoint(
    data    = function(d) subset(d, label %in% this_study_ids),
    mapping = aes(x = x, y = y),
    shape   = 17,
    size    = 2,
    color   = "green"
  ) +
  # Tip labels
  geom_tiplab(
    aes(label = label),
    size     = 3.0,
    offset   = 0.001,
    hjust    = 0,
    color    = "Black"
      )

# ============================================================
# PASO 7 – HEATMAP COLUMNA 1: SOURCE
# ============================================================
max_x              <- max(p$data$x, na.rm = TRUE)
offset_1           <- max_x * 0.07       # espacio para tip labels
source_col_width_abs <- max_x * 0.03      # ancho absoluto columna Source
gap                <- max_x * 0.04        # espacio entre columnas
study_offset       <- offset_1 + source_col_width_abs + gap  # offset Study

col_width   <- max_x * 0.03    # ancho de cada columna
gap         <- max_x * 0.04    # espacio entre columnas
p <- gheatmap(
  p,
  mat_source,
  offset            = offset_1,          # distancia desde tip labels
  width             = 0.05,
  colnames          = TRUE,
  colnames_position = "top",
  colnames_angle    = 0,
  colnames_offset_y = 1,
  font.size         = 3.5,
  color             = "white"
) +
  scale_fill_manual(
    name     = "Source",
    values   = source_palette,
    na.value = "grey90",
    guide    = guide_legend(order = 1,
                            override.aes = list(color = NA))
  )

# ============================================================
# PASO 8 – HEATMAP COLUMNA 2: STUDY
# ⚠️ new_scale_fill() es obligatorio entre columnas de gheatmap
# ============================================================
# Punto fantasma invisible para forzar entrada en leyenda
p <- p +
  geom_point(
    data        = data.frame(x = -Inf, y = -Inf,
                             Study = "Wild bird isolates (this study)"),
    aes(x = x, y = y, shape = Study),
    color       = "green",
    size        = 4,
    inherit.aes = FALSE
  ) +
  scale_shape_manual(
    name   = "Study",
    values = c("Wild bird isolates (this study)" = 17),  # triángulo
    guide  = guide_legend(
      order        = 2,
      override.aes = list(color = "green", size = 4)
    )
  )

# ============================================================
# PASO 9 – TEMA PUBLICABLE
# ============================================================

p <- p +
  theme_tree2() +
  labs(
    x = "SNP Distances"    # ← título del eje X
  ) +
  theme(
    legend.position    = "right",
    legend.title       = element_text(size = 12,  face = "bold"),
    legend.text        = element_text(size = 8),
    legend.key.size    = unit(0.5, "cm"),
    legend.spacing.y   = unit(0.2, "cm"),
    legend.box.spacing = unit(0.3, "cm"),
    axis.text.x        = element_text(size = 10, color = "black"),
    axis.title.x       = element_text(size = 12, face = "bold",   # ← estilo título
                                      margin = margin(t = 6)),
    plot.margin        = margin(15, 10, 10, 10)
  )
p
# ============================================================
# PASO 10 – EXPORTAR
# ============================================================

ggsave(figure_file("Salmonella_fig3_heatmap.pdf"), plot = p,
       width = 22, height = 28, units = "cm", device = cairo_pdf)

ggsave(figure_file("Salmonella_fig3_heatmap.png"), plot = p,
       width = 22, height = 28, units = "cm", dpi = 600, bg = "white")
ggsave(figure_file("Salmonella_fig3_heatmap.jpg"), plot = p,
       width = 22, height = 28, units = "cm", dpi = 600, bg = "white")
################
cat("Filas con Wild bird en mat_study:\n")
print(sum(mat_study$Study == "Wild bird isolates (this study)", na.rm = TRUE))

cat("\nIDs con Wild bird en metadata:\n")
print(metadata[metadata$Study == "Wild bird isolates (this study)",
               c("Isolate_ID", "Study")])

