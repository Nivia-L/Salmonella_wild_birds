
load("C:/Users/Nivia/Desktop/Samonella WildBirds/Salmonella_wild_birds/Data-RSalmonella.RData")
# ============================================================
# FIGURA 3b – Dendrograma rectangular + heatmap de metadatos
# ggtree + gheatmap – versión publicable
# ============================================================

library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(dplyr)
install.packages("phangorn")
library(phangorn)

# ============================================================
# PASO 1 – CARGAR DATOS
# ============================================================

tree <- read.tree("snp_tree_2023-09-04.nwk")

metadata <- read.csv("Supplementary_Table_S1_metadata.csv",
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
# PASO 2b – ORDENAR POR CLADOS (ladderize)
# Reorganiza las ramas para que los clados queden agrupados
# visualmente de forma más limpia
# ============================================================

# Ladderize: ramas más cortas arriba, más largas abajo
tree <- ladderize(tree, right = FALSE)
#============================================================
  # PASO 2c – DIAGNÓSTICO DE NODOS  ← aquí
  # Correr UNA sola vez para ver qué nodos rotar
  # Luego se elimina del script final
  # ============================================================

library(phangorn)

p_temp <- ggtree(tree) %<+% meta_aligned

cat("=== NODOS INTERNOS Y SUS SOURCES ===\n")
for (nd in (n_tips + 1):(n_tips + tree$Nnode)) {
  tips <- p_temp$data[
    p_temp$data$node %in%
      Descendants(tree, nd, type = "tips")[[1]] &
      p_temp$data$isTip, ]
  srcs <- sort(unique(tips$Source))
  cat(sprintf("Nodo %3d | n=%2d | %s\n",
              nd, nrow(tips), paste(srcs, collapse = ", ")))
}

# ============================================================
# PASO 2d – ROTAR NODOS CORRECTAMENTE
# rotate() de ggtree trabaja sobre el objeto p, no sobre tree
# ============================================================

# Construir p_temp primero con el árbol ladderizado
p_temp <- ggtree(tree) %<+% meta_aligned

# Rotar nodos sobre p_temp
p_temp <- rotate(p_temp, 76)
p_temp <- rotate(p_temp, 81)
p_temp <- rotate(p_temp, 83)
p_temp <- rotate(p_temp, 85)
p_temp <- rotate(p_temp, 90)
p_temp <- rotate(p_temp, 116)
p_temp <- rotate(p_temp, 119)
p_temp <- rotate(p_temp, 127)
p_temp <- rotate(p_temp, 136)

# Verificar orden de puntas
tip_order <- p_temp$data[p_temp$data$isTip, ]
tip_order <- tip_order[order(tip_order$y, decreasing = TRUE), ]
cat("=== ORDEN DE PUNTAS DESPUÉS DE ROTACIÓN ===\n")
print(tip_order[, c("label", "Source")])

# Extraer el árbol rotado para usarlo en el script principal
tree <- as.phylo(p_temp)
# ============================================================
# PASO 2c – ROTAR NODOS PARA AGRUPAR POR SOURCE
# Objetivo: Heavy_breeders abajo, Wild Bird agrupado,
#           Layers juntos, Dogs+Pigeons juntos
# ============================================================

library(ggtree)

# Estrategia de rotación:
# Nodo 76 (raíz): rotar para poner Heavy_breeders (nodo 77) abajo
# Nodo 81: rotar para separar clados mixtos
# Nodo 83: rotar para subir Wild Bird
# Nodo 116: rotar para agrupar Dogs+Pigeons con Layers
# Nodo 119: rotar para poner Dogs+Pigeons juntos
# Nodo 127: rotar para agrupar Layers
# Nodo 136: rotar para subir Wild Bird

tree_rot <- tree  # trabajar sobre copia

# Rotar nodos clave
tree_rot <- rotate(tree_rot, 76)   # Heavy_breeders al fondo
tree_rot <- rotate(tree_rot, 81)   # separa bloque mixto
tree_rot <- rotate(tree_rot, 83)   # sube Wild Bird + Layers
tree_rot <- rotate(tree_rot, 85)   # agrupa Wild Bird
tree_rot <- rotate(tree_rot, 90)   # agrupa Wild Bird
tree_rot <- rotate(tree_rot, 116)  # Dogs+Pigeons juntos
tree_rot <- rotate(tree_rot, 119)  # Dogs+Pigeons arriba
tree_rot <- rotate(tree_rot, 127)  # Layers agrupados
tree_rot <- rotate(tree_rot, 136)  # Wild Bird al final del clado

# Verificar resultado
p_check <- ggtree(tree_rot) %<+% meta_aligned
cat("=== ORDEN DE PUNTAS DESPUÉS DE ROTACIÓN ===\n")
tip_order <- p_check$data[p_check$data$isTip, ]
tip_order <- tip_order[order(tip_order$y), ]
print(tip_order[, c("label", "Source")], n = 75)
################################
# Encontrar el nodo ancestral común de los 4 Wild Bird
library(phangorn)

wb_ids <- c("CN-Y-11", "FC-Y-11", "U154s", "U168s")
wb_nodes <- which(tree$tip.label %in% wb_ids)
mrca_wb  <- getMRCA(tree, wb_nodes)
cat("Nodo MRCA de Wild Bird:", mrca_wb, "\n")

# Ver nodos en el camino entre U168s y el resto de Wild Bird
u168_node  <- which(tree$tip.label == "U168s")
mrca_check <- getMRCA(tree, c(which(tree$tip.label == "U168s"),
                              which(tree$tip.label == "CN-Y-11")))
cat("Nodo MRCA U168s + CN-Y-11:", mrca_check, "\n")

# Ver ancestros de U168s
cat("Ancestros de U168s:\n")
print(Ancestors(tree, u168_node, type = "all"))
###############################
# ============================================================
# ROTACIONES ADICIONALES para juntar los 4 Wild Bird
# U168s está en rama 83→84→85→86→90→91→92→93
# CN-Y-11, FC-Y-11, U154s están en rama 85→90→91→92→93
# Rotar nodo 83 para subir el subclado Wild Bird
# ============================================================

p_temp <- rotate(p_temp, 83)   # sube Wild Bird hacia CN-Y-11
p_temp <- rotate(p_temp, 84)   # acerca U168s
p_temp <- rotate(p_temp, 86)   # agrupa Wild Bird + Broilers

# Verificar
tip_order2 <- p_temp$data[p_temp$data$isTip, ]
tip_order2 <- tip_order2[order(tip_order2$y, decreasing = TRUE), ]

cat("=== Wild Bird en nuevo orden ===\n")
wb_pos <- which(tip_order2$Source == "Wild Bird")
print(tip_order2[max(1, min(wb_pos)-2):min(75, max(wb_pos)+2),
                 c("label","Source")])

# Extraer árbol rotado
tree <- as.phylo(p_temp)
#####################################################
p_temp <- rotate(p_temp, 82)

# Verificar Wild Bird
tip_order3 <- p_temp$data[p_temp$data$isTip, ]
tip_order3 <- tip_order3[order(tip_order3$y, decreasing = TRUE), ]

cat("=== Posiciones Wild Bird ===\n")
print(tip_order3[tip_order3$Source == "Wild Bird", c("label","Source","y")])

# Ver contexto completo
wb_pos <- which(tip_order3$Source == "Wild Bird")
cat("\n=== Contexto alrededor de Wild Bird ===\n")
print(tip_order3[max(1, min(wb_pos)-3):min(75, max(wb_pos)+3),
                 c("label","Source")], n = 30)

tree <- as.phylo(p_temp)


#########################
# Rotar toda la cadena de ancestros de U168s
p_temp <- rotate(p_temp, 93)
p_temp <- rotate(p_temp, 92)
p_temp <- rotate(p_temp, 91)

# Verificar
tip_order4 <- p_temp$data[p_temp$data$isTip, ]
tip_order4 <- tip_order4[order(tip_order4$y, decreasing = TRUE), ]

cat("=== Posiciones Wild Bird ===\n")
print(tip_order4[tip_order4$Source == "Wild Bird",
                 c("label","Source","y")])

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

ggsave("Salmonella_fig3_heatmap.pdf", plot = p,
       width = 22, height = 28, units = "cm", device = cairo_pdf)

ggsave("Salmonella_fig3_heatmap.png", plot = p,
       width = 22, height = 28, units = "cm", dpi = 600, bg = "white")
ggsave("Salmonella_fig3_heatmap.jpg", plot = p,
       width = 22, height = 28, units = "cm", dpi = 600, bg = "white")
################
cat("Filas con Wild bird en mat_study:\n")
print(sum(mat_study$Study == "Wild bird isolates (this study)", na.rm = TRUE))

cat("\nIDs con Wild bird en metadata:\n")
print(metadata[metadata$Study == "Wild bird isolates (this study)",
               c("Isolate_ID", "Study")])

