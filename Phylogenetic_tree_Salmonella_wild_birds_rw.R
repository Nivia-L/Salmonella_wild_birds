# ============================================================
# ÁRBOL FILOGENÉTICO CIRCULAR
# Estilo: ramas por Source + anillo externo por Study
# ============================================================

# PASO 0 – LIBRERÍAS
# Instalar solo la primera vez:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("ape","treeio","ggtree","ggtreeExtra"))
# install.packages(c("ggplot2","ggnewscale","dplyr"))

library(ape)
library(treeio)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(dplyr)

# ============================================================
# PASO 1 – CARGAR DATOS
# ============================================================

tree <- read.tree("snp_tree_2023-09-04.nwk")  # ← cambia este nombre

metadata <- read.csv("Supplementary_Table_S1_metadata.csv",
                     stringsAsFactors = FALSE)

metadata$Isolate_ID[metadata$Isolate_ID == "U2106s (b)"] <- "b"

# ============================================================
# PASO 1b – ENRIQUECER CON POBLACIONES DEL .snp
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
# PASO 1c – AISLAMIENTOS DE ESTE ESTUDIO
# ============================================================

this_study_ids <- c("CN-Y-11", "FC-Y-11", "U154s", "U168s")
metadata$Source[metadata$Isolate_ID %in% this_study_ids] <- "Wild Bird"
metadata$Study[metadata$Isolate_ID %in% this_study_ids]  <- "Wild bird isolates (this study)"

cat("Sources únicos:\n")
print(sort(unique(metadata$Source)))

# ============================================================
# PASO 2 – PALETAS
# ============================================================

# Colores de ramas por Source
source_palette <- c(
  "Broilers"                = "#1f78b4",
  "Heavy_breeders"          = "#ff69b4",
  "Layers"                  = "yellow",
  "Pigeons"                 = "#66c2d4",
  "Pigs"                    = "#ff7f00",
  "Human"                   = "#e31a1c",
  "Wild_Duck"               = "#8B4513",
  "Dogs"                    = "#808000",
  "Supermarket_environment" = "#fdbf6f",
  "Wild Bird"               = "green",
  "Mixed"                   = "grey75"
)

# Color del anillo externo por Study
study_palette <- c(
  "Wild bird isolates (this study)" = "green")
 

# ============================================================
# PASO 3 – VERIFICAR COINCIDENCIA
# ============================================================

cat(sprintf("Puntas en metadata: %d / %d\n",
            sum(tree$tip.label %in% metadata$Isolate_ID),
            length(tree$tip.label)))

missing <- setdiff(tree$tip.label, metadata$Isolate_ID)
if (length(missing) > 0) { cat("Sin metadata:\n"); print(missing) }

# ============================================================
# PASO 4 – PROPAGAR SOURCE A NODOS INTERNOS
# ============================================================

n_tips  <- length(tree$tip.label)
n_nodes <- tree$Nnode
edge    <- tree$edge

tip_source <- metadata$Source[match(tree$tip.label, metadata$Isolate_ID)]
tip_study  <- metadata$Study[ match(tree$tip.label, metadata$Isolate_ID)]
tip_study[is.na(tip_study)] <- "Other isolates"

node_source <- c(tip_source, rep(NA_character_, n_nodes))

for (i in nrow(edge):1) {
  p_idx <- edge[i, 1]
  c_idx <- edge[i, 2]
  cs    <- node_source[c_idx]
  ps    <- node_source[p_idx]
  if (is.na(ps)) {
    node_source[p_idx] <- cs
  } else if (!is.na(cs) && ps != cs) {
    node_source[p_idx] <- "Mixed"
  }
}

node_df <- data.frame(
  node     = 1:(n_tips + n_nodes),
  Source   = node_source,
  Study    = c(tip_study, rep(NA_character_, n_nodes)),
  is_study = c(tree$tip.label %in% this_study_ids, rep(FALSE, n_nodes)),
  stringsAsFactors = FALSE
)

cat(sprintf("Puntas sin Source: %d\n", sum(is.na(tip_source))))

# ============================================================
# PASO 5 – ÁRBOL BASE
# Ramas coloreadas por Source usando I() para evitar conflictos
# ============================================================

p <- ggtree(tree, layout = "circular", linewidth = 0.6) %<+% node_df +
  aes(color = Source) +
  scale_color_manual(
    name     = "Source",
    values   = c(source_palette, "Mixed" = "grey75"),
    na.value = "grey80",
    breaks   = names(source_palette)[names(source_palette) != "Mixed"],
    guide    = guide_legend(
      title.position = "top",
      ncol           = 1,
      override.aes   = list(linewidth = 2.5, shape = NA)
    )
  )
# ============================================================
# PASO 6 – LABELS DE PUNTAS
# Van primero, el anillo externo irá encima
# ============================================================

# ============================================================
# PASO 6 – LABELS (van ANTES del anillo externo)
# ============================================================

p <- p +
  geom_tiplab(
    data     = subset(p$data, isTip & !is_study),
    aes(label = label),
    color    = "grey20",
    size     = 2.8,
    offset   = 0.05,
    hjust    = 0
  ) +
  geom_tiplab(
    data     = subset(p$data, isTip & is_study),
    aes(label = label),
    color    = "green",
    fontface = "bold",
    size     = 3.0,
    offset   = 0.1,
    hjust    = 0
  )

# ============================================================
# PASO 7 – ANILLO EXTERNO (va DESPUÉS de los labels)
# ── offset alto para que quede fuera de los labels
# ============================================================

ring_df <- data.frame(
  Isolate_ID  = tree$tip.label,
  Source_ring = tip_source,
  stringsAsFactors = FALSE
)
ring_df$Source_ring[is.na(ring_df$Source_ring)] <- "Unknown"

p <- p +
  new_scale_fill() +
  suppressWarnings(
    geom_fruit(
      data    = ring_df,
      geom    = geom_tile,
      mapping = aes(y = Isolate_ID, fill = Source_ring),
      offset  = 0.5,    # ← alto para quedar fuera de los labels
      pwidth  = 0.5,   # ← grosor del anillo
      color   = NA
    )
  ) +
  scale_fill_manual(
    name     = "Source ",   # espacio para no duplicar título en leyenda
    values   = c(source_palette, "Unknown" = "grey90"),
    guide    = "none"       # sin leyenda extra, la de ramas ya la cubre
  )

# ============================================================
# PASO 8 – LEYENDA STUDY
# ============================================================

p <- p +
  scale_shape_manual(
    name   = "Study",
    values = c("Wild bird isolates (this study)" = 15),
    guide  = guide_legend(
      title.position = "top",
      ncol           = 1,
      override.aes   = list(color = "green", size = 4)
    )
  ) +
  geom_point(
    data        = data.frame(x = Inf, y = Inf,
                             Study = "Wild bird isolates (this study)"),
    aes(x = x, y = y, shape = Study),
    color       = "white",   # ← invisible en el plot
    size        = 0.001,
    inherit.aes = FALSE
  )

# ============================================================
# PASO 9 – TEMA
# ============================================================
print(p)
# ============================================================
# PASO 9 – GUARDAR
# ============================================================

ggsave("phylo_circular_tree.png",
       plot = p, width = 16, height = 16,
       units = "in", dpi = 300, bg = "white")

ggsave("phylo_circular_tree.jpg",
       plot = p, width = 16, height = 16,
       units = "in", dpi = 300, bg = "white")
