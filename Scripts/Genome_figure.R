# Project paths only; no analytical changes.
setup_candidates <- c(file.path("R", "00_project_setup.R"), file.path("..", "R", "00_project_setup.R"))
setup_file <- setup_candidates[file.exists(setup_candidates)][1]
if (is.na(setup_file)) stop("Cannot find R/00_project_setup.R. Open Salmonella_wild_birds.Rproj before running.")
source(setup_file)
##############
library(ggtree)
library(ggplot2)
library(ape)
library(treeio)
############
tree <- read.tree(data_file("snp_tree_2023-09-04.nwk"))
plot(tree)
metadata <- read.csv(data_file("Supplementary_Table_S1_metadata.csv"))
metadata$label <- metadata$Isolate_ID

# 1. PREPARACIÓN RADICAL
# Eliminamos el outlier 'b' que distorsiona la escala del círculo central
tree_final <- drop.tip(tree, "b") 
tree_final <- ladderize(tree_final)

metadata_clean <- metadata
names(metadata_clean)[10] <- "Label_Original"
metadata_clean <- metadata_clean[!is.na(metadata_clean$Source), ]

# 2. CÁLCULO DE COORDENADAS (Sin solapamientos)
max_x <- max(fortify(tree_final)$x, na.rm = TRUE)

# Definimos distancias claras para cada "pista"
anillo_interno <- max_x * 1.08  # Banda de color pegada al árbol
distancia_texto <- max_x * 1.20  # Nombres de muestras
anillo_externo <- max_x * 1.70   # Segunda banda (opcional)
limite_canvas   <- max_x * 2.4   # Espacio total para que no se corte la leyenda

# 3. PALETA DE COLORES (Tu selección)
source_palette <- c(
  "Broilers" = "#1f78b4", "Heavy breeders" = "#4daf4a", "Layers" = "#33a02c",
  "Pigeons" = "#66c2d4", "Pigs" = "#ff8c00", "Wild Duck" = "#b15928",
  "Human" = "#e31a1c", "Dogs" = "#808000", "Supermarket environment" = "#f0e6b6"
)

# 4. CONSTRUCCIÓN POR CAPAS
p <- ggtree(tree_final, layout = "circular", linewidth = 0.5, color = "grey40") %<+% metadata_clean +
  
  # CAPA A: Círculos en las puntas (Transparentes para no saturar)
  geom_tippoint(aes(color = Source), size = 2.5, alpha = 0.8) +
  
  # CAPA B: LA BANDA DE COLOR (Simulando un anillo sólido)
  # Usamos una barra vertical muy ancha para que se toquen entre sí
  geom_tippoint(aes(color = Source, x = anillo_interno), 
                shape = 16, size = 4, na.rm = TRUE) +
  
  # CAPA C: ETIQUETAS DE TEXTO (Alineadas y limpias)
  # 'align = TRUE' es lo que crea las líneas guía punteadas hacia el nombre
  geom_tiplab(size = 3.2, offset = distancia_texto - max_x, 
              color = "black", fontface = "plain", 
              align = TRUE, linetype = "dotted", linesize = 0.3) +
  
  # CAPA D: CONFIGURACIÓN VISUAL
  scale_color_manual(values = source_palette, na.translate = FALSE) +
  theme_tree() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(family = "serif", size = 10),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  
  # CAPA E: AJUSTE DE LIENZO
  xlim(0, limite_canvas)

# 5. RENDERIZADO LIMPIO
print(p)
###############
# Este paso es fundamental para que no se vea pixelado al imprimir
ggsave(figure_file("Salmonella_fig3.png"), plot = p, width = 12, height = 12, 
       units = "in", dpi = 600, bg = "white")

print(p)

