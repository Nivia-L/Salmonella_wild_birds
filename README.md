# Genomic Evidence of Multidrug-Resistant *Salmonella* in Wild Waterbirds from High-Andean Lakes of Ecuador

This repository contains the R scripts and supporting files used for the genomic analyses reported in the manuscript **“Genomic Evidence of Multidrug-Resistant Salmonella in Wild Waterbirds from High-Andean Lakes of Ecuador.”**

## Repository organization

```text
Salmonella_wild_birds/
├── Salmonella_wild_birds.Rproj
├── README.md
├── .gitignore
├── R/
│   └── 00_project_setup.R
├── Data/
│   ├── Data_PCA2.snp
│   ├── Data-RSalmonella.RData
│   ├── Supplementary_Table_S1_metadata.csv
│   ├── Supplementary_Table_S2_Source_distribution.csv
│   ├── Supplementary_Table_S3_Genomes_used_in_genomic_analysis.csv
│   └── snp_tree_2023-09-04.nwk
├── Scripts/
│   ├── Dendr_rect_nodos_or.R
│   ├── Dendrog_and_heatmap.R
│   ├── Genome_figure.R
│   ├── PCA2_raw_rw.R
│   ├── Phylogenetic_tree_Salmonella_wild_birds_rw.R
│   └── Table_and_SNPdistance_matrix.R
└── Outputs/
    ├── figures/
    └── tables/
```

## Reproducibility

The scripts retain the analytical procedures used for the manuscript. The organization was changed only to make file paths portable and outputs reproducible on another computer.

The scripts no longer install packages automatically and do not contain computer-specific absolute paths. Input files are read from `Data/`, and generated figures and tables are written to `Outputs/figures/` and `Outputs/tables/`.

### Requirements

R and RStudio are recommended. The following packages are required:

- `ape`
- `adegenet`
- `dplyr`
- `ggplot2`
- `ggnewscale`
- `ggtree`
- `ggtreeExtra`
- `phangorn`
- `reshape2`
- `treeio`
- `BiocManager` (for installation of Bioconductor packages)

`ggtree`, `ggtreeExtra`, and `treeio` are Bioconductor packages.

Install packages once, outside the analysis scripts. For example:

```r
install.packages(c("ape", "adegenet", "dplyr", "ggplot2", "ggnewscale",
                  "phangorn", "reshape2", "BiocManager"))
BiocManager::install(c("ggtree", "ggtreeExtra", "treeio"))
```

### Running the analysis

1. Download or clone the repository.
2. Open `Salmonella_wild_birds.Rproj` in RStudio.
3. Confirm that the required input files are present in `Data/`.
4. Run the scripts in `Scripts/`.
5. Figures will be written to `Outputs/figures/` and tables to `Outputs/tables/`.

The scripts are also designed to be run individually from an RStudio project session.

### Input files required

The following files must be present in `Data/` for the complete workflow:

- `Data-RSalmonella.RData`
- `Data_PCA2.snp`
- `Supplementary_Table_S1_metadata.csv`
- `Supplementary_Table_S2_Source_distribution.csv`
- `Supplementary_Table_S3_Genomes_used_in_genomic_analysis.csv`
- `snp_tree_2023-09-04.nwk`

`Data_PCA2.snp` is included in this project package. The remaining project data files should be copied from the final local second-revision folder before uploading the repository. This package cannot include local files that were not provided in the upload.

### Main analyses

- PCA of SNP genotypes: `PCA2_raw_rw.R`
- Circular phylogenetic tree: `Phylogenetic_tree_Salmonella_wild_birds_rw.R`
- Rectangular dendrogram: `Dendr_rect_nodos_or.R`
- Dendrogram with metadata heatmap: `Dendrog_and_heatmap.R`
- Genome/tree figure: `Genome_figure.R`
- Mean pairwise SNP-distance tables and heatmap: `Table_and_SNPdistance_matrix.R`

## Data and output policy

The repository should contain the input tables and files required to reproduce the reported analyses, together with the scripts. Personal R session files (`.Rhistory`, `.RData`, `.Rproj.user/`) are excluded by `.gitignore` unless a specific data file is intentionally renamed and placed in `Data/`.

## Manuscript repository

GitHub repository: https://github.com/Nivia-L/Salmonella_wild_birds
