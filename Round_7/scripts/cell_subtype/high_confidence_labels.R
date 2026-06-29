library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(shiny)
library(plotly)
library(ggplot2)
library(data.table)
library(FNN)
library(qs)
library(presto)
library(UpSetR)
library(pheatmap)
library(gridExtra)
library(shiny)
library(ggplot2)
library(pheatmap)

group_levels <- c(
  # Adult
  "Adult Control",
  "Adult RSV infected",
  "Adult reinfected",
  
  # Neonate - baseline
  "Neonate control",
  
  # Neonate - IFN upon infection
  "Neonate RSV infected (NO IFN)",
  "Neonate IFN and RSV infected",
  
  # Neonate - IFN upon reinfection
  "Neonate (NO IFN) reinfected",
  "Neonate IFN and RSV reinfected"
)


# Gene names matching the 1k mouse gene panel - from Diego
cell_types_markers_diego <- list(
  "Dendritic_Cells" = list(
    "Plasmacytoid_dendritic_cells_(pDCs)" = list(
      "Markers" = c("Itgax", "Bst2", "Siglech") # CD11c, BST2, Siglec-H
    ),
    "Conventional_dendritic_cells_(cDCs)" = list(
      "cDC1" = c("Itgax", "Cd8a"),   # CD11c, CD8a
      "cDC2" = c("Itgax", "Itgam")   # CD11c, CD11b
    )
  ),
  "Lymphocytes" = list(
    "T_CD4_Cells" = list(
      "Markers" = c("Ptprc", "Cd3e", "Cd4") # CD45, CD3, CD4
    ),
    "T_CD8_Cells" = list(
      "Markers" = c("Ptprc", "Cd3e", "Cd8a") # CD45, CD3, CD8a
    ),
    "Activated_B_Cells" = list(
      "Markers" = c("Ptprc", "Cd19", "Cd69"), # B220=CD45, CD19, CD69
      "Ligands_on_pDC" = c("Tnfsf13b", "Tnfsf13"), # BAFF, APRIL
      "Receptors_on_B_cells" = c("Tnfrsf13c", "Tnfrsf13b", "Tnfrsf17"), # BAFF-R, TACI, BCMA
      "Memory_B_cells" = "Cd27",
      "Plasmablast_and_Plasma_Cells" = "Sdc1" # CD138
    )
  ),
  "Epithelial_Cells" = list(
    "General_Marker" = "Epcam",
    "Subtypes" = list(
      "Club_Cells" = c("Scgb1a1"), # CC10/SCGB1A1
      "AT1" = c("Pdpn", "Cav1"), # podoplanin, caveolin
      "AT2" = "Sftpc",
      "Goblet_cells" = "Muc5ac",
      "Ciliated_cells" = c("Tuba1a", "Foxj1"), # α-Ac-Tub = Tuba1a (alpha-tubulin)
      "Basal_cells" = "Krt5",
      "Neuroendocrine_cells" = c("Calca", "Chga") # CGRP, Chromogranin A
    ),
    "ILCs" = list(
      "Markers" = c("Ptprc", "Il7r", "Cd3e"), # CD45, CD127, CD3 (lineage negative = Cd3e-)
      "ILC1" = NULL,
      "ILC2" = c("Il1rl1"), # ST2
      "ILC3" = NULL
    )
  )
)


source("/mnt/scratch2/Maycon/RuiPatricia_Medullo/Round_2/scripts/utils/utils.R")


run_tma_clustering <- function(
  obj,
  assay = "Nanostring",
  dims_use = 1:30,
  resolutions = c(0.3, 0.5),
  scale_factor = 1000,
  reduction_name = "PCA_regular",
  umap_name = "UMAP_regular",
  graph_name = "regular_snn",
  cluster_prefix = "cluster_regular_res",
  diet = TRUE,
  verbose = TRUE
) {
  DefaultAssay(obj) <- assay

  if (diet) {
    obj <- DietSeurat(
      obj,
      assays = assay,
      counts = TRUE,
      data = FALSE,
      scale.data = FALSE,
      dimreducs = NULL,
      graphs = NULL
    )
  }

  DefaultAssay(obj) <- assay

  obj <- NormalizeData(obj, assay = assay, scale.factor = scale_factor, verbose = verbose)

  features_use <- rownames(obj[[assay]])
  VariableFeatures(obj) <- features_use

  obj <- ScaleData(obj, assay = assay, features = features_use, verbose = verbose)

  obj <- RunPCA(
    obj,
    assay = assay,
    features = features_use,
    npcs = max(dims_use),
    reduction.name = reduction_name,
    verbose = verbose
  )

  obj <- FindNeighbors(
    obj,
    reduction = reduction_name,
    dims = dims_use,
    graph.name = graph_name,
    verbose = verbose
  )

  for (res in resolutions) {
    res_name <- gsub("\\.", "", as.character(res))
    cluster_col <- paste0(cluster_prefix, res_name)

    obj <- FindClusters(
      obj,
      graph.name = graph_name,
      resolution = res,
      cluster.name = cluster_col,
      verbose = verbose
    )
  }

  obj <- RunUMAP(
    obj,
    reduction = reduction_name,
    dims = dims_use,
    reduction.name = umap_name,
    verbose = verbose
  )

  return(obj)
}


get_top_markers <- function(
  obj,
  cluster_col,
  assay = "Nanostring",
  top_n = 10,
  logfc_min = 0,
  padj_max = 0.05
) {
  obj@assays$RNA <- obj@assays[[assay]]

  markers <- wilcoxauc(obj, cluster_col)

  top_markers <- markers %>%
    group_by(group) %>%
    dplyr::filter(logFC > logfc_min, padj <= padj_max) %>%
    dplyr::top_n(wt = logFC, n = top_n) %>%
    ungroup()

  return(list(
    markers = markers,
    top_markers = top_markers
  ))
}


plot_marker_dotplot <- function(
  obj,
  cluster_col,
  features,
  dot_shape = 21,
  dot_border_color = "black",
  dot_border_width = 0.5,
  color_palette = "magma",
  color_direction = 1,
  min_dot_size = 0,
  max_dot_size = 6,
  x_text_angle = 45,
  x_text_hjust = 1,
  x_text_vjust = 1,
  y_text_size = 11,
  x_text_size = 11,
  legend_title_size = 11,
  legend_text_size = 10,
  plot_title = "",
  base_size = 12,
  theme_base = ggplot2::theme_classic()
) {
  Idents(obj) <- cluster_col

  DotPlot(
    obj,
    cluster.idents = TRUE,
    features = unique(features)
  ) +
    geom_point(
      aes(size = pct.exp),
      shape = dot_shape,
      colour = dot_border_color,
      stroke = dot_border_width
    ) +
    scale_colour_viridis_c(
      option = color_palette,
      direction = color_direction
    ) +
    scale_size(
      range = c(min_dot_size, max_dot_size)
    ) +
    guides(
      size = guide_legend(
        override.aes = list(
          shape = dot_shape,
          colour = dot_border_color,
          fill = "white",
          stroke = dot_border_width
        )
      )
    ) +
    theme_base +
    theme(
      axis.text.x = element_text(
        angle = x_text_angle,
        hjust = x_text_hjust,
        vjust = x_text_vjust,
        size = x_text_size
      ),
      axis.text.y = element_text(size = y_text_size),
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size),
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle(plot_title)
}


### Load data ------
SeuObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/data/Processed/merged_Sobj_UMAP.qs")
SeuObj@assays$niche <- NULL

### Lymphocytes - subtypes (Tcd4, Tcd8, Bcell) ---------

# Extract raw counts matrix from Nanostring assay
counts <- SeuObj@assays$Nanostring$counts
gene_panel <- rownames(counts)
# Tcd8 markers rule 
gene_panel[grepl("Ptprc", gene_panel)]  #Ptprc
gene_panel[grepl("Cd3", gene_panel)] #Cd3d, Cd3e, Cd3g
gene_panel[grepl("Cd8", gene_panel)] #Cd8a, Cd8b1
# negative to Cd19

# Tcd4 markers rule 
gene_panel[grepl("Ptprc", gene_panel)]  #Ptprc
gene_panel[grepl("Cd3", gene_panel)] #Cd3d, Cd3e, Cd3g
gene_panel[grepl("Cd4", gene_panel)] #Cd4
# negative to Cd19

# Bcell markers rule
gene_panel[grepl("Ptprc", gene_panel)]  #Ptprc
gene_panel[grepl("Cd19", gene_panel)] #Cd19
# negative to Cd3


## Loose call: any of the subunits, exp >= 1
# Marker sets
ptprc_marker <- "Ptprc"
cd3_markers  <- c("Cd3d", "Cd3e", "Cd3g")
cd8_markers  <- c("Cd8a", "Cd8b1")
cd4_markers  <- c("Cd4")
b_markers    <- c("Cd19")

# Keep only markers present in panel
cd3_markers <- intersect(cd3_markers, gene_panel)

cd8_markers <- intersect(cd8_markers, gene_panel)

cd4_markers <- intersect(cd4_markers, gene_panel)

b_markers   <- intersect(b_markers, gene_panel)

# Helper: TRUE if any marker in set has counts > 0

has_any <- function(markers, counts) {

  if (length(markers) == 0) {

    return(rep(FALSE, ncol(counts)))

  }

  Matrix::colSums(counts[markers, , drop = FALSE] > 0) > 0

}

# Marker positivity
has_Ptprc <- counts["Ptprc", ] >= 1

has_Cd3   <- has_any(cd3_markers, counts)

has_Cd8   <- has_any(cd8_markers, counts)

has_Cd4   <- has_any(cd4_markers, counts)

has_Cd19  <- has_any(b_markers, counts)

# Initialize annotation

annot <- rep("other", ncol(counts))

# Annotation rules

annot[has_Ptprc & has_Cd3 & has_Cd4 & !has_Cd19 & !has_Cd8] <- "Tcd4"

annot[has_Ptprc & has_Cd3 & has_Cd8 & !has_Cd19 & !has_Cd4] <- "Tcd8"

annot[has_Ptprc & has_Cd19 & !has_Cd3] <- "Bcell"


# Optional: flag ambiguous T cells
#annot[has_Ptprc & has_Cd3 & has_Cd4 & has_Cd8 & !has_Cd19] <- "Tcell_CD4_CD8_double_positive"

# Add annotation to Seurat metadata
SeuObj@meta.data$lyn_minor_celltype <- annot

# Percent within each lyn_minor_celltype row
tab <- table(SeuObj$lyn_minor_celltype, SeuObj$major_celltype)
tab_pct <- prop.table(tab, margin = 1) * 100

pheatmap::pheatmap(
  tab_pct,
  scale = "none",
  cluster_rows = T,
  cluster_cols = T,
  display_numbers = T,
  number_format = "%.1f",
  main = "% of cells per lyn_minor_celltype"
)

table(annot)
annot_1 <- annot


## Strict call: all of the subunits, exp >= 1
# Marker sets
ptprc_marker <- "Ptprc"
cd3_markers  <- c("Cd3d", "Cd3e", "Cd3g")
cd8_markers  <- c("Cd8a", "Cd8b1")
cd4_markers  <- c("Cd4")
b_markers    <- c("Cd19")

# Require markers to exist in gene panel
required_markers <- c(ptprc_marker, cd3_markers, cd8_markers, cd4_markers, b_markers)
missing_markers <- setdiff(required_markers, gene_panel)

if (length(missing_markers) > 0) {
  stop("Missing required markers: ", paste(missing_markers, collapse = ", "))
}

# Helper: TRUE only if all markers are expressed >= 1 transcript
has_all <- function(markers, counts) {
  Matrix::colSums(counts[markers, , drop = FALSE] >= 1) == length(markers)
}

# Marker positivity: concomitant expression
has_Ptprc <- counts["Ptprc", ] >= 1
has_Cd3   <- has_all(cd3_markers, counts)
has_Cd8   <- has_all(cd8_markers, counts)
has_Cd4   <- has_all(cd4_markers, counts)
has_Cd19  <- has_all(b_markers, counts)

# Initialize annotation
annot <- rep("other", ncol(counts))

# Annotation rules
annot[has_Ptprc & has_Cd3 & has_Cd4 & !has_Cd19 & !has_Cd8] <- "Tcd4"

annot[has_Ptprc & has_Cd3 & has_Cd8 & !has_Cd19 & !has_Cd4] <- "Tcd8"

annot[has_Ptprc & has_Cd19 & !has_Cd3] <- "Bcell"

# Add annotation
SeuObj@meta.data$lyn_minor_celltype <- annot

# Percent within each lyn_minor_celltype row
tab <- table(SeuObj$lyn_minor_celltype, SeuObj$major_celltype)
tab_pct <- prop.table(tab, margin = 1) * 100

pheatmap::pheatmap(
  tab_pct,
  scale = "none",
  cluster_rows = T,
  cluster_cols = T,
  display_numbers = T,
  number_format = "%.1f",
  main = "% of cells per lyn_minor_celltype"
)

table(annot)
annot_2 <- annot

SeuObj@meta.data$lyn_minor_celltype_LOOSE <-  annot_1
SeuObj@meta.data$lyn_minor_celltype_STRICT <-  annot_2


obj <- subset(SeuObj, lyn_minor_celltype_LOOSE %in% c("Bcell", "Tcd4", "Tcd8") | lyn_minor_celltype_STRICT %in% c("Bcell", "Tcd4", "Tcd8"))

assay <- "Nanostring"
scale_factor <- 1000
obj <- run_tma_clustering(
  obj = obj,
  assay = assay,
  dims_use = 1:15,
  resolutions = c(0.3, 0.5),
  scale_factor = scale_factor
)

# lyn_color <- c("Bcell" = "#F0E442", 
#                    "Tcd4" = "#E69F22", 
#                    "Tcd8" = "#D55E00", 
#                    "other" = "black")
lyn_color <- c("Bcell" = "gold",
                   "Tcd4" = "blue4",
                   "Tcd8" = "#D55E90",
                   "other" = "black")

DimPlot(obj, reduction = "UMAP_regular", group.by = "cluster_regular_res05") |
DimPlot(subset(obj, lyn_minor_celltype_LOOSE %in% "other", invert = T), 
        reduction = "UMAP_regular", 
        group.by = "lyn_minor_celltype_LOOSE", 
        cols = lyn_color,
        pt.size = 0.5
) |
DimPlot(subset(obj, lyn_minor_celltype_STRICT %in% "other", invert = T), 
        reduction = "UMAP_regular",
        group.by = "lyn_minor_celltype_STRICT",
        cols = lyn_color,
        pt.size = 0.5)

library(viridis)
FeaturePlot(
  subset(obj, lyn_minor_celltype_STRICT %in% "other", invert = TRUE),
  reduction = "UMAP_regular",
  features = c("nCount", "nFeature"),
  max.cutoff = "q99",
  order = TRUE,
  pt.size = 0.5
) & 
  scale_color_viridis(option = "viridis")




library(UCell)
geneSets <- list(
  Tcd4  = c("Cd3d", "Cd3e", "Cd3g", "Cd4"),
  Tcd8  = c("Cd3d", "Cd3e", "Cd3g", "Cd8a", "Cd8b1"),
  Bcell = c("Cd19", "Ms4a1", "Cd79a")
)

# Keep only genes present
geneSets <- lapply(geneSets, function(gs) intersect(gs, rownames(obj)))

obj <- AddModuleScore_UCell(
  obj = obj,
  #assay = "Nanostring",
  features = geneSets
)

library(viridis)
FeaturePlot(
  obj,
  #subset(obj, lyn_minor_celltype_STRICT %in% "other", invert = TRUE),
  reduction = "UMAP_regular",
  features = c("Tcd4_UCell", "Tcd8_UCell", "Bcell_UCell"),
  max.cutoff = "q99",
  order = TRUE,
  pt.size = 0.5
) & 
  scale_color_viridis(option = "viridis")





obj <- subset(SeuObj, lyn_minor_celltype_LOOSE %in% c("Tcd8", "Tcd4"))

assay <- "Nanostring"
scale_factor <- 1000
obj <- run_tma_clustering(
  obj = obj,
  assay = assay,
  dims_use = 1:50,
  resolutions = c(0.3, 0.5),
  scale_factor = scale_factor
)

DimPlot(obj, reduction = "UMAP_regular", group.by = "cluster_regular_res05") |
DimPlot(subset(obj, lyn_minor_celltype_LOOSE %in% "other", invert = T), 
        reduction = "UMAP_regular", 
        group.by = "lyn_minor_celltype_LOOSE", 
        cols = lyn_color,
        pt.size = 0.5
) 
  
  
  

# Keep only genes actually present in your panel
immune_celltype_genes <- c(
  # Pan-immune / leukocyte
  "Ptprc", "Ptprcap",

  # T cells
  "Cd2", "Cd3d", "Cd3e", "Cd3g", "Cd5", "Cd6",
  "Cd4", "Cd8a", "Cd8b1",
  "Trac", "Cd28", "Icos", "Il7r", "Sell", "Tcf7",

  # T cell activation / exhaustion / regulatory
  "Cd44", "Cd69", "Ctla4", "Pdcd1", "Lag3", "Tigit",
  "Foxp3", "Il2ra", "Havcr2",

  # Cytotoxic T / NK
  "Nkg7", "Prf1", "Gzma", "Gzmb", "Gzmk",
  "Ncr1", "Klrb1c", "Klrk1", "Eomes", "Tbx21",

  # B cells
  "Cd19", "Ms4a1", "Cd79a", "Cd22", "Cd37",
  "Pax5", "Blk", "Btk", "Cr2", "Ighm", "Ighd",
  "Igkc", "Fcrla", "Vpreb3",

  # Plasma cells
  "Sdc1", "Mzb1", "Jchain", "Xbp1", "Prdm1",

  # Myeloid / macrophage / monocyte
  "Cd14", "Cd68", "Adgre5", "Itgam", "Csf1r",
  "Lyz3", "Aif1", "C1qa", "C1qb", "C1qc",
  "Apoe", "Mrc1", "Marco", "Mertk", "Ms4a4a",
  "Tyrobp", "Fcgr4", "Fcer1g",

  # Dendritic cells / antigen presentation
  "Itgax", "Flt3", "Clec10a", "Clec12a", "Clec4a1",
  "Cd209a", "Cd209e", "Lamp3", "Ciita",
  "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-DMa", "H2-DMb1", "H2-DMb2",
  "Cd74",

  # Neutrophils / granulocytes
  "S100a8", "S100a9", "Mpo", "Elane", "Prtn3",
  "Csf3r", "Fpr1", "Camp",

  # Mast cells
  "Kit", "Cpa3", "Tpsab1", "Tpsb2", "Mcpt8",

  # General immune signaling / chemokine receptors
  "Ccr2", "Ccr5", "Ccr7", "Cxcr3", "Cxcr4", "Cxcr5",
  "Cx3cr1", "Il2rg", "Il7r", "Ifng", "Tnf"
)

immune_celltype_genes <- intersect(immune_celltype_genes, gene_panel)

obj <- subset(SeuObj, lyn_minor_celltype_LOOSE %in% c("Bcell", "Tcd4", "Tcd8") | lyn_minor_celltype_STRICT %in% c("Bcell", "Tcd4", "Tcd8"))
obj <- subset(
  obj,
  features = immune_celltype_genes
)


assay <- "Nanostring"
scale_factor <- 1000
obj <- run_tma_clustering(
  obj = obj,
  assay = assay,
  dims_use = 1:15,
  resolutions = c(0.3, 0.5),
  scale_factor = scale_factor
)

# lyn_color <- c("Bcell" = "#F0E442", 
#                    "Tcd4" = "#E69F22", 
#                    "Tcd8" = "#D55E00", 
#                    "other" = "black")
lyn_color <- c("Bcell" = "gold",
                   "Tcd4" = "blue4",
                   "Tcd8" = "#D55E90",
                   "other" = "black")

DimPlot(obj, reduction = "UMAP_regular", group.by = "cluster_regular_res05") |
DimPlot(subset(obj, lyn_minor_celltype_LOOSE %in% "other", invert = T), 
        reduction = "UMAP_regular", 
        group.by = "lyn_minor_celltype_LOOSE", 
        cols = lyn_color,
        pt.size = 0.5
) |
DimPlot(subset(obj, lyn_minor_celltype_STRICT %in% "other", invert = T), 
        reduction = "UMAP_regular",
        group.by = "lyn_minor_celltype_STRICT",
        cols = lyn_color,
        pt.size = 0.5)

library(viridis)
FeaturePlot(
  subset(obj, lyn_minor_celltype_STRICT %in% "other", invert = TRUE),
  reduction = "UMAP_regular",
  features = c("nCount", "nFeature"),
  max.cutoff = "q99",
  order = TRUE,
  pt.size = 0.5
) & 
  scale_color_viridis(option = "viridis")











