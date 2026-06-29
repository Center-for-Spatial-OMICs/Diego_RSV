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


source("/mnt/scratch2/Maycon/Utils/R_codabase/Utils.R")


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


### Decided on  -------------------------------------------------------
## 1. Let's keep only these cell types 
# Dendritic cells: Cd74; H2-Ab1; B2m; Ppia; Ciita; Itgb2
# Epithelial cells: Tubb4b, Sod1, Cd24a, Hsp90aa1, Krt8, Cd9
# B Cells: Jchain; Igkc; Ighm; Ly6d; Cd37; Blk
# T cells: Cd274; Cd4; Gata3; Tnfrsf9; Ly75; Cd59a
# Myeloid cells (Monocyte progenitor): Elane; Ctsg; Srgn; Cd1d1; Cd48; Clec10a
# Endothelial cells: Cdh5; Tie1; Esam; Pecam1; Cd93; Icam2

## 2. For major cell types
# Using a cluster-based approach 
# Lymphocytes, Myeloid, Dendritic, Epithelial, Endothelial

## 3. For subcell-type (Tcells, Bcells, pDC, cDC, epithelials ...)
# Using a cell-based approach 
# Let's use all sub-units as "any" - so any cell that express at least one sub-unit (eg: Cd3d and Cd3e) would be called a Tcd4 cell.


## Adding/deriving columns -------
SeuObj$Group_2 <- SeuObj$Group
SeuObj$Group_2[
  SeuObj$Group %in% c(
    "Adult reinfected",
    "Adult RSV infected"
  )
] <- "Adult"

SeuObj$Group_2[
  SeuObj$Group %in% c(
    "Neonate IFN and RSV infected",
    "Neonate IFN and RSV reinfected",
    "Neonate RSV infected (NO IFN)",
    "Neonate (NO IFN) reinfected"
  )
] <- "Neonate RSV"

SeuObj$Group_2[
  SeuObj$Group %in% c(
    "Adult Control",
    "Neonate control"
  )
] <- "Neonate"


# Removing "cells not in any of the FOV maps I got from Hannah 
SeuObj@meta.data[SeuObj@meta.data$TMA_4 %in% "NA", ]
SeuObj@meta.data[SeuObj@meta.data$TMA_4 %in% "NA.1", ]

SeuObj <- subset(SeuObj,TMA_4 %in% c("NA", "NA.1"), invert = T)


## Check current QC state --------
qc_results <- run_qc_and_plot(
    seurat_obj = SeuObj,
    assay = "Nanostring",
    sample_col = "TMA_4",
    group_col = "Group",
    fill_col = "Group_2",
    qc_col = "qcFlag",
    ncount_thresh = 60,
    nfeature_thresh = 30, 
    bins = 100,
    ncount_xlim = c(0, 1000),
    nfeature_xlim = c(0, 1000)
)

# qc_results$qc_plot
# qc_results$ncount_plot
# qc_results$nfeature_plot
# qc_results$bin_plot

SeuObj <- qc_results$seurat_obj


# Note: let's keep all cells regardless QC - we can layover the QC later when deciding on cell types. That way we can be more certain whether we can toss a TMA


### Cluster-based annotation --------
assay <- "Nanostring"
scale_factor <- 1000


obj <- run_tma_clustering(
  obj = SeuObj,
  assay = assay,
  dims_use = 1:30,
  resolutions = c(0.3, 0.6),
  scale_factor = scale_factor
)


#qsave(obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_all467158cells_clustered.qs")

SeuObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_all467158cells_clustered.qs")


cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- SeuObj$cluster_regular_res03 %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)


DimPlot(SeuObj,
        reduction = "UMAP_regular",
        group.by = "cluster_regular_res03", cols = celltype_cols)

mk_general <- get_top_markers(
  obj = SeuObj,
  cluster_col = "cluster_regular_res03",
  assay = assay,
  top_n = 20
)

# top cluster markers
plot_marker_dotplot(
  SeuObj,
  cluster_col = "cluster_regular_res03",
  features = unique(mk_general$top_markers$feature),
  dot_shape = 21,
  dot_border_color = "black",
  dot_border_width = 0.5,
  color_palette = "magma",
  color_direction = 1,
  min_dot_size = 0,
  max_dot_size = 10,
  x_text_angle = 90,
  x_text_hjust = 1,
  x_text_vjust = 1,
  y_text_size = 22,
  x_text_size = 16,
  legend_title_size = 11,
  legend_text_size = 10,
  plot_title = "",
  base_size = 12,
  theme_base = ggplot2::theme_classic()
)

write.csv(
  data.frame(feature = mk_general$top_markers),
  file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/top_cluster_markers.csv",
  row.names = FALSE
)


cluster_to_celltype <- list(
  `0`  = "Lymphocytes", # Cd274, Cd4, Ly75
  `1`  = "Dendritic",   # Cd74, H2-Ab1, B2m
  `2`  = "Endothelial", # Cdh5, Tie1, Esam, Pecam1, Cd93, Icam2
  `3`  = "Myeloid",     # weak/unclear from provided markers; no strong anchor detected
  `4`  = "Dendritic",   # Cd74, H2-Ab1, Ppia
  `5`  = "Epithelial",  # Krt7/Krt19-like epithelial signal; weak match to provided marker list
  `6`  = "Epithelial",  # Tubb4b, Sod1, Cd24a, Hsp90aa1, Krt8, Cd9
  `7`  = "Lymphocytes", # Igkc, Ighm, Cd37
  `8`  = "Myeloid",     # Ptprc, Itgb2, Tyrobp, S100a8, S100a9
  `9`  = "Epithelial",  # weak/unclear from provided categories; no strong anchor detected
  `10` = "Lymphocytes", # Tnfrsf9
  `11` = "Endothelial"  # Flt1, Efnb2, Cxcl12, vascular/perivascular-like signal
)


Idents(SeuObj) <- "cluster_regular_res03"
SeuObj@meta.data$major_celltype <- cluster_to_celltype[as.character(Idents(SeuObj))]
table(SeuObj$major_celltype)
SeuObj$major_celltype <- as.character(unlist(SeuObj$major_celltype))


cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- SeuObj$major_celltype %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)


DimPlot(SeuObj,
        reduction = "UMAP_regular",
        group.by = "major_celltype", 
        cols = celltype_cols)


### Cell-based annotation --------
library(Seurat)
library(Matrix)

## lyn_minor_celltype (Tcd4, Tcd8, Bcell)
# Marker definitions
marker_list <- list(
  Ptprc = c("Ptprc"),
  Cd3   = c("Cd3d", "Cd3e", "Cd3g"),
  Cd8   = c("Cd8a", "Cd8b1"),
  Cd4   = c("Cd4"),
  Cd19  = c("Cd19")
)

# Annotation rules
annotation_rules <- list(

  Tcd4 = expression(
    Ptprc & Cd3 & Cd4 & !Cd19 & !Cd8
  ),

  Tcd8 = expression(
    Ptprc & Cd3 & Cd8 & !Cd19 & !Cd4
  ),

  Bcell = expression(
    Ptprc & Cd19 & !Cd3
  )

)

# Run annotation
SeuObj <- cell_based_anno(
  seurat_obj = SeuObj,
  assay = "Nanostring",
  slot = "counts",
  mode = "any",  # choose: "any" or "all"

  anno_column = c("lyn_minor_celltype__cell_based_any"),

  marker_list = marker_list,

  annotation_rules = annotation_rules,

  # # Optional custom thresholds
  # thresholds = list(
  #   Ptprc = 2,
  #   Cd4   = 1,
  #   Cd19  = 1 ....
  # ),

  default_threshold = 1
)


SeuObj <- cell_based_anno(
  seurat_obj = SeuObj,
  assay = "Nanostring",
  slot = "counts",
  mode = "all",  # choose: "any" or "all"

  anno_column = c("lyn_minor_celltype__cell_based_all"),

  marker_list = marker_list,

  annotation_rules = annotation_rules,

  # # Optional custom thresholds
  # thresholds = list(
  #   Ptprc = 2,
  #   Cd4   = 1,
  #   Cd19  = 1 ....
  # ),

  default_threshold = 1
)


## Cell Type Annotation Agreement (Jaccard %) Across annotations
library(pheatmap)

# Build annotation dataframe
anno_df <- data.frame(
  cell_id = colnames(SeuObj),
  anno_1 = SeuObj$lyn_minor_celltype__cell_based_any,
  anno_2 = SeuObj$lyn_minor_celltype__cell_based_all,
  stringsAsFactors = FALSE
)

# Remove NA annotations
anno_df <- anno_df[
  !is.na(anno_df$anno_1) &
    !is.na(anno_df$anno_2),
]

# Confusion table
tab <- table(
  anno_df$anno_1,
  anno_df$anno_2
)

# Jaccard similarity matrix
jaccard <- matrix(
  0,
  nrow = nrow(tab),
  ncol = ncol(tab),
  dimnames = dimnames(tab)
)

for (i in rownames(tab)) {

  for (j in colnames(tab)) {

    intersection <- tab[i, j]

    union <- (
      sum(tab[i, ]) +
      sum(tab[, j]) -
      intersection
    )

    jaccard[i, j] <- intersection / union * 100

  }
}

# Add label sizes
rownames(jaccard) <- paste0(
  rownames(jaccard),
  " (n=",
  rowSums(tab),
  ")"
)

colnames(jaccard) <- paste0(
  colnames(jaccard),
  " (n=",
  colSums(tab),
  ")"
)

# Plot
pheatmap(
  mat = jaccard,

  color = colorRampPalette(
    c("white", "gold", "darkred")
  )(100),

  border_color = "grey90",

  cluster_rows = TRUE,
  cluster_cols = TRUE,

  display_numbers = round(jaccard, 1),

  number_color = "black",
  fontsize_number = 10,

  fontsize_row = 12,
  fontsize_col = 12,

  angle_col = 45,

  main = "Cell Type Annotation Agreement (Jaccard %)"
)


### Cell-based annotation under Cluster-based annotation --------
library(Seurat)
library(Matrix)

## lyn_minor_celltype (Tcd4, Tcd8, Bcell)
# Marker definitions
marker_list <- list(
  Ptprc = c("Ptprc"),
  Cd3   = c("Cd3d", "Cd3e", "Cd3g"),
  Cd8   = c("Cd8a", "Cd8b1"),
  Cd4   = c("Cd4"),
  Cd19  = c("Cd19")
)

# Annotation rules
annotation_rules <- list(

  Tcd4 = expression(
    Ptprc & Cd3 & Cd4 & !Cd19 & !Cd8
  ),

  Tcd8 = expression(
    Ptprc & Cd3 & Cd8 & !Cd19 & !Cd4
  ),

  Bcell = expression(
    Ptprc & Cd19 & !Cd3
  )

)

obj <- subset(SeuObj, major_celltype %in% "Lymphocytes")
# Run annotation
obj <- cell_based_anno(
  seurat_obj = obj,
  assay = "Nanostring",
  slot = "counts",
  mode = "all",  # choose: "any" or "all"

  anno_column = c("lyn_minor_celltype__cell_based_any__cluster_gated"),

  marker_list = marker_list,

  annotation_rules = annotation_rules,

  # # Optional custom thresholds
  # thresholds = list(
  #   Ptprc = 2,
  #   Cd4   = 1,
  #   Cd19  = 1 ....
  # ),

  default_threshold = 1
)

# get anno back to the main obj 
SeuObj$lyn_minor_celltype__cell_based_any__cluster_gated <- "other"
SeuObj@meta.data[
  colnames(obj),
  "lyn_minor_celltype__cell_based_any__cluster_gated"
] <- obj$lyn_minor_celltype__cell_based_any__cluster_gated


## Cell Type Annotation Agreement (Jaccard %) Across annotations
library(pheatmap)

# Build annotation dataframe
anno_df <- data.frame(
  cell_id = colnames(SeuObj),
  anno_1 = SeuObj$lyn_minor_celltype__cell_based_any,
  anno_2 = SeuObj$lyn_minor_celltype__cell_based_any__cluster_gated,
  stringsAsFactors = FALSE
)

# Remove NA annotations
anno_df <- anno_df[
  !is.na(anno_df$anno_1) &
    !is.na(anno_df$anno_2),
]

# Confusion table
tab <- table(
  anno_df$anno_1,
  anno_df$anno_2
)

# Jaccard similarity matrix
jaccard <- matrix(
  0,
  nrow = nrow(tab),
  ncol = ncol(tab),
  dimnames = dimnames(tab)
)

for (i in rownames(tab)) {

  for (j in colnames(tab)) {

    intersection <- tab[i, j]

    union <- (
      sum(tab[i, ]) +
      sum(tab[, j]) -
      intersection
    )

    jaccard[i, j] <- intersection / union * 100

  }
}

# Add label sizes
rownames(jaccard) <- paste0(
  rownames(jaccard),
  " (n=",
  rowSums(tab),
  ")"
)

colnames(jaccard) <- paste0(
  colnames(jaccard),
  " (n=",
  colSums(tab),
  ")"
)

# Plot
pheatmap(
  mat = jaccard,

  color = colorRampPalette(
    c("white", "gold", "darkred")
  )(100),

  border_color = "grey90",

  cluster_rows = TRUE,
  cluster_cols = TRUE,

  display_numbers = round(jaccard, 1),

  number_color = "black",
  fontsize_number = 10,

  fontsize_row = 12,
  fontsize_col = 12,

  angle_col = 45,

  main = "Cell Type Annotation Agreement (Jaccard %)"
)



## Cell Type Annotation Agreement (Jaccard %) Across annotations
library(pheatmap)

# Build annotation dataframe
anno_df <- data.frame(
  cell_id = colnames(SeuObj),
  anno_1 = SeuObj$lyn_minor_celltype__cell_based_any,
  anno_2 = SeuObj$major_celltype,
  stringsAsFactors = FALSE
)

# Remove NA annotations
anno_df <- anno_df[
  !is.na(anno_df$anno_1) &
    !is.na(anno_df$anno_2),
]

# Confusion table
tab <- table(
  anno_df$anno_1,
  anno_df$anno_2
)

# Jaccard similarity matrix
jaccard <- matrix(
  0,
  nrow = nrow(tab),
  ncol = ncol(tab),
  dimnames = dimnames(tab)
)

for (i in rownames(tab)) {

  for (j in colnames(tab)) {

    intersection <- tab[i, j]

    union <- (
      sum(tab[i, ]) +
      sum(tab[, j]) -
      intersection
    )

    jaccard[i, j] <- intersection / union * 100

  }
}

# Add label sizes
rownames(jaccard) <- paste0(
  rownames(jaccard),
  " (n=",
  rowSums(tab),
  ")"
)

colnames(jaccard) <- paste0(
  colnames(jaccard),
  " (n=",
  colSums(tab),
  ")"
)

# Plot
pheatmap(
  mat = jaccard,

  color = colorRampPalette(
    c("white", "gold", "darkred")
  )(100),

  border_color = "grey90",

  cluster_rows = TRUE,
  cluster_cols = TRUE,

  display_numbers = round(jaccard, 1),

  number_color = "black",
  fontsize_number = 10,

  fontsize_row = 12,
  fontsize_col = 12,

  angle_col = 45,

  main = "Cell Type Annotation Agreement (Jaccard %)"
)







## Cell Type Annotation Agreement (Jaccard %) Across annotations
library(pheatmap)

# Build annotation dataframe
anno_df <- data.frame(
  cell_id = colnames(SeuObj),
  anno_1 = SeuObj$lyn_minor_celltype__cell_based_any,
  anno_2 = SeuObj$lyn_minor_celltype__cell_based_any__cluster_gated,
  stringsAsFactors = FALSE
)

# Remove NA annotations
anno_df <- anno_df[
  !is.na(anno_df$anno_1) &
    !is.na(anno_df$anno_2),
]

# Confusion table
tab <- table(
  anno_df$anno_1,
  anno_df$anno_2
)

# Jaccard similarity matrix
jaccard <- matrix(
  0,
  nrow = nrow(tab),
  ncol = ncol(tab),
  dimnames = dimnames(tab)
)

for (i in rownames(tab)) {

  for (j in colnames(tab)) {

    intersection <- tab[i, j]

    union <- (
      sum(tab[i, ]) +
      sum(tab[, j]) -
      intersection
    )

    jaccard[i, j] <- intersection / union * 100

  }
}

# Add label sizes
rownames(jaccard) <- paste0(
  rownames(jaccard),
  " (n=",
  rowSums(tab),
  ")"
)

colnames(jaccard) <- paste0(
  colnames(jaccard),
  " (n=",
  colSums(tab),
  ")"
)

# Plot
pheatmap(
  mat = jaccard,

  color = colorRampPalette(
    c("white", "gold", "darkred")
  )(100),

  border_color = "grey90",

  cluster_rows = TRUE,
  cluster_cols = TRUE,

  display_numbers = round(jaccard, 1),

  number_color = "black",
  fontsize_number = 10,

  fontsize_row = 12,
  fontsize_col = 12,

  angle_col = 45,

  main = "Cell Type Annotation Agreement (Jaccard %)"
)



## Cell Type Annotation Agreement (Jaccard %) Across annotations
library(pheatmap)

# Build annotation dataframe
anno_df <- data.frame(
  cell_id = colnames(SeuObj),
  anno_1 = SeuObj$lyn_minor_celltype__cell_based_any__cluster_gated,
  anno_2 = SeuObj$major_celltype,
  stringsAsFactors = FALSE
)

# Remove NA annotations
anno_df <- anno_df[
  !is.na(anno_df$anno_1) &
    !is.na(anno_df$anno_2),
]

# Confusion table
tab <- table(
  anno_df$anno_1,
  anno_df$anno_2
)

# Jaccard similarity matrix
jaccard <- matrix(
  0,
  nrow = nrow(tab),
  ncol = ncol(tab),
  dimnames = dimnames(tab)
)

for (i in rownames(tab)) {

  for (j in colnames(tab)) {

    intersection <- tab[i, j]

    union <- (
      sum(tab[i, ]) +
      sum(tab[, j]) -
      intersection
    )

    jaccard[i, j] <- intersection / union * 100

  }
}

# Add label sizes
rownames(jaccard) <- paste0(
  rownames(jaccard),
  " (n=",
  rowSums(tab),
  ")"
)

colnames(jaccard) <- paste0(
  colnames(jaccard),
  " (n=",
  colSums(tab),
  ")"
)

# Plot
pheatmap(
  mat = jaccard,

  color = colorRampPalette(
    c("white", "gold", "darkred")
  )(100),

  border_color = "grey90",

  cluster_rows = TRUE,
  cluster_cols = TRUE,

  display_numbers = round(jaccard, 1),

  number_color = "black",
  fontsize_number = 10,

  fontsize_row = 12,
  fontsize_col = 12,

  angle_col = 45,

  main = "Cell Type Annotation Agreement (Jaccard %)"
)



## Cell Type Annotation Agreement (Jaccard %) Across annotations
# Evaluating QC ~ cell anno limits 


SeuObj$qcFlag_2 <- "Good quality"

SeuObj$qcFlag_2[
  SeuObj$nCount_RNA < 180 |
  SeuObj$nFeature_RNA < 90
] <- "Low quality"

table(SeuObj$qcFlag_2)


library(pheatmap)

anno_df <- data.frame(
  cell_id = colnames(SeuObj),
  anno_1 = SeuObj$lyn_minor_celltype__cell_based_any,
  anno_2 = SeuObj$major_celltype,
  qcFlag = SeuObj$qcFlag, # or qcFlag_2
  stringsAsFactors = FALSE
)

anno_df <- anno_df[
  !is.na(anno_df$anno_1) &
    !is.na(anno_df$anno_2) &
    !is.na(anno_df$qcFlag),
]

anno_df$anno_2_qc <- paste0(
  anno_df$anno_2,
  " | ",
  anno_df$qcFlag
)

tab <- table(
  anno_df$anno_1,
  anno_df$anno_2_qc
)

jaccard <- matrix(
  0,
  nrow = nrow(tab),
  ncol = ncol(tab),
  dimnames = dimnames(tab)
)

for (i in rownames(tab)) {
  for (j in colnames(tab)) {
    intersection <- tab[i, j]
    union <- sum(tab[i, ]) + sum(tab[, j]) - intersection

    jaccard[i, j] <- ifelse(
      union == 0,
      0,
      intersection / union * 100
    )
  }
}

col_info <- data.frame(
  colname = colnames(jaccard),
  celltype = sub(" \\| .*", "", colnames(jaccard)),
  qcFlag = sub(".*\\| ", "", colnames(jaccard)),
  stringsAsFactors = FALSE
)

qc_order <- c("Good quality", "Low quality")

col_info$qcFlag <- factor(
  col_info$qcFlag,
  levels = qc_order
)

col_info <- col_info[
  order(col_info$celltype, col_info$qcFlag),
]

jaccard <- jaccard[, col_info$colname, drop = FALSE]

row_labels <- paste0(
  rownames(jaccard),
  " (n=",
  rowSums(tab),
  ")"
)

col_labels <- paste0(
  col_info$colname,
  " (n=",
  colSums(tab)[col_info$colname],
  ")"
)

rownames(jaccard) <- row_labels
colnames(jaccard) <- col_labels

qc_annotation <- data.frame(
  qcFlag = col_info$qcFlag
)

rownames(qc_annotation) <- col_labels

annotation_colors <- list(
  qcFlag = c(
    "Good quality" = "forestgreen",
    "Low quality" = "firebrick"
  )
)

pheatmap(
  mat = jaccard,
  color = colorRampPalette(c("white", "gold", "darkred"))(100),
  border_color = "grey90",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = qc_annotation,
  annotation_colors = annotation_colors,
  display_numbers = round(jaccard, 1),
  number_color = "black",
  fontsize_number = 10,
  fontsize_row = 12,
  fontsize_col = 10,
  angle_col = 45,
  main = "Cell Type Annotation Agreement (Jaccard %) with QC"
)





library(mclust)
library(aricode)
library(dplyr)
library(ggplot2)

anno_df <- data.frame(
  cell_id = colnames(SeuObj),
  anno_1 = SeuObj$lyn_minor_celltype__cell_based_any,
  anno_2 = SeuObj$major_celltype,
  qcFlag = SeuObj$qcFlag_2,
  stringsAsFactors = FALSE
)

anno_df <- anno_df[
  complete.cases(anno_df),
]

anno_df <- anno_df[!anno_df$anno_1 %in% "other" , ]

# Overall contingency test
cont_table <- table(
  anno_df$anno_1,
  anno_df$anno_2,
  anno_df$qcFlag
)

chisq_result <- chisq.test(
  apply(cont_table, c(1,2), sum)
)

print(chisq_result)

# Metrics per QC group
qc_stats <- lapply(
  unique(anno_df$qcFlag),
  function(qc) {

    sub <- anno_df[
      anno_df$qcFlag == qc,
    ]

    ari <- adjustedRandIndex(
      sub$anno_1,
      sub$anno_2
    )

    nmi <- NMI(
      sub$anno_1,
      sub$anno_2
    )

    data.frame(
      qcFlag = qc,
      ARI = ari,
      NMI = nmi,
      n_cells = nrow(sub)
    )
  }
)

qc_stats <- bind_rows(qc_stats)

print(qc_stats)

# ARI ~ 0 or negative
# 
# Adjusted Rand Index:
# 
# * 1 = perfect agreement
# * 0 = random agreement
# * <0 = worse than random

# If ~ 0, it says that "qcFlag_2" is not the main driver of disagreement

# Per-cell-type Jaccard comparison
calc_jaccard <- function(df) {

  tab <- table(df$anno_1, df$anno_2)

  out <- data.frame()

  for (i in rownames(tab)) {
    for (j in colnames(tab)) {

      intersection <- tab[i,j]

      union <- sum(tab[i,]) +
        sum(tab[,j]) -
        intersection

      jaccard <- ifelse(
        union == 0,
        NA,
        intersection / union
      )

      out <- rbind(
        out,
        data.frame(
          anno_1 = i,
          anno_2 = j,
          jaccard = jaccard
        )
      )
    }
  }

  out
}

jaccard_df <- lapply(
  unique(anno_df$qcFlag),
  function(qc) {

    sub <- anno_df[
      anno_df$qcFlag == qc,
    ]

    res <- calc_jaccard(sub)

    res$qcFlag <- qc

    res
  }
)

jaccard_df <- bind_rows(jaccard_df)

# Plot QC impact
ggplot(
  jaccard_df,
  aes(
    x = interaction(anno_1, anno_2),
    y = jaccard,
    fill = qcFlag
  )
) +
  geom_bar(
    stat = "identity",
    position = "dodge"
  ) +
  coord_flip() +
  theme_classic() +
  labs(
    x = "Annotation overlap",
    y = "Jaccard index",
    title = "Effect of QC on annotation agreement"
  )



### Keeping the independent cell type methods 
SeuObj$lyn_minor_celltype__cell_based_any__cluster_gated <- NULL
SeuObj$lyn_minor_celltype__cell_based_cluster_gated <- NULL
SeuObj$lyn_minor_celltype__cell_based_all <- NULL
SeuObj$major_celltype # ok 
SeuObj$lyn_minor_celltype__cell_based_any # ok 

cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- SeuObj$major_celltype %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)


DimPlot(SeuObj,
        reduction = "UMAP_regular",
        group.by = "major_celltype", 
        cols = celltype_cols)


cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- SeuObj$lyn_minor_celltype__cell_based_any %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)


DimPlot(subset(SeuObj, lyn_minor_celltype__cell_based_any %in% "other", invert = T),
        reduction = "UMAP_regular",
        group.by = "lyn_minor_celltype__cell_based_any", 
        cols = celltype_cols)



#qsave(SeuObj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_all467158cells_clustered.qs")



### Bcell, Tcd4, Tcd8 -------
# Just running it again to keep this code together with the cDC/pDC
SeuObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_all467158cells_clustered.qs")
library(Seurat)
library(Matrix)

## lyn_minor_celltype (Tcd4, Tcd8, Bcell)
# Marker definitions
marker_list <- list(
  Ptprc = c("Ptprc"),
  Cd3   = c("Cd3d", "Cd3e", "Cd3g"),
  Cd8   = c("Cd8a", "Cd8b1"),
  Cd4   = c("Cd4"),
  Cd19  = c("Cd19")
)

# Annotation rules
annotation_rules <- list(

  Tcd4 = expression(
    Ptprc & Cd3 & Cd4 & !Cd19 & !Cd8
  ),

  Tcd8 = expression(
    Ptprc & Cd3 & Cd8 & !Cd19 & !Cd4
  ),

  Bcell = expression(
    Ptprc & Cd19 & !Cd3
  )

)

# Run annotation
SeuObj <- cell_based_anno(
  seurat_obj = SeuObj,
  assay = "Nanostring",
  slot = "counts",
  mode = "any",  # choose: "any" or "all"

  anno_column = c("lyn_minor_celltype__cell_based_any"),

  marker_list = marker_list,

  annotation_rules = annotation_rules,

  # # Optional custom thresholds
  # thresholds = list(
  #   Ptprc = 2,
  #   Cd4   = 1,
  #   Cd19  = 1 ....
  # ),

  default_threshold = 1
)


SeuObj <- cell_based_anno(
  seurat_obj = SeuObj,
  assay = "Nanostring",
  slot = "counts",
  mode = "all",  # choose: "any" or "all"

  anno_column = c("lyn_minor_celltype__cell_based_all"),

  marker_list = marker_list,

  annotation_rules = annotation_rules,

  # # Optional custom thresholds
  # thresholds = list(
  #   Ptprc = 2,
  #   Cd4   = 1,
  #   Cd19  = 1 ....
  # ),

  default_threshold = 1
)

table(SeuObj$lyn_minor_celltype__cell_based_all, useNA = "always")


### pDC and cDC -----------
SeuObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_all467158cells_clustered.qs")

counts <- SeuObj@assays$Nanostring$counts
gene_panel <- rownames(counts)

# pDC
gene_panel[grepl("Ptprc", gene_panel)]  #Ptprc
gene_panel[grepl("Itgax", gene_panel)] #Itgax
gene_panel[grepl("Bst2", gene_panel)] #Bst2

# cDC
gene_panel[grepl("Ptprc", gene_panel)]  #Ptprc
gene_panel[grepl("Itgax", gene_panel)] #Itgax
gene_panel[grepl("Cd8", gene_panel)] #Cd8a Cd8b1
gene_panel[grepl("Itgam", gene_panel)] #Itgam



library(Seurat)
library(Matrix)

## dend_minor_celltype (pDC, cDC)

marker_list <- list(
  Ptprc = c("Ptprc"),
  Itgax = c("Itgax"),
  Bst2  = c("Bst2"),
  Cd8   = c("Cd8a", "Cd8b1") #, Itgam = c("Itgam")
 
)

annotation_rules <- list(

  pDC = expression(
    Ptprc & Itgax & Bst2
  ),

  cDC = expression(
    Ptprc & Itgax & Cd8 #& Itgam
  )

)

SeuObj <- cell_based_anno(
  seurat_obj = SeuObj,
  assay = "Nanostring",
  slot = "counts",
  mode = "any",  # choose: "any" or "all"

  anno_column = "dend_minor_celltype__cell_based_any",

  marker_list = marker_list,

  annotation_rules = annotation_rules,

  default_threshold = 1
)

table(SeuObj$dend_minor_celltype__cell_based_any, useNA = "always")


cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- SeuObj$dend_minor_celltype__cell_based_any %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)


DimPlot(subset(SeuObj, dend_minor_celltype__cell_based_any %in% "other", invert = T),
        reduction = "UMAP_regular",
        group.by = "dend_minor_celltype__cell_based_any", 
        cols = celltype_cols)


table(SeuObj$dend_minor_celltype__cell_based_any,
      SeuObj$lyn_minor_celltype__cell_based_any)


### Epithelial -----------
counts <- SeuObj@assays$Nanostring$counts
gene_panel <- rownames(counts)

# Club_Cells (CC10/SCGB1A1)
gene_panel[grepl("Epcam", gene_panel)]  
gene_panel[grepl("Scgb1a1", gene_panel)] # not exist 
# AT1 
gene_panel[grepl("Epcam", gene_panel)]  
gene_panel[grepl("Pdpn", gene_panel)] # not exist 
gene_panel[grepl("Cav1", gene_panel)] 
# AT2
gene_panel[grepl("Epcam", gene_panel)]  
gene_panel[grepl("Sftpc", gene_panel)] # not exist
# Goblet_cells
gene_panel[grepl("Epcam", gene_panel)]  
gene_panel[grepl("Muc5ac", gene_panel)] # not exist
# Ciliated_cells
gene_panel[grepl("Epcam", gene_panel)]  
gene_panel[grepl("Tuba1a", gene_panel)] # not exist
gene_panel[grepl("Foxj1", gene_panel)] # not exist
# Basal_cells
gene_panel[grepl("Epcam", gene_panel)]  
gene_panel[grepl("Krt5", gene_panel)] # not exist
# Neuroendocrine_cells
gene_panel[grepl("Epcam", gene_panel)]  
gene_panel[grepl("Calca", gene_panel)] # not exist
gene_panel[grepl("Chga", gene_panel)] # not exist


# that's why we cant subtype epthilal cells only with their markers


### Major-guided (not gated !) final cell-based annotation --------
minor_to_major <- list(
  # Tcd4  = c("Lymphocytes", "Dendritic"),
  # Tcd8  = c("Lymphocytes", "Dendritic"),
  # Bcell = c("Lymphocytes", "Dendritic"),
  Tcd4  = c("Lymphocytes"),
  Tcd8  = c("Lymphocytes"),
  Bcell = c("Lymphocytes"),
  pDC   = c("Dendritic"),
  cDC   = c("Dendritic")
)

SeuObj <- resolve_cell_based_conflicts(
  seurat_obj = SeuObj,
  minor_cols = c(
    "lyn_minor_celltype__cell_based_any",
    "dend_minor_celltype__cell_based_any"
  ),
  other_label = c("other"),
  major_col = "major_celltype",
  minor_to_major = minor_to_major,
  minor_preference = c("pDC", "cDC", "Tcd4", "Tcd8", "Bcell"),
  output_col = "cell_based_final"
)

table(SeuObj$cell_based_final, useNA = "always")
table(SeuObj$lyn_minor_celltype__cell_based_any, useNA = "always")
table(SeuObj$dend_minor_celltype__cell_based_any, useNA = "always")

tab <- table(SeuObj$cell_based_final, SeuObj$major_celltype)
# Row-normalized proportions
spec <- prop.table(tab, margin = 1)
pheatmap(spec)


cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- SeuObj$cell_based_final %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)


DimPlot(subset(SeuObj, cell_based_final %in% c("other", "Ambiguous_overlap_across_diff_cell_based"), invert = T),
        reduction = "UMAP_regular",
        group.by = "cell_based_final", 
        cols = celltype_cols)


table(SeuObj$cell_based_final)
table(SeuObj$dend_minor_celltype__cell_based_any)
table(SeuObj$lyn_minor_celltype__cell_based_any)


# Sum non-"other" annotations for multiple columns
cols_use <- c(
  "lyn_minor_celltype__cell_based_any",
  "dend_minor_celltype__cell_based_any",
  "cell_based_final"
)

sapply(cols_use, function(col) {

  tab <- table(SeuObj[[col]])

  sum(
    tab[names(tab) != "other"]
  )

})





### Check number of cells across groups ----
table(SeuObj@meta.data$Group, SeuObj@meta.data$cell_based_final)
library(ggplot2)

df <- as.data.frame(
  table(
    SeuObj@meta.data$Group,
    SeuObj@meta.data$cell_based_final
  )
)


colnames(df) <- c("Group", "Celltype", "Count")
df <- df[!df$Celltype %in% "other", ]
ggplot(df, aes(x = Group, y = Count, fill = Celltype)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(
    x = "Group",
    y = "Cell Count"
  )


# comparisons 
# primary infection 
table(SeuObj@meta.data$Group)

prim_inf_comp <- list (
  c("Neonate RSV infected (NO IFN)", 
  "Neonate IFN and RSV infected"),
  
  c("Neonate RSV infected (NO IFN)",
     "Adult RSV infected"),
  
  c("Neonate IFN and RSV infected",
    "Adult RSV infected")
)

table(SeuObj@meta.data$Group)

reinf_comp <- list (
  c("Neonate (NO IFN) reinfected",
    "Neonate IFN and RSV reinfected"),
  
  c("Neonate (NO IFN) reinfected",
    "Adult reinfected"),
  
  c("Neonate IFN and RSV reinfected",
    "Adult reinfected")
)



library(dplyr)
library(ggplot2)
library(patchwork)

meta <- SeuObj@meta.data

# comparisons
prim_inf_comp <- list(
  c("Neonate RSV infected (NO IFN)", "Neonate IFN and RSV infected"),
  c("Neonate RSV infected (NO IFN)", "Adult RSV infected"),
  c("Neonate IFN and RSV infected", "Adult RSV infected")
)

reinf_comp <- list(
  c("Neonate (NO IFN) reinfected", "Neonate IFN and RSV reinfected"),
  c("Neonate (NO IFN) reinfected", "Adult reinfected"),
  c("Neonate IFN and RSV reinfected", "Adult reinfected")
)

all_comp <- c(prim_inf_comp, reinf_comp)

names(all_comp) <- c(
  "Primary: Neo RSV no IFN vs Neo IFN+RSV",
  "Primary: Neo RSV no IFN vs Adult RSV",
  "Primary: Neo IFN+RSV vs Adult RSV",
  "Reinfected: Neo no IFN vs Neo IFN+RSV",
  "Reinfected: Neo no IFN vs Adult",
  "Reinfected: Neo IFN+RSV vs Adult"
)

# build comparison dataframe
plot_df <- lapply(names(all_comp), function(comp_name) {
  
  groups_use <- all_comp[[comp_name]]
  
  meta %>%
    dplyr::filter(
      Group %in% groups_use,
      !cell_based_final %in% c(
        "other",
        "Ambiguous_overlap_across_diff_cell_based"
      )
    ) %>%
    dplyr::count(Group, cell_based_final, name = "n_cells") %>%
    dplyr::mutate(comparison = comp_name)
  
}) %>%
  dplyr::bind_rows()

# barplot
ggplot(plot_df, aes(x = cell_based_final, y = n_cells, fill = Group)) +
  geom_col(position = "dodge") +
  facet_wrap(~ comparison, scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 10),
    legend.title = element_blank()
  ) +
  labs(
    x = "cell_based_final",
    y = "Number of cells"
  )


#qsave(SeuObj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_all467158cells_clustered.qs")



### Correcting cell annotation (major cell type) - all the other major cell types are pretty good ---

# Myeloid  ----
mk_general <- get_top_markers(
  obj = SeuObj,
  cluster_col = "cluster_regular_res06",
  assay = assay,
  top_n = 5
)
table(SeuObj$major_celltype, 
      SeuObj$cluster_regular_res06) #4, 7, 0, 11, 12 -> figure out which one is actually myeloid 

mk_general$top_markers[mk_general$top_markers$group %in% c("4", "7", "0", "11", "12"), ] %>% data.frame

# 7 is our myeloid 
# Myeloid / macrophage-like APCs (Cd74, Tyrobp, Apoe, Ctsd, Psap)
# 
SeuObj@meta.data$major_celltype <- as.character(
  SeuObj@meta.data$major_celltype
)
SeuObj@meta.data[
  SeuObj@meta.data$major_celltype %in% "Myeloid",
  "major_celltype"
] <- "Stroma"

SeuObj@meta.data[
  SeuObj@meta.data$cluster_regular_res06 %in% "7",
  "major_celltype"
] <- "Myeloid"



# Dendritic ----
table(SeuObj$major_celltype, 
      SeuObj$cluster_regular_res06) #1, 14, 5, 8 

mk_general$top_markers[mk_general$top_markers$group %in% c("1", "14", "5", "8"), ] %>% data.frame()

SeuObj@meta.data$major_celltype <- as.character(
  SeuObj@meta.data$major_celltype
)
SeuObj@meta.data[
  SeuObj@meta.data$major_celltype %in% "Dendritic",
  "major_celltype"
] <- "other"

SeuObj@meta.data[
  SeuObj@meta.data$cluster_regular_res06 %in% "1",
  "major_celltype"
] <- "Dendritic"
SeuObj@meta.data[
  SeuObj@meta.data$cluster_regular_res06 %in% "8",
  "major_celltype"
] <- "Lymphocytes"
SeuObj@meta.data[
  SeuObj@meta.data$cluster_regular_res06 %in% "5",
  "major_celltype"
] <- "Epithelial"


DimPlot(SeuObj, 
        group.by = "major_celltype", 
        reduction = "UMAP_regular") & NoAxes() 


# Lymphocytes ----
table(SeuObj$major_celltype, 
      SeuObj$cluster_regular_res06) #1, 10,  2, 3, 4. 

mk_general$top_markers[mk_general$top_markers$group %in% c("1", "10", "2", "3", "4"), ] %>% data.frame()

SeuObj@meta.data$major_celltype <- as.character(
  SeuObj@meta.data$major_celltype
)
SeuObj@meta.data[
  SeuObj@meta.data$major_celltype %in% "Lymphocytes",
  "major_celltype"
] <- "other"

SeuObj@meta.data[
  SeuObj@meta.data$cluster_regular_res06 %in% "10",
  "major_celltype"
] <- "Lymphocytes"

SeuObj@meta.data[
  SeuObj@meta.data$cluster_regular_res06 %in% c("2", "3"),
  "major_celltype"
] <- "Lymphocytes"


SeuObj@meta.data[
  SeuObj@meta.data$cluster_regular_res06 %in% c("4"),
  "major_celltype"
] <- "Stroma"


# Epithelial ----
table(SeuObj$major_celltype, 
      SeuObj$cluster_regular_res06) #15, 5, 6, 9

mk_general$top_markers[mk_general$top_markers$group %in% c("15","5","6","9"), ] %>% data.frame()

SeuObj@meta.data$major_celltype <- as.character(
  SeuObj@meta.data$major_celltype
)
SeuObj@meta.data[
  SeuObj@meta.data$major_celltype %in% "Epithelial",
  "major_celltype"
] <- "other"

SeuObj@meta.data[
  SeuObj@meta.data$cluster_regular_res06 %in% c("5", "6", "9"),
  "major_celltype"
] <- "Epithelial"

# Endothelial ----
table(SeuObj$major_celltype, 
      SeuObj$cluster_regular_res06) #0, 11, 12,

mk_general$top_markers[mk_general$top_markers$group %in% c("0","11","12"), ] %>% data.frame()

SeuObj@meta.data$major_celltype <- as.character(
  SeuObj@meta.data$major_celltype
)
SeuObj@meta.data[
  SeuObj@meta.data$major_celltype %in% "Endothelial",
  "major_celltype"
] <- "other"

SeuObj@meta.data[
  SeuObj@meta.data$cluster_regular_res06 %in% c("0"),
  "major_celltype"
] <- "Endothelial"

table(SeuObj$major_celltype)

DimPlot(SeuObj, 
        group.by = "major_celltype", 
        reduction = "UMAP_regular") & NoAxes() 

#qsave(SeuObj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_all467158cells_clustered.qs")


### Correcting cell annotation (major/minor cell type) - all the other major cell types are pretty good ---
# Dendritic cells
assay <- "Nanostring"
scale_factor <- 1000

obj <- run_tma_clustering(
  obj = subset(SeuObj, major_celltype %in% c("Dendritic")),
  assay = assay,
  dims_use = 1:15,
  resolutions = c(0.3, 0.6),
  scale_factor = scale_factor
)

DimPlot(obj,
        group.by = "cluster_regular_res06",
        reduction = "UMAP_regular") & NoAxes()



# just checking annotation - looks good! no action needed here
mk_general <- get_top_markers(
  obj = obj,
  cluster_col = "cluster_regular_res06",
  assay = assay,
  top_n = 5
)


# from chatGPT
cDC_genes <- c(
  "Flt3",

  "Itgax",

  "Cd74",

  "H2-Ab1",

  "H2-Eb1",

  "Ciita",

  "Cd83",

  "Cd86",

  "Cd80"
)

pDC_genes <- c(
  "Bst2",
  "Il3ra",
  "Ly6d",
  "Tlr7",
  "Tlr8",
  "Ifi44l",
  "Ifit1",
  "Ifit3/b",
  "Oas1a/g",
  "Oas2",
  "Oas3",
  "Irf3",
  "Gzmb",
  "Ccr9",
  "Jchain"
)

dendritic_diego <- c("Cd74","H2-Ab1", "B2m", "Ppia", "Ciita", "Itgb2")
plot_marker_dotplot(
  obj,
  cluster_col = "cluster_regular_res06",
  #features = unique(mk_general$top_markers[, ]$feature),
  features = dendritic_diego, #pDC_genes, #cDC_genes,
  dot_shape = 21,
  dot_border_color = "black",
  dot_border_width = 0.5,
  color_palette = "magma",
  color_direction = 1,
  min_dot_size = 0,
  max_dot_size = 10,
  x_text_angle = 90,
  x_text_hjust = 1,
  x_text_vjust = 1,
  y_text_size = 22,
  x_text_size = 16,
  legend_title_size = 11,
  legend_text_size = 10,
  plot_title = "",
  base_size = 12,
  theme_base = ggplot2::theme_classic()
)

mk_general$top_markers[, ] %>% data.frame()
mk_general$markers[mk_general$markers$feature %in% cDC_genes, ] %>% data.frame()
mk_general$markers[mk_general$markers$feature %in% pDC_genes, ] %>% data.frame()


table(obj$cluster_regular_res06)

obj@meta.data$major_celltype
obj@meta.data[obj@meta.data$cluster_regular_res06 %in% "3",  ]$major_celltype <- "Endothelial"
obj@meta.data[obj@meta.data$cluster_regular_res06 %in% "1",  ]$major_celltype <- "Stroma"
obj@meta.data[obj@meta.data$cluster_regular_res06 %in% "0",  ]$major_celltype <- "other"

SeuObj$major_celltype[rownames(obj@meta.data)] <-
  obj$major_celltype

DimPlot(subset(SeuObj, major_celltype %in% "other", invert = T), 
        group.by = "major_celltype", 
        reduction = "UMAP_regular",
        cols = celltype_cols) & NoAxes()


assay <- "Nanostring"
scale_factor <- 1000

obj <- run_tma_clustering(
  obj = subset(SeuObj, major_celltype %in% c("Dendritic")),
  assay = assay,
  dims_use = 1:15,
  resolutions = c(0.3, 0.6),
  scale_factor = scale_factor
)

DimPlot(obj,
        group.by = "cluster_regular_res03",
        reduction = "UMAP_regular") & NoAxes()

 
mk_general <- get_top_markers(
  obj = obj,
  cluster_col = "cluster_regular_res03",
  assay = assay,
  top_n = 20
)

plot_marker_dotplot(
  obj,
  cluster_col = "cluster_regular_res03",
  features = unique(mk_general$top_markers[, ]$feature),
  #features = dendritic_diego, #pDC_genes, #cDC_genes,
  dot_shape = 21,
  dot_border_color = "black",
  dot_border_width = 0.5,
  color_palette = "magma",
  color_direction = 1,
  min_dot_size = 0,
  max_dot_size = 10,
  x_text_angle = 90,
  x_text_hjust = 1,
  x_text_vjust = 1,
  y_text_size = 22,
  x_text_size = 16,
  legend_title_size = 11,
  legend_text_size = 10,
  plot_title = "",
  base_size = 12,
  theme_base = ggplot2::theme_classic()
)


mk_general$top_markers[, ] %>% data.frame()
#0 Lamp3+ activated dendritic cells (mregDC-like) #cDCs
#1 Apoe+/C1q+ macrophage/APC Medium-high
# no pDCs ...

obj@meta.data$minor_celltype <- NA
obj@meta.data[obj@meta.data$cluster_regular_res03 %in% "0",  ]$minor_celltype <- "cDC (Dendritic)"
obj@meta.data[obj@meta.data$cluster_regular_res03 %in% "1",  ]$minor_celltype <- "Apoe+/C1q+ macrophage/APC (Dendritic)"

SeuObj@meta.data$minor_celltype <- NA
SeuObj$minor_celltype[rownames(obj@meta.data)] <-
  obj$minor_celltype

table(SeuObj$minor_celltype)

# qsave(SeuObj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_all467158cells_clustered.qs")