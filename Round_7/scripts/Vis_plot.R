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

# Merged obj 
SeuObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/data/Processed/merged_Sobj_UMAP.qs")

# Seurat obj list (spatial layers)
obj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/Sobj_list_TMA_ID_fixed.qs")
obj_list$`NA` <- NULL
obj_list$`NA.1` <- NULL

# SeuObj_BTcell <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/lymphocytes/BcellTcell_Sobj.qs")
# 
# SeuObj_DCcell <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/dendritic/pDCcDC_Sobj.qs")

### Settings --------
cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- SeuObj$major_celltype %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)

celltype_BTcell_cols <- c("Bcell" = "#F0E442", 
                   "Tcd4" = "#E69F00", 
                   "Tcd8" = "#D55E00")


celltype_levels <- c(
  # Epithelial lineage
  "Alveolar_bipotent_progenitor",
  "AT2_cell",
  "AT1_cell",
  "Epithelial_cell",
  
  # Myeloid / innate immune
  "Monocyte_progenitor",
  "Interstitial_macrophage",
  "Alveolar_macrophage",
  "Dendritic_cell",
  "Granulocyte",
  "Nuocyte",
  
  # Lymphoid / adaptive immune
  "Lymphocyte",
  
  # Structural / vascular
  "Endothelial_cell",
  "Stromal_cell"
)

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


### Figures --------------------------------------------------


### Figure 1A ----
# Diego: No edits 
DimPlot(SeuObj, 
        group.by = "major_celltype", 
        cols = celltype_cols) & NoAxes()

### Figure 1B ----
# Diego: For CD4 and CD8, we would like clarification to ensure that these signals represent true T cell subsets (i.e., CD3+CD4+ helper T cells and CD3+CD8+ cytotoxic T cells), rather than any cell expressing CD4 or CD8 individually. Please confirm whether these markers are restricted to annotated T cell populations or based on single-marker expression alone. Our goal is to ensure that this panel reflects lineage-defined cell types accurately.

# Maycon: I won't change this figure cause it's made to show generic markers. I could do it on figure 2 and 3 where we're talking about specific cell types (Tcd4, pDCs ...)


### Figure 1C ----
# From Diego - matching into our current major_celltypee markers (from May 12, 2026)
geneSets <- list(
  Dendritic = c("Cd74", "H2-Ab1", "B2m", "Ppia", "Ciita", "Itgb2"), #not too great
  
  Epithelial = c("Tubb4b", "Sod1", "Cd24a", "Hsp90aa1", "Krt8", "Cd9"), #okay
  
  Bcell = c("Jchain", "Igkc", "Ighm", "Ly6d", "Cd37", "Blk"), # okay
  
  Tcell = c("Cd274", "Cd4", "Gata3", "Tnfrsf9", "Ly75", "Cd59a"), # okay
  
  Myeloid = c("Elane", "Ctsg", "Srgn", "Cd1d1", "Cd48", "Clec10a"), #okay
  
  Endothelial = c("Cdh5", "Tie1", "Esam", "Pecam1", "Cd93", "Icam2") #okay
)

# Keep only genes present in your dataset
geneSets <- lapply(geneSets, function(gs) intersect(gs, rownames(obj)))

obj <- SeuObj
obj@assays$Nanostring$data <- NULL
obj@assays$Nanostring$scale.data <- NULL
assay <- "Nanostring"
scale_factor <- 1000
obj <- NormalizeData(obj, assay = assay, scale.factor = scale_factor)
# mk_general <- get_top_markers(
#   obj = obj,
#   cluster_col = "celltype_final",
#   assay = assay,
#   top_n = 5
# )

plot_marker_dotplot(
  obj,
  cluster_col = "major_celltype",
  features = unique(unlist(geneSets)),
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



### Figure 2A ----
# Diego: No edits 
DimPlot(SeuObj_BTcell, 
        group.by = "lyn_minor_celltype", 
        cols = celltype_BTcell_cols) & NoAxes()

### Figure 2B ----
gene_panel <- SeuObj_BTcell@assays$Nanostring$counts %>% rownames()
gene_panel[grepl("Cd3", gene_panel)]

## Tcell 
genes <- c("Cd4", "Cd3d")
celltype_col <- "lyn_minor_celltype"
celltype_val <- "Tcd4"

counts <- SeuObj_BTcell@assays$Nanostring$counts[genes, ]
expr <- t(as.matrix(counts))

SeuObj_BTcell$gene1_gene2_pos <- expr[, genes[1]] > 0 & expr[, genes[2]] > 0

SeuObj_sub <- SeuObj_BTcell[
  ,
  SeuObj_BTcell$gene1_gene2_pos &
    SeuObj_BTcell@meta.data[[celltype_col]] %in% celltype_val
]

# blend expression of the two genes in the subset
p1 <- FeaturePlot(
  SeuObj_sub,
  features = genes,
  blend = TRUE,
  order = TRUE,
  cols = c("grey90", "blue", "red")
) & NoAxes()



## Tcell 
genes <- c("Cd4", "Cd3e")
celltype_col <- "lyn_minor_celltype"
celltype_val <- "Tcd4"

counts <- SeuObj_BTcell@assays$Nanostring$counts[genes, ]
expr <- t(as.matrix(counts))

SeuObj_BTcell$gene1_gene2_pos <- expr[, genes[1]] > 0 & expr[, genes[2]] > 0

SeuObj_sub <- SeuObj_BTcell[
  ,
  SeuObj_BTcell$gene1_gene2_pos &
    SeuObj_BTcell@meta.data[[celltype_col]] %in% celltype_val
]

# blend expression of the two genes in the subset
p2 <- FeaturePlot(
  SeuObj_sub,
  features = genes,
  blend = TRUE,
  order = TRUE,
  cols = c("grey90", "blue", "red")
) & NoAxes()



## Tcell 
genes <- c("Cd4", "Cd3g")
celltype_col <- "lyn_minor_celltype"
celltype_val <- "Tcd4"

counts <- SeuObj_BTcell@assays$Nanostring$counts[genes, ]
expr <- t(as.matrix(counts))

SeuObj_BTcell$gene1_gene2_pos <- expr[, genes[1]] > 0 & expr[, genes[2]] > 0

SeuObj_sub <- SeuObj_BTcell[
  ,
  SeuObj_BTcell$gene1_gene2_pos &
    SeuObj_BTcell@meta.data[[celltype_col]] %in% celltype_val
]

# blend expression of the two genes in the subset
p3 <- FeaturePlot(
  SeuObj_sub,
  features = genes,
  blend = TRUE,
  order = TRUE,
  cols = c("grey90", "blue", "red")
) & NoAxes()


#p1 | p2 | p3

p1 / p2 / p3


library(UpSetR)
gene1 <- "Cd4"
cd3_genes <- c("Cd3d", "Cd3e", "Cd3g")
all_genes <- c(gene1, cd3_genes)

celltype_col <- "major_celltype"
celltype_val <- "Lymphocyte"
obj <- SeuObj

counts <- obj@assays$Nanostring$counts[all_genes, ]
expr <- t(as.matrix(counts))

cells_keep <- obj@meta.data[[celltype_col]] %in% celltype_val
expr_sub <- expr[cells_keep, ]

# Cell IDs called by each Cd4 + Cd3x rule
calls <- lapply(cd3_genes, function(g) {
  rownames(expr_sub)[expr_sub[, gene1] > 0 & expr_sub[, g] > 0]
})
names(calls) <- cd3_genes

# Total unique cells called by ANY Cd3 gene
cells_any <- Reduce(union, calls)

# Cells called by ALL Cd3 genes
cells_all <- Reduce(intersect, calls)

# Summary per strategy
summary_df <- data.frame(
  strategy = c(cd3_genes, "ANY_Cd3", "ALL_Cd3"),
  n_cells_called = c(
    sapply(calls, length),
    length(cells_any),
    length(cells_all)
  )
)

summary_df

gene_panel[grepl("Ptprc", gene_panel)]
gene_panel[grepl("Cd4", gene_panel)]
gene_panel[grepl("Cd3", gene_panel)]
gene_panel[grepl("Cd8", gene_panel)]
gene_panel[grepl("Cd19", gene_panel)]

gene_panel[grepl("Itgax", gene_panel)]
gene_panel[grepl("Bst2", gene_panel)]
gene_panel[grepl("Itgam", gene_panel)]


# writeLines(gene_panel, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/gene_panel_1k.txt")

library(presto)
group_col <- "major_celltype"
SeuObj@assays$RNA <- SeuObj@assays$Nanostring
markers <- wilcoxauc(SeuObj, group_col)
markers %>% head()
dim(markers)
top_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  dplyr::top_n(wt = logFC, n = 30) %>%
  ungroup() 
top_markers <- data.frame(top_markers)
top_markers <- top_markers[, c("feature", "group",  "logFC", "padj")]
names(top_markers) <- c("gene", "general_celltype", "logFC", "padj")
head(top_markers)
# write.csv(top_markers, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/top30_markers_majorcelltype.csv")



### I notice there some "gene name variants" for cd3 ... I need to ask them to annotate the cell types based on our panel!! 
### I'm asking Diego to go over the genes we actually have in the panel + the current cell type markers we have. 