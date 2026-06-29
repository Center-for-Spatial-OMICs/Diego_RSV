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

### Loading variables -----
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



### Load data ------
SeuObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/data/Processed/merged_Sobj_UMAP.qs")
SeuObj@assays$niche <- NULL

obj <- SeuObj
assay <- "Nanostring"
scale_factor <- 1000
obj <- run_tma_clustering(
  obj = obj,
  assay = assay,
  dims_use = 1:30,
  resolutions = c(0.3, 0.5),
  scale_factor = scale_factor
)


cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- obj$major_celltype %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)
DimPlot(obj, 
        reduction = "UMAP_regular", 
        group.by = "major_celltype", 
        cols = celltype_cols) 

cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- obj$cluster_regular_res03 %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)


DimPlot(obj, 
        reduction = "UMAP_regular", 
        group.by = "cluster_regular_res03", 
        cols = celltype_cols) 




mk_general <- get_top_markers(
  obj = obj,
  cluster_col = "cluster_regular_res03",
  assay = assay,
  top_n = 10
)

# Diego's cell type markers
plot_marker_dotplot(
  obj,
  cluster_col = "cluster_regular_res03",
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

# top cluster markers
plot_marker_dotplot(
  obj,
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
  file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/cell_subtype/top_cluster_markers.csv",
  row.names = FALSE
)


cluster_labels <- c(

  "0"  = "Cd274_Cd68_immune",

  "1"  = "Dendritic",

  "2"  = "Endothelial",

  "3"  = "Fibroblast",

  "4"  = "Lamp3_Cd74_dendritic_like",

  "5"  = "Vegfa_Icam1_epithelial_stress",

  "6"  = "Epithelial",

  "7"  = "Bcell",

  "8"  = "Macrophage",

  "9"  = "Myeloid",

  "10" = "Cardiomyocyte",

  "11" = "Pericyte"

)

Idents(obj) <- "cluster_regular_res03"
obj@meta.data$celltype_marker_based <- cluster_labels[as.character(Idents(obj))]

cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- obj$celltype_marker_based %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)

table(obj$celltype_marker_based, useNA = "always")
table(obj$celltype_marker_based, useNA = "always")
DimPlot(obj, 
        reduction = "UMAP_regular", 
        group.by = "celltype_marker_based", 
        cols = celltype_cols) 


# Diego's cell type markers
plot_marker_dotplot(
  obj,
  cluster_col = "celltype_marker_based",
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

# top cluster markers
mk_general <- get_top_markers(
  obj = obj,
  cluster_col = "celltype_marker_based",
  assay = assay,
  top_n = 5
)
plot_marker_dotplot(
  obj,
  cluster_col = "celltype_marker_based",
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


# top cluster markers %in% chatBot anchored genes to cell type it
cluster_markers_anchors <- list(
  Cd274_Cd68_immune = c("Adm2", "Parp1", "Cd274", "Ror1", "Hdac11"),
  
  Dendritic = c("Cd74", "B2m", "Ccl5", "H2-Eb1", "H2-D1"),
  
  Endothelial = c("Cdh5", "Ramp2", "Eng", "Ace", "Col4a1"),
  
  Fibroblast = c("Gsn", "Col1a2", "Bgn", "Col1a1", "Itga8"),
  
  Lamp3_Cd74_dendritic_like = c("Cd74", "Lamp3", "Chil1", "Etv5", "Fgfr2"),
  
  Vegfa_Icam1_epithelial_stress = c("Vegfa", "Spock2", "Icam1", "S100a6", "Igfbp7"),
  
  Epithelial = c("Sod1", "Tubb4b", "Ppia", "Tpt1", "Cd24a"),
  
  Bcell = c("Igkc", "Ighm", "Cd74", "Cd37", "H2-DMb2"),
  
  Macrophage = c("Psap", "Ctsd", "Tyrobp", "Gpx1", "Cotl1"),
  
  Myeloid = c("S100a9", "S100a8", "Csf3r", "Tyrobp", "Srgn"),
  
  Cardiomyocyte = c("Myh6", "Myl7", "Tpm1", "Tcap", "Myl4"),
  
  Pericyte = c("Efnb2", "Cxcl12", "Pdgfrb", "Notch3", "Adgrg6")
)

plot_marker_dotplot(
  obj,
  cluster_col = "celltype_marker_based",
  features = unique(unlist(cluster_markers_anchors)),
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


# Diego's cell type markers on our ref based cell type
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

# top cluster markers 
mk_general <- get_top_markers(
  obj = obj,
  cluster_col = "major_celltype",
  assay = assay,
  top_n = 5
)
plot_marker_dotplot(
  obj,
  cluster_col = "major_celltype",
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