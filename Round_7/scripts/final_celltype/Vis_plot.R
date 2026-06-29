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

source("/mnt/scratch2/Maycon/Utils/R_codabase/Utils.R")

### Loading data ------
# Merged obj 
SeuObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_all467158cells_clustered.qs")

# Seurat obj list (spatial layers)
obj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/Sobj_list_TMA_ID_fixed.qs")
obj_list$`NA` <- NULL
obj_list$`NA.1` <- NULL

### Settings --------
# cols <-
#   c(
#     '#8DD3C7', '#BEBADA', '#FB8072', '#FDB462', '#80B1D3', '#B3DE69', '#FCCDE5',
#     '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
#     '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
#   )
# factor_levels <- SeuObj$major_celltype %>% table() %>% names()
# celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)

celltype_cols <- c(
  "Dendritic"    = "#8DD3C7",
  "Endothelial"  = "#E6C7A8",
  "Epithelial"   = "#C28E5B",
  "Lymphocytes"  = "gold",
  "Stroma"       = "#8C6D5A",
  "Myeloid"      = "#80B1D3"
)

celltype_BTcell_cols <- c("Bcell" = "#F0E442", 
                   "Tcd4" = "#E69F00", 
                   "Tcd8" = "#D55E00")

celltype_DCs_cols <- c("cDC (Dendritic)" = "#7f3f00",
                       "Apoe+/C1q+ macrophage/APC (Dendritic)" = "#ff7f00")

# Major cell types supervised-markers
marker_list <- list(
"Dendritic"=c("Cd74","H2-Ab1","B2m","Ppia","Ciita","Itgb2"),
"Epithelial"=c("Tubb4b","Sod1","Cd24a","Hsp90aa1","Krt8","Cd9"), #Epcam would work too
"B_Cells"=c("Jchain","Igkc","Ighm","Ly6d","Cd37","Blk"),
"T_Cells"=c("Cd274","Cd4","Gata3","Tnfrsf9","Ly75","Cd59a"),
"Myeloid"=c("Elane","Ctsg","Srgn","Cd1d1","Cd48","Clec10a"),
"Endothelial"=c("Cdh5","Tie1","Esam","Pecam1","Cd93","Icam2")
)

# Lymphocytes subtype supervised-markers 
marker_list <- list(
  Ptprc = c("Ptprc"),
  Cd3   = c("Cd3d", "Cd3e", "Cd3g"),
  Cd8   = c("Cd8a", "Cd8b1"),
  Cd4   = c("Cd4"),
  Cd19  = c("Cd19")
)

# Dendritic subtype supervised-markers
marker_list <- list(
  Ptprc = c("Ptprc"),
  Itgax = c("Itgax"),
  Bst2  = c("Bst2"),
  Cd8   = c("Cd8a", "Cd8b1") #, Itgam = c("Itgam")
)

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

### Figures --------------------------------------------------


### Figure 1A ----
# Diego: No edits 
DimPlot(subset(SeuObj, major_celltype %in% "other", invert = T), 
        group.by = "major_celltype", 
        reduction = "UMAP_regular",
        cols = celltype_cols) & NoAxes() 

# DimPlot(SeuObj, 
#         group.by = "cluster_regular_res06", 
#         reduction = "UMAP_regular") & NoAxes() 

### Figure 1B ----
DefaultAssay(SeuObj) <- "Nanostring"
Seurat::FeaturePlot(
  object = SeuObj,
  features = "Tubb4b",
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) +
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) +
NoAxes()


DefaultAssay(SeuObj) <- "Nanostring"
Seurat::FeaturePlot(
  object = SeuObj,
  features = "Pecam1",
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) +
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) +
NoAxes()



DefaultAssay(SeuObj) <- "Nanostring"
Seurat::FeaturePlot(
  object = SeuObj,
  features = c("Clec10a"),
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) &
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) &
NoAxes()










### Figure 1C ----
# From Diego - matching into our current major_celltypee markers (from May 12, 2026)

# Keep only genes present in your dataset
mk_general <- get_top_markers(
  obj = subset(SeuObj, major_celltype %in% "other", invert = T),
  cluster_col = "major_celltype",
  assay = assay,
  top_n = 5
)

Stroma_markers <- mk_general$top_markers[mk_general$top_markers$group %in% c("Stroma"), ]$feature

marker_list <- list(
"Dendritic"=c("Cd74","H2-Ab1","B2m","Ppia","Ciita","Itgb2"),
"Epithelial"=c("Tubb4b","Sod1","Cd24a","Hsp90aa1","Krt8","Cd9"), #Epcam would work too
"B_Cells"=c("Jchain","Igkc","Ighm","Ly6d","Cd37","Blk"),
"T_Cells"=c("Cd274","Cd4","Gata3","Tnfrsf9","Ly75","Cd59a"),
"Myeloid"= c("Elane","Ctsg","Srgn","Cd1d1","Cd48","Clec10a"),
"Endothelial"=c("Cdh5","Tie1","Esam","Pecam1","Cd93","Icam2")
)
geneSets <- lapply(marker_list, function(gs) intersect(gs, rownames(SeuObj)))

plot_marker_dotplot(
  subset(SeuObj, major_celltype %in% "other", invert = T),
  cluster_col = "major_celltype",
  features = c(unlist(geneSets) %>% as.vector(), Stroma_markers),
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




### Figure 1D -----------
## Make shorter group labels
library(dittoSeq)
library(ggplot2)
library(dplyr)

## ----------------------------------
## Shared plotting theme
## ----------------------------------

prop_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),

    axis.title.x = element_blank(),
    axis.title.y = element_blank(),

    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1
    ),

    axis.text.y = element_text(size = 14),

    axis.ticks = element_blank(),

    legend.title = element_blank(),
    legend.text = element_text(size = 12),

    strip.text.x = element_text(
      size = 16,
      face = "bold"
    )
  )

## ==================================
## Primary infection
## ==================================

primary_map <- c(
  "Neonate control" = "Neo Ctrl",
  "Neonate RSV infected (NO IFN)" = "Neo RSV",
  "Neonate IFN and RSV infected" = "Neo IFN+RSV",
  "Adult RSV infected" = "Adult RSV"
)

primary_order <- c(
  "Neo Ctrl",
  "Neo RSV",
  "Neo IFN+RSV",
  "Adult RSV"
)

SeuObj_primary <- subset(
  SeuObj,
  subset =
    major_celltype != "other" &
    Group %in% names(primary_map)
)

SeuObj_primary$Group_short <- recode(
  SeuObj_primary$Group,
  !!!primary_map
)

SeuObj_primary$Group_short <- factor(
  SeuObj_primary$Group_short,
  levels = primary_order
)

p_primary <- dittoBarPlot(
  SeuObj_primary,
  var = "major_celltype",
  group.by = "Group_short",
  color.panel = celltype_cols,
  main = "Primary infection"
) +
  scale_x_discrete(limits = primary_order) +
  prop_theme

p_primary


## ==================================
## Secondary infection
## ==================================

reinf_map <- c(
  "Neonate control" = "Neo Ctrl",
  "Neonate (NO IFN) reinfected" = "Neo Reinf",
  "Neonate IFN and RSV reinfected" = "Neo IFN Reinf",
  "Adult reinfected" = "Adult Reinf"
)

reinf_order <- c(
  "Neo Ctrl",
  "Neo Reinf",
  "Neo IFN Reinf",
  "Adult Reinf"
)

SeuObj_reinf <- subset(
  SeuObj,
  subset =
    major_celltype != "other" &
    Group %in% names(reinf_map)
)

SeuObj_reinf$Group_short <- recode(
  SeuObj_reinf$Group,
  !!!reinf_map
)

SeuObj_reinf$Group_short <- factor(
  SeuObj_reinf$Group_short,
  levels = reinf_order
)

p_reinf <- dittoBarPlot(
  SeuObj_reinf,
  var = "major_celltype",
  group.by = "Group_short",
  color.panel = celltype_cols,
  main = "Secondary infection"
) +
  scale_x_discrete(limits = reinf_order) +
  prop_theme

p_reinf


### Figure 1E -----------
# Just when we're sure about final fogure - it takes too much time to crop and plot they way people want it




### Figure 2A -----------

# From Bcell, Tcd4, Tcd8, recluster and plot its umap
tab <- table(SeuObj$cell_based_final, SeuObj$major_celltype)
# Row-normalized proportions
spec <- prop.table(tab, margin = 1)
pheatmap(spec)

assay <- "Nanostring"
scale_factor <- 1000

obj <- run_tma_clustering(
  obj = subset(SeuObj, cell_based_final %in% c("Bcell",
                                               "Tcd4",
                                               "Tcd8")),
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

plot_marker_dotplot(
  obj,
  cluster_col = "cluster_regular_res06",
  features = unique(mk_general$top_markers[, ]$feature),
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

table(obj$major_celltype,
      obj$cluster_regular_res06)
table(obj$cell_based_final,
      obj$cluster_regular_res06)

DimPlot(obj, 
        group.by = "cell_based_final", 
        reduction = "UMAP_regular", 
        cols = celltype_BTcell_cols) & NoAxes() 



#qsave(obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_Lyn4650cells_clustered.qs")

obj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_Lyn4650cells_clustered.qs")

DimPlot(obj, 
        group.by = "cell_based_final", 
        reduction = "UMAP_regular", 
        cols = celltype_BTcell_cols) & NoAxes() 

### Figure 2B -----------
DefaultAssay(obj) <- "Nanostring"
Seurat::FeaturePlot(
  object = obj,
  features = c("Cd19"),
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) &
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) &
NoAxes()


DefaultAssay(obj) <- "Nanostring"
Seurat::FeaturePlot(
  object = obj,
  features = c("Cd8a"), #and/or Cd8b1
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) &
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) &
NoAxes()


DefaultAssay(obj) <- "Nanostring"
Seurat::FeaturePlot(
  object = obj,
  features = c("Cd4"), #and/or Cd8b1
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) &
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) &
NoAxes()


### Figure 2C -----------
## ----------------------------------
## Shared plotting theme
## ----------------------------------

prop_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),

    axis.title.x = element_blank(),
    axis.title.y = element_blank(),

    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1
    ),

    axis.text.y = element_text(size = 14),

    axis.ticks = element_blank(),

    legend.title = element_blank(),
    legend.text = element_text(size = 12),

    strip.text.x = element_text(
      size = 16,
      face = "bold"
    )
  )

## ==================================
## Primary infection
## ==================================

primary_map <- c(
  "Neonate control" = "Neo Ctrl",
  "Neonate RSV infected (NO IFN)" = "Neo RSV",
  "Neonate IFN and RSV infected" = "Neo IFN+RSV",
  "Adult RSV infected" = "Adult RSV"
)

primary_order <- c(
  "Neo Ctrl",
  "Neo RSV",
  "Neo IFN+RSV",
  "Adult RSV"
)

obj_primary <- subset(
  obj,
  subset =
    cell_based_final != "other" &
    Group %in% names(primary_map)
)

obj_primary$Group_short <- recode(
  obj_primary$Group,
  !!!primary_map
)

obj_primary$Group_short <- factor(
  obj_primary$Group_short,
  levels = primary_order
)

p_primary <- dittoBarPlot(
  obj_primary,
  var = "cell_based_final",
  group.by = "Group_short",
  color.panel = celltype_BTcell_cols,
  main = "Primary infection"
) +
  scale_x_discrete(limits = primary_order) +
  prop_theme

p_primary


## ==================================
## Secondary infection
## ==================================

reinf_map <- c(
  "Neonate control" = "Neo Ctrl",
  "Neonate (NO IFN) reinfected" = "Neo Reinf",
  "Neonate IFN and RSV reinfected" = "Neo IFN Reinf",
  "Adult reinfected" = "Adult Reinf"
)

reinf_order <- c(
  "Neo Ctrl",
  "Neo Reinf",
  "Neo IFN Reinf",
  "Adult Reinf"
)

obj_reinf <- subset(
  obj,
  subset =
    cell_based_final != "other" &
    Group %in% names(reinf_map)
)

obj_reinf$Group_short <- recode(
  obj_reinf$Group,
  !!!reinf_map
)

obj_reinf$Group_short <- factor(
  obj_reinf$Group_short,
  levels = reinf_order
)

p_reinf <- dittoBarPlot(
  obj_reinf,
  var = "cell_based_final",
  group.by = "Group_short",
  color.panel = celltype_BTcell_cols,
  main = "Secondary infection"
) +
  scale_x_discrete(limits = reinf_order) +
  prop_theme

p_reinf


### Figure 2D -----------
# spatial plots - save it for later 

### Figure 2E -----------
# these volcano plots and barplots are all coming from here /mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/scripts/final_celltype/DEGs.R




### Figure 3A -----------

# From Bcell, Tcd4, Tcd8, recluster and plot its umap
# Usage
assay <- "Nanostring"
scale_factor <- 1000

obj <- run_tma_clustering(
  obj = subset(SeuObj, minor_celltype %in% c("Apoe+/C1q+ macrophage/APC (Dendritic)",
                                               "cDC (Dendritic)")),
  assay = assay,
  dims_use = 1:15,
  resolutions = c(0.3, 0.6),
  scale_factor = scale_factor
)

DimPlot(obj, 
        group.by = "cluster_regular_res03", 
        reduction = "UMAP_regular") & NoAxes() 



# just checking annotation - looks good! no action needed here
mk_general <- get_top_markers(
  obj = obj,
  cluster_col = "cluster_regular_res03",
  assay = assay,
  top_n = 5
)

plot_marker_dotplot(
  obj,
  cluster_col = "cluster_regular_res06",
  features = unique(mk_general$top_markers[, ]$feature),
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

table(obj$major_celltype,
      obj$cluster_regular_res06)
table(obj$cell_based_final,
      obj$cluster_regular_res06)

DimPlot(obj, 
        group.by = "minor_celltype", 
        reduction = "UMAP_regular", 
        cols = celltype_DCs_cols) & NoAxes() 


# qsave(obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_Den21940cells_clustered.qs")

obj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_Den21940cells_clustered.qs")
DimPlot(obj, 
        group.by = "minor_celltype", 
        reduction = "UMAP_regular", 
        cols = celltype_DCs_cols) & NoAxes() 

### Figure 3B -----------
DefaultAssay(obj) <- "Nanostring"
Seurat::FeaturePlot(
  object = obj,
  features = c("Itgax"), #pDC
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) &
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) &
NoAxes()

DefaultAssay(obj) <- "Nanostring"
Seurat::FeaturePlot(
  object = obj,
  features = c("Bst2"), #pDC
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) &
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) &
NoAxes()


DefaultAssay(obj) <- "Nanostring"
Seurat::FeaturePlot(
  object = obj,
  features = c("Cd8a"), #and/or Cd8b1
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) &
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) &
NoAxes()


DefaultAssay(obj) <- "Nanostring"
Seurat::FeaturePlot(
  object = obj,
  features = c("Lamp3", "Lcn2"), 
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) &
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) &
NoAxes()

DefaultAssay(obj) <- "Nanostring"
Seurat::FeaturePlot(
  object = obj,
  features = c("Apoe", "C1qa"), 
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) &
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) &
NoAxes()



### Figure 3C -----------
## ----------------------------------
## Shared plotting theme
## ----------------------------------

prop_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),

    axis.title.x = element_blank(),
    axis.title.y = element_blank(),

    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1
    ),

    axis.text.y = element_text(size = 14),

    axis.ticks = element_blank(),

    legend.title = element_blank(),
    legend.text = element_text(size = 12),

    strip.text.x = element_text(
      size = 16,
      face = "bold"
    )
  )

## ==================================
## Primary infection
## ==================================

primary_map <- c(
  "Neonate control" = "Neo Ctrl",
  "Neonate RSV infected (NO IFN)" = "Neo RSV",
  "Neonate IFN and RSV infected" = "Neo IFN+RSV",
  "Adult RSV infected" = "Adult RSV"
)

primary_order <- c(
  "Neo Ctrl",
  "Neo RSV",
  "Neo IFN+RSV",
  "Adult RSV"
)

obj_primary <- subset(
  obj,
  subset =
    minor_celltype != "other" &
    Group %in% names(primary_map)
)

obj_primary$Group_short <- recode(
  obj_primary$Group,
  !!!primary_map
)

obj_primary$Group_short <- factor(
  obj_primary$Group_short,
  levels = primary_order
)

p_primary <- dittoBarPlot(
  obj_primary,
  var = "minor_celltype",
  group.by = "Group_short",
  color.panel = celltype_DCs_cols,
  main = "Primary infection"
) +
  scale_x_discrete(limits = primary_order) +
  prop_theme

p_primary


## ==================================
## Secondary infection
## ==================================

reinf_map <- c(
  "Neonate control" = "Neo Ctrl",
  "Neonate (NO IFN) reinfected" = "Neo Reinf",
  "Neonate IFN and RSV reinfected" = "Neo IFN Reinf",
  "Adult reinfected" = "Adult Reinf"
)

reinf_order <- c(
  "Neo Ctrl",
  "Neo Reinf",
  "Neo IFN Reinf",
  "Adult Reinf"
)

obj_reinf <- subset(
  obj,
  subset =
    minor_celltype != "other" &
    Group %in% names(reinf_map)
)

obj_reinf$Group_short <- recode(
  obj_reinf$Group,
  !!!reinf_map
)

obj_reinf$Group_short <- factor(
  obj_reinf$Group_short,
  levels = reinf_order
)

p_reinf <- dittoBarPlot(
  obj_reinf,
  var = "minor_celltype",
  group.by = "Group_short",
  color.panel = celltype_DCs_cols,
  main = "Secondary infection"
) +
  scale_x_discrete(limits = reinf_order) +
  prop_theme

p_reinf


### Figure 3D -----------
# spatial plots - save it for later 

### Figure 3E -----------
# these volcano plots and barplots are all coming from here /mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/scripts/final_celltype/DEGs.R




# TESTING 
counts <- SeuObj@assays$Nanostring$counts
gene_panel <- rownames(counts)

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

DefaultAssay(SeuObj) <- "Nanostring"
Seurat::FeaturePlot(
  object = SeuObj,
  features = c("Cd8a"),
  slot = "data",
  reduction = "UMAP_regular",
  order = TRUE,
  raster = FALSE,
  #max.cutoff = "q99",
  pt.size = 0.4
  
) &
scale_color_gradientn(
  colors = c(
    "white",
    "#DEEBF7",
    "#9ECAE1",
    "#6BAED6",
    "#4292C6",
    "#FB8D59",
    "#FB6A4A",
    "#D7301F"
  )
) &
NoAxes()

