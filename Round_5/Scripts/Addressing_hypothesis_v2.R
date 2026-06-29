## ====== Date: Sep 04, 2025  ===== 


### Load and define functions ---------

library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(qs)

options(future.globals.maxSize = 99999 * 1024^2)
Sys.setenv("VROOM_CONNECTION_SIZE" = 999999 * 2)

### Load and define variables ---------
Sobj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj_list_celltype_niche.qs") #objects are separated by slides, not TMAs. Although I have a TMA columns to retrieve this info on their @meta.data 


comparisons <- list(
  "Comparison" = list(
    group_1 = "Neonate RSV infected (NO IFN)",
    group_2 = "Neonate IFN and RSV infected",
    group_3 = "Adult RSV infected"
  )
)

# markers from Diego 
# Gene names matching the 1k mouse gene panel
cell_types <- list(
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


cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- Sobj_list[[1]]$major_celltype %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)

library(viridis)
niche_cols <- viridis(5)
factor_levels <- Sobj_list[[1]]$niches_5 %>% table() %>% names()
niche_cols <- setNames(niche_cols[1:length(factor_levels)], factor_levels)


### Basic vis. of this mouse-cohort -------
Sobj_list$Adult_1$TMA_2
Sobj_list$Adult_1$TMA
Sobj_list$Adult_1$major_celltype

# adult 1
table(Sobj_list$Adult_1$TMA, useNA = "always") #okay 
dittoBarPlot(
  Sobj_list$Adult_1,
  var = "major_celltype",
  group.by = "TMA",
  color.panel = celltype_cols,
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)

table(Sobj_list$Adult_1$TMA_2, useNA = "always") #okay 
dittoBarPlot(
  Sobj_list$Adult_1,
  var = "major_celltype",
  group.by = "TMA_2",
  color.panel = celltype_cols,
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
) & NoLegend()


# adult 2
table(Sobj_list$Adult_2$TMA, useNA = "always") #okay 
dittoBarPlot(
  Sobj_list$Adult_2,
  var = "major_celltype",
  group.by = "TMA",
  color.panel = celltype_cols,
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)& NoLegend()

table(Sobj_list$Adult_2$TMA_2, useNA = "always") #okay 
dittoBarPlot(
  Sobj_list$Adult_2,
  var = "major_celltype",
  group.by = "TMA_2",
  color.panel = celltype_cols,
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
) & NoLegend()



# neonate 1
table(Sobj_list$Neonate_1$TMA, useNA = "always") #okay 
dittoBarPlot(
  Sobj_list$Neonate_1,
  var = "major_celltype",
  group.by = "TMA",
  color.panel = celltype_cols,
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
) & NoLegend()

table(Sobj_list$Neonate_1$TMA_2, useNA = "always") #okay 
dittoBarPlot(
  Sobj_list$Neonate_1,
  var = "major_celltype",
  group.by = "TMA_2",
  color.panel = celltype_cols,
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
) & NoLegend()


# neonate 2
table(Sobj_list$Neonate_2$TMA, useNA = "always") #okay 
dittoBarPlot(
  Sobj_list$Neonate_2,
  var = "major_celltype",
  group.by = "TMA",
  color.panel = celltype_cols,
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)& NoLegend()

table(Sobj_list$Neonate_2$TMA_2, useNA = "always") #okay 
dittoBarPlot(
  Sobj_list$Neonate_2,
  var = "major_celltype",
  group.by = "TMA_2",
  color.panel = celltype_cols,
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
) & NoLegend()




all_meta <- lapply(Sobj_list, function(obj) {
  meta <- obj@meta.data
  return(meta)
})
all_meta <- do.call(plyr::rbind.fill, all_meta)

# nFeature_RNA ~ gropups 
ggplot(all_meta, aes(x = nFeature_RNA, y = TMA, fill = TMA_2)) +
  geom_violin() +  # Violin plot without trimming tails
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 3,
    color = "black",
    fill = "white"
  ) +
  xlim(0, 500) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, color = "black", face = "bold")) + NoLegend() +
  facet_grid(~ Group)


# nCount_RNA ~ gropups 
ggplot(all_meta, aes(x = nCount_RNA, y = TMA, fill = TMA_2)) +
  geom_violin() +  # Violin plot without trimming tails
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 3,
    color = "black",
    fill = "white"
  ) +
  xlim(0, 1000) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, color = "black", face = "bold")) + NoLegend() +
  facet_grid(~ Group)




### Spatial feature plot 

# Separate slides into TMAs-wise objects

# adult 1
Sobj_list_TMAwise <- SplitObject(Sobj_list$Adult_1, split.by = "TMA")

library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)

# Apply DefaultBoundary setting
Sobj_list_TMAwise <- lapply(Sobj_list_TMAwise, function(x) {
  DefaultBoundary(x[["FOV"]]) <- "centroids"
  return(x)
})

# Get unique groups from the combined metadata of all objects
groups <- unique(sapply(Sobj_list_TMAwise, function(x) x@meta.data$Group[1]))

plot_list <- list()

for (grp in groups) {
  # Filter objects that belong to this group
  objs_in_group <- Sobj_list_TMAwise[sapply(Sobj_list_TMAwise, function(x) x@meta.data$Group[1] == grp)]
  
  # For each object in the group, create a plot for the specified features
  plots_per_obj <- lapply(objs_in_group, function(obj) {
    ImageFeaturePlot(obj, features = c('nCount_RNA', 'nFeature_RNA'), max.cutoff = c("q99")) &
      scale_fill_viridis() &
      coord_flip() &
      ggtitle(obj@meta.data$Group[1]) &
      theme(plot.title = element_text(size = 14, hjust = 0))
  })
  
  # Store the list of plots per group
  plot_list[[grp]] <- plots_per_obj
}

# To view plots for a particular group, e.g.:
wrap_plots(plot_list[["Adult Control"]], ncol = 5)

wrap_plots(plot_list[["Adult RSV infected"]], ncol = 5)

wrap_plots(plot_list[["Adult reinfected"]], ncol = 5)




# adult 2
Sobj_list_TMAwise <- SplitObject(Sobj_list$Adult_2, split.by = "TMA")

library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)

# Apply DefaultBoundary setting
Sobj_list_TMAwise <- lapply(Sobj_list_TMAwise, function(x) {
  DefaultBoundary(x[["FOV"]]) <- "centroids"
  return(x)
})

# Get unique groups from the combined metadata of all objects
groups <- unique(sapply(Sobj_list_TMAwise, function(x) x@meta.data$Group[1]))

plot_list <- list()

for (grp in groups) {
  # Filter objects that belong to this group
  objs_in_group <- Sobj_list_TMAwise[sapply(Sobj_list_TMAwise, function(x) x@meta.data$Group[1] == grp)]
  
  # For each object in the group, create a plot for the specified features
  plots_per_obj <- lapply(objs_in_group, function(obj) {
    ImageFeaturePlot(obj, features = c('nCount_RNA', 'nFeature_RNA'), max.cutoff = c("q99")) &
      scale_fill_viridis() &
      coord_flip() &
      ggtitle(obj@meta.data$Group[1]) &
      theme(plot.title = element_text(size = 14, hjust = 0))
  })
  
  # Store the list of plots per group
  plot_list[[grp]] <- plots_per_obj
}

# To view plots for a particular group, e.g.:
wrap_plots(plot_list[["Adult Control"]], ncol = 4)

wrap_plots(plot_list[["Adult RSV infected"]], ncol = 5)

wrap_plots(plot_list[["Adult reinfected"]], ncol = 4)





# neonate 1
Sobj_list_TMAwise <- SplitObject(Sobj_list$Neonate_1, split.by = "TMA")

library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)

# Apply DefaultBoundary setting
Sobj_list_TMAwise <- lapply(Sobj_list_TMAwise, function(x) {
  DefaultBoundary(x[["FOV"]]) <- "centroids"
  return(x)
})

# Get unique groups from the combined metadata of all objects
groups <- unique(sapply(Sobj_list_TMAwise, function(x) x@meta.data$Group[1]))

plot_list <- list()

for (grp in groups) {
  # Filter objects that belong to this group
  objs_in_group <- Sobj_list_TMAwise[sapply(Sobj_list_TMAwise, function(x) x@meta.data$Group[1] == grp)]
  
  # For each object in the group, create a plot for the specified features
  plots_per_obj <- lapply(objs_in_group, function(obj) {
    ImageFeaturePlot(obj, features = c('nCount_RNA', 'nFeature_RNA'), max.cutoff = c("q99")) &
      scale_fill_viridis() &
      coord_flip() &
      ggtitle(obj@meta.data$Group[1]) &
      theme(plot.title = element_text(size = 14, hjust = 0))
  })
  
  # Store the list of plots per group
  plot_list[[grp]] <- plots_per_obj
}

# To view plots for a particular group, e.g.:
wrap_plots(plot_list[["Neonate control"]], ncol = 4)

wrap_plots(plot_list[["Neonate RSV infected (NO IFN)"]], ncol = 3)

wrap_plots(plot_list[["Neonate IFN and RSV infected"]], ncol = 5)

wrap_plots(plot_list[["Neonate (NO IFN) reinfected"]], ncol = 4)

wrap_plots(plot_list[["Neonate IFN and RSV reinfected"]], ncol = 3)




# neonate 2
Sobj_list_TMAwise <- SplitObject(Sobj_list$Neonate_2, split.by = "TMA")

library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)

# Apply DefaultBoundary setting
Sobj_list_TMAwise <- lapply(Sobj_list_TMAwise, function(x) {
  DefaultBoundary(x[["FOV"]]) <- "centroids"
  return(x)
})

# Get unique groups from the combined metadata of all objects
groups <- unique(sapply(Sobj_list_TMAwise, function(x) x@meta.data$Group[1]))

plot_list <- list()

for (grp in groups) {
  # Filter objects that belong to this group
  objs_in_group <- Sobj_list_TMAwise[sapply(Sobj_list_TMAwise, function(x) x@meta.data$Group[1] == grp)]
  
  # For each object in the group, create a plot for the specified features
  plots_per_obj <- lapply(objs_in_group, function(obj) {
    ImageFeaturePlot(obj, features = c('nCount_RNA', 'nFeature_RNA'), max.cutoff = c("q99")) &
      scale_fill_viridis() &
      coord_flip() &
      ggtitle(obj@meta.data$Group[1]) &
      theme(plot.title = element_text(size = 14, hjust = 0))
  })
  
  # Store the list of plots per group
  plot_list[[grp]] <- plots_per_obj
}

# To view plots for a particular group, e.g.:
wrap_plots(plot_list[["Neonate control"]], ncol = 4)

wrap_plots(plot_list[["Neonate RSV infected (NO IFN)"]], ncol = 2)

wrap_plots(plot_list[["Neonate IFN and RSV infected"]], ncol = 1)

wrap_plots(plot_list[["Neonate (NO IFN) reinfected"]], ncol = 4)

wrap_plots(plot_list[["Neonate IFN and RSV reinfected"]], ncol = 5)

wrap_plots









