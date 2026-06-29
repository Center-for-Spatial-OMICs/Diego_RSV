## Description: Running donwstream analysis up to clustering for all slides (adult 1, neonate 1, adult 2, neoante 2)

## Load packages --------
library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)


## Help functions  
process_seurat_object <- function(Seu_obj,
                                  Downsize_to_test = FALSE,
                                  Downsize_n_cells = 100,
                                  n_insitutype_clusters = 10, 
                                  Metadata_dir,
                                  run_insitutype = TRUE,
                                  run_seuratworkflow = TRUE) {
  # Downsize to test the function
  if(Downsize_to_test == TRUE) {
    Seu_obj <- subset(x = Seu_obj, downsample = Downsize_n_cells)
    
  }
  
  # Load metadata based on conditions
  metadata <- fread(Metadata_dir)
  # Match cell ID in metadata
  metadata$cell = paste0(as.character(metadata$cell_ID), "_", metadata$fov)
  
  # Prepare Seurat object
  seurat_obj <- Seu_obj
  seurat_obj@meta.data$cell <- rownames(seurat_obj@meta.data)
  
  # Update metadata to only include matching cells
  metadata <- metadata[metadata$cell %in% seurat_obj@meta.data$cell, ]
  
  # First unsupervised clustering if run_insitutype is TRUE
  if (run_insitutype) {
    unsup <- insitutype(
      x = as.matrix(t(seurat_obj@assays$Nanostring$counts)),
      neg = metadata$nCount_negprobes,
      reference_profiles = NULL,
      bg = NULL,
      n_clusts = n_insitutype_clusters,
      n_phase1 = 200,
      n_phase2 = 500,
      n_phase3 = 2000,
      n_starts = 1,
      max_iters = 5
    )
    
    seurat_obj$insitutype_cluster <- 'other'
    seurat_obj$insitutype_cluster[names(unsup$clust)] <- unsup$clust
    names(seurat_obj@meta.data)[names(seurat_obj@meta.data) == "insitutype_cluster"] <- paste0("insitutype_cluster_", n_insitutype_clusters, "_clusters")
  }
  
  # Normalize, scale, and perform PCA and clustering if run_seuratworkflow is TRUE
  if (run_seuratworkflow) {
    seurat_obj <- seurat_obj %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      FindNeighbors(dims = 1:15) %>%
      FindClusters() %>% 
      RunUMAP(dims = 1:15)
    
    # Rename clusters if needed (optional)
    # names(seurat_obj@meta.data)[names(seurat_obj@meta.data) == "seurat_clusters"] <- "seurat_clusters_new_ROI"
  }
  
  # Initialize return list
  result <- list(seurat_object = seurat_obj)
  
  # Conditionally add unsupervised results if insitutype was run
  if (run_insitutype) {
    result$unsupervised_results <- unsup
  } else {
    result$unsupervised_results <- NULL
  }
  
  return(result)
}



## Define variables -----------
# Data dir paths
Adult_1_dir <- "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult1831439_08_11_2024_9_49_46_243/flatFiles/831439_IND_CR_FFPE_TMA_Adult1_STJ_N_R1"

Neonate_1_dir <- "/mnt/scratch2/CosMX/Diego/DiegoRSVNeonate1831439_08_11_2024_11_18_07_116/flatFiles/831439_IND_CR_FFPE_TMA_Neonate1_STJ_N_R1"

Adult_2_dir <- "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult2831439_08_11_2024_12_58_14_961/flatFiles/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1"

Neonate_2_dir <- "/mnt/scratch2/CosMX/Diego/DIegoRSVNeonate2831439_08_11_2024_13_26_50_655/flatFiles/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1"

# Seurat objects
sobj_list <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult12_Neonate12_seuobj_list.rds")

Adult_1_obj <- sobj_list$Adult_1
Neonate_1_obj <- sobj_list$Neonate_1

Adult_2_obj <- sobj_list$Adult_2
Neonate_2_obj <- sobj_list$Neonate_2

# DEG from GeoMx 
library(readr)
Adu1_Neo1_GeoMx_DEG_YL <- read_csv("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adu1_Neo1_GeoMx_DEG_YL.csv") %>% data.frame()

names(Adu1_Neo1_GeoMx_DEG_YL)[names(Adu1_Neo1_GeoMx_DEG_YL) == "...1"] <- "gene"

dim(Adu1_Neo1_GeoMx_DEG_YL)
head(Adu1_Neo1_GeoMx_DEG_YL)

deg_filt_neonate <- Adu1_Neo1_GeoMx_DEG_YL %>% filter(adj.P.Val <= 0.05 &
                                                   logFC > 0)
deg_filt_adult <- Adu1_Neo1_GeoMx_DEG_YL %>% filter(adj.P.Val <= 0.05 &
                                                     logFC < 0)


## Running clustering (Adult 1) ----------
Adult_1_obj <- subset(Adult_1_obj, TMA %in% "", invert = TRUE)
meta <- list.files(path = Adult_1_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)

Out_list <- process_seurat_object(Adult_1_obj, #already have run_seuratworkflow ran on it
                                    Metadata_dir = meta,
                                    # Downsize_to_test = TRUE,
                                    # Downsize_n_cells = 100,
                                    n_insitutype_clusters = 10, 
                                    run_insitutype = TRUE,
                                    run_seuratworkflow = FALSE)


Adult_1_obj <- Out_list$seurat_object
unsup_Adult_1 <- Out_list$unsupervised_results


## Running clustering (Neonate 1) ----------
Neonate_1_obj <- subset(Neonate_1_obj, TMA %in% "", invert = TRUE)
meta <- list.files(path = Neonate_1_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)

Out_list <- process_seurat_object(Neonate_1_obj, #already have run_seuratworkflow ran on it
                                  Metadata_dir = meta,
                                  # Downsize_to_test = TRUE,
                                  # Downsize_n_cells = 100,
                                  n_insitutype_clusters = 10, 
                                  run_insitutype = TRUE,
                                  run_seuratworkflow = FALSE)


Neonate_1_obj <- Out_list$seurat_object
unsup_Neonate_1 <- Out_list$unsupervised_results


## Running clustering (Adult 2) ----------
Adult_2_obj <- subset(Adult_2_obj, TMA %in% "", invert = TRUE)
meta <- list.files(path = Adult_2_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)

Out_list <- process_seurat_object(Adult_2_obj, #already have run_seuratworkflow ran on it
                                  Metadata_dir = meta,
                                  # Downsize_to_test = TRUE,
                                  # Downsize_n_cells = 100,
                                  n_insitutype_clusters = 10, 
                                  run_insitutype = TRUE,
                                  run_seuratworkflow = FALSE)


Adult_2_obj <- Out_list$seurat_object
unsup_Adult_2 <- Out_list$unsupervised_results


## Running clustering (Neonate 2) ----------
Neonate_2_obj <- subset(Neonate_2_obj, TMA %in% "", invert = TRUE)
meta <- list.files(path = Neonate_2_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)

Out_list <- process_seurat_object(Neonate_2_obj, #already have run_seuratworkflow ran on it
                                  Metadata_dir = meta,
                                  # Downsize_to_test = TRUE,
                                  # Downsize_n_cells = 100,
                                  n_insitutype_clusters = 10, 
                                  run_insitutype = TRUE,
                                  run_seuratworkflow = FALSE)


Neonate_2_obj <- Out_list$seurat_object
unsup_Neonate_2 <- Out_list$unsupervised_results


# ## In case we have not run the seurat workflow yet - just get the "data" slot for the downstream plots
# Adult_1_obj <- Adult_1_obj %>% NormalizeData()
# Adult_2_obj <- Adult_2_obj %>% NormalizeData()
# Neonate_2_obj <- Neonate_2_obj %>% NormalizeData() 
# Neonate_1_obj <- Neonate_2_obj %>% NormalizeData() 


# ## Plotting report (Unstructured Code) -----------------
# 
# ## Adult 1 example
# Start <- Sys.time()
# 
# obj <- Adult_1_obj %>% NormalizeData()
# unsup <- unsup_Adult_1
# Groups <- obj@meta.data$Group %>% table() %>% names()
# Groups <- Groups[c(1, 3, 2)]
# 
# obj@meta.data$insitutype_cluster_10_clusters <- factor(obj@meta.data$insitutype_cluster_10_clusters, levels =c(names(table(obj@meta.data$insitutype_cluster_10_clusters)))) 
# 
# obj@meta.data$Group <- factor(obj@meta.data$Group, levels =c(Groups)) 
# 
# pdf(file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Plots/Adult_1_report_GroupsAdded.pdf", width = 12, height = 6, )
# 
# ## Plot: Insisutype cluster quality 
# # cols <-
# #   c(
# #     '#8DD3C7',
# #     '#BEBADA',
# #     '#FB8072',
# #     '#80B1D3',
# #     '#FDB462',
# #     '#B3DE69',
# #     '#FCCDE5',
# #     '#D9D9D9',
# #     '#BC80BD',
# #     '#CCEBC5',
# #     '#FFED6F',
# #     '#E41A1C',
# #     '#377EB8',
# #     '#4DAF4A',
# #     '#984EA3',
# #     '#FF7F00',
# #     '#FFFF33',
# #     '#A65628',
# #     '#F781BF',
# #     '#999999'
# #   )
# cols = brewer.pal(10, 'Paired')
# cols <- cols[seq_along(unique(unsup$clust))]
# names(cols) <- unique(unsup$clust)
# fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup, col = cols[unsup$clust])
# #class(fp)
# print(fp)
# 
# 
# ## Plot: Cluster over spatial plot
# for (group in Groups) {
#   plt <- ImageDimPlot(subset(obj, Group == group), 
#                       group.by = 'insitutype_cluster_10_clusters', 
#                       cols = cols,
#                       axes = TRUE) +
#     ggtitle(paste0("Adult_1 - Group: ", group))
#   
#   print(plt)
# }
# 
# ## Plot: Clusters proportion by TMA
# plt <- dittoBarPlot(obj, 
#                     "insitutype_cluster_10_clusters", 
#                     group.by = "TMA", 
#                     color.panel = cols,
#                     main = "Cluster Proportion By TMA") 
# print(plt)
# 
# plt <- dittoBarPlot(obj, 
#                     "insitutype_cluster_10_clusters", 
#                     group.by = "Group", 
#                     color.panel = cols,
#                     main = "Cluster Proportion By Group") 
# print(plt)
# 
# ## Plot: Markers Dotplot 
# library(presto)
# obj@assays$RNA <- obj@assays$Nanostring
# markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
# 
# top_markers <- markers %>%
#   group_by(group) %>%
#   dplyr::filter(logFC > 0 & padj <= 0.05) %>%
#   dplyr::top_n(wt = logFC, n = 20) %>%
#   ungroup() %>% data.frame()
# 
# # write_csv(top_markers, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult_1_top20_markers_per_cluster.csv")
# 
# top_markers <- markers %>%
#   group_by(group) %>%
#   dplyr::filter(logFC > 0 & padj <= 0.05) %>%
#   dplyr::top_n(wt = logFC, n = 5) %>%
#   ungroup() %>%
#   pull(feature) %>%
#   unique()
# 
# insitutype_cluters_top_markers <- top_markers
# 
# 
# 
# ## Cluster Markers across clusters
# Idents(obj) <- "insitutype_cluster_10_clusters"
# plt <- DotPlot(obj, features = unique(c(insitutype_cluters_top_markers))) +
#   geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
#   scale_colour_viridis(option = "magma") +
#   guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle("Adult_1 - Cluster Markers across clusters")
# 
# print(plt)
# 
# 
# ## Cluster Markers across Groups
# Idents(obj) <- "Group"
# plt <- DotPlot(obj, features = unique(c(insitutype_cluters_top_markers))) +
#   geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
#   scale_colour_viridis(option = "magma") +
#   guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle("Adult_1 - Cluster Markers across groups")
# 
# print(plt)
# 
# 
# 
# ## Highly Specific Cluster makers in spatial 
# top_5_auc_markers <- markers[order(markers$auc, decreasing = TRUE), ]$feature[1:5]
# 
# ## Plot: Spatial plot 
# for (i in 1:length(top_5_auc_markers)) {
#   for (group in Groups) {
#     plt <- ImageFeaturePlot(subset(obj, Group == group), 
#                             border.color = 'white', 
#                             size = 0.5,
#                             #min.cutoff = 'q10',
#                             # max.cutoff = 'q90',
#                             features = c(top_5_auc_markers[i])) +
#       ggtitle(paste0("Adult_1 - Group: ", group, " - Top Cluster Markers ", top_5_auc_markers[i]))
#     
#     print(plt)
#   }
# }
# 
# ## Plot: Markers Vlnplot
# library(presto)
# obj@assays$RNA <- obj@assays$Nanostring
# markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
# insitutype_cluters_top_auc_markers <- markers[markers$auc > 0.70, ]$feature
# 
# # VlnPlot(obj, features = insitutype_cluters_top_auc_markers[1:6], group.by = "TMA", pt.size = 0)
# features <- insitutype_cluters_top_auc_markers
# num_features <- length(features)
# for (i in seq(1, num_features, by = 6)) {
#   # Select the current batch of features (up to 6)
#   current_features <- features[i:min(i + 5, num_features)]
#   
#   # Create the violin plot for the current batch
#   plt <- VlnPlot(obj, features = current_features, group.by = "Group", pt.size = 0)
#   
#   print(plt)
# }
# 
# 
# 
# ## GeoMx_DEG  %in% CosMx ---------
# # All 1k genes %in% GeoMx DEG
# rownames(obj) %in% deg_filt_adult$gene %>% table() 
# rownames(obj) %in% deg_filt_neonate$gene %>% table() 
# 
# # Insitutype cl. markers %in% GeoMx DEG
# insitutype_cluters_top_markers %in% deg_filt_neonate$gene %>% table() 
# insitutype_cluters_top_markers %in% deg_filt_adult$gene %>% table() 
# 
# CosMx_GeoMx_markers <- intersect(insitutype_cluters_top_markers, deg_filt_adult$gene)
# 
# CosMx_GeoMx_markers <- CosMx_GeoMx_markers[!grepl("Mt", CosMx_GeoMx_markers)][1:5]
# 
# ## Plot: Markers Dotplot
# Idents(obj) <- "insitutype_cluster_10_clusters"
# plt <- DotPlot(obj, features = unique(c(CosMx_GeoMx_markers))) +
#   geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
#   scale_colour_viridis(option = "magma") +
#   guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle("GeoMx DEG in CosMx (Adult_1)")
# 
# print(plt)
# 
# ## Plot: Spatial plot 
# for (i in 1:length(CosMx_GeoMx_markers)) {
#   for (group in Groups) {
#     plt <- ImageFeaturePlot(subset(obj, Group == group), 
#                             border.color = 'white', 
#                             size = 0.5,
#                             #min.cutoff = 'q10',
#                             # max.cutoff = 'q90',
#                             features = c(CosMx_GeoMx_markers[i])) +
#       ggtitle(paste0("Adult_1 - Group: ", group, " - GeoMx DEG in CosMx: ", CosMx_GeoMx_markers[i]))
#     
#     print(plt)
#   }
# }
# 
# 
# dev.off()
# End <- Sys.time()
# plot_time_Adult_1 <- print(End - Start)













## Plotting report (Refectored Code) -----------------

## Define {plot_report} function 
## Note: "deg_filt_adult" and "deg_filt_neonate" variables have being used but not declared in the arguments. They're defined early in this code


plot_report <- function(obj, unsup, Groups, obj_name) {
  Start <- Sys.time()
  
  obj <- obj %>% NormalizeData()
  
  # Ordering clusters
  obj@meta.data$insitutype_cluster_10_clusters <- factor(obj@meta.data$insitutype_cluster_10_clusters, levels = names(table(obj@meta.data$insitutype_cluster_10_clusters)))
  # Ordering groups
  obj@meta.data$Group <- factor(obj@meta.data$Group, levels = Groups)
  
  pdf(file = paste0("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Plots/", obj_name, "_report_GroupsAdded.pdf"), width = 12, height = 10)
  
  cols = brewer.pal(10, 'Paired')
  cols <- cols[seq_along(unique(unsup$clust))]
  names(cols) <- unique(unsup$clust)
  fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup, col = cols[unsup$clust])
  print(fp)
  
  ## Plot: Cluster over spatial plot
  for (group in Groups) {
    plt <- ImageDimPlot(subset(obj, Group == group), 
                        group.by = 'insitutype_cluster_10_clusters', 
                        cols = cols,
                        axes = TRUE) +
      ggtitle(paste0(obj_name, " - Group: ", group))
    
    print(plt)
  }
  
  ## Plot: Clusters proportion by TMA~Group
  for (group in Groups) {
  plt <- dittoBarPlot(subset(obj, Group == group), 
                      var = "insitutype_cluster_10_clusters", 
                      group.by = "TMA",
                      color.panel = cols,
                      main = paste0("Cluster Proportion By TMA - Group: ", group) )
  print(plt)
  }
  

  
  ## Plot: Markers Dotplot 
  library(presto)
  obj@assays$RNA <- obj@assays$Nanostring
  markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
  
  top_markers <- markers %>%
    group_by(group) %>%
    dplyr::filter(logFC > 0 & padj <= 0.05) %>%
    dplyr::top_n(wt = logFC, n = 20) %>%
    ungroup() %>% data.frame()
  
  # write_csv(top_markers, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult_1_top20_markers_per_cluster.csv")
  
  top_markers <- markers %>%
    group_by(group) %>%
    dplyr::filter(logFC > 0 & padj <= 0.05) %>%
    dplyr::top_n(wt = logFC, n = 5) %>%
    ungroup() %>%
    pull(feature) %>%
    unique()
  
  insitutype_cluters_top_markers <- top_markers
  
  
  
  ## Cluster Markers across clusters
  Idents(obj) <- "insitutype_cluster_10_clusters"
  plt <- DotPlot(obj, features = unique(c(insitutype_cluters_top_markers))) +
    geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
    scale_colour_viridis(option = "magma") +
    guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(obj_name," - Cluster Markers across clusters")
  
  print(plt)
  
  
  ## Cluster Markers across Groups
  Idents(obj) <- "Group"
  plt <- DotPlot(obj, features = unique(c(insitutype_cluters_top_markers))) +
    geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
    scale_colour_viridis(option = "magma") +
    guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(obj_name," - Cluster Markers across groups")
  
  print(plt)
  
  
  
  ## Highly Specific Cluster makers in spatial 
  top_5_auc_markers <- markers[order(markers$auc, decreasing = TRUE), ]$feature[1:5]
  
  ## Plot: Spatial plot 
  for (i in 1:length(top_5_auc_markers)) {
    for (group in Groups) {
      plt <- ImageFeaturePlot(subset(obj, Group == group), 
                              border.color = 'white', 
                              size = 0.5,
                              #min.cutoff = 'q10',
                              # max.cutoff = 'q90',
                              features = c(top_5_auc_markers[i])) +
        ggtitle(paste0(obj_name," - Group: ", group, " - Top Cluster Markers ", top_5_auc_markers[i]))
      
      print(plt)
    }
  }
  
  ## Plot: Markers Vlnplot
  library(presto)
  obj@assays$RNA <- obj@assays$Nanostring
  markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
  insitutype_cluters_top_auc_markers <- markers[markers$auc > 0.70, ]$feature
  
  # VlnPlot(obj, features = insitutype_cluters_top_auc_markers[1:6], group.by = "TMA", pt.size = 0)
  features <- insitutype_cluters_top_auc_markers
  num_features <- length(features)
  for (i in seq(1, num_features, by = 6)) {
    # Select the current batch of features (up to 6)
    current_features <- features[i:min(i + 5, num_features)]
    
    # Create the violin plot for the current batch
    plt <- VlnPlot(obj, features = current_features, group.by = "Group", pt.size = 0)
    
    print(plt)
  }
  
  
  
  ## GeoMx_DEG  %in% CosMx ---------
  # All 1k genes %in% GeoMx DEG
  rownames(obj) %in% deg_filt_adult$gene %>% table() 
  rownames(obj) %in% deg_filt_neonate$gene %>% table() 
  
  # Insitutype cl. markers %in% GeoMx DEG
  insitutype_cluters_top_markers %in% deg_filt_neonate$gene %>% table() 
  insitutype_cluters_top_markers %in% deg_filt_adult$gene %>% table() 
  
  CosMx_GeoMx_markers <- intersect(insitutype_cluters_top_markers, deg_filt_adult$gene)
  
  CosMx_GeoMx_markers <- CosMx_GeoMx_markers[!grepl("Mt", CosMx_GeoMx_markers)][1:5]
  
  ## Plot: Markers Dotplot
  Idents(obj) <- "insitutype_cluster_10_clusters"
  plt <- DotPlot(obj, features = unique(c(CosMx_GeoMx_markers))) +
    geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
    scale_colour_viridis(option = "magma") +
    guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(paste0("GeoMx DEG in CosMx (", obj_name,")"))
  
  print(plt)
  
  ## Plot: Spatial plot 
  for (i in 1:length(CosMx_GeoMx_markers)) {
    for (group in Groups) {
      plt <- ImageFeaturePlot(subset(obj, Group == group), 
                              border.color = 'white', 
                              size = 0.5,
                              #min.cutoff = 'q10',
                              # max.cutoff = 'q90',
                              features = c(CosMx_GeoMx_markers[i])) +
        ggtitle(paste0(obj_name, " - Group: ", group, " - GeoMx DEG in CosMx: ", CosMx_GeoMx_markers[i]))
      
      print(plt)
    }
  }
  
  dev.off()
  End <- Sys.time()
  plot_time <- print(End - Start)
  return(plot_time)
}

# Process Adult_1
Groups <- Adult_1_obj@meta.data$Group %>% table() %>% names()
Groups <- Groups[c(1, 3, 2)]
plot_time_Adult_1 <- plot_report(obj = Adult_1_obj, 
                                 unsup = unsup_Adult_1, 
                                 Groups = Groups,
                                 obj_name = "Adult_1")

# Process Adult_2
Groups <- Adult_2_obj@meta.data$Group %>% table() %>% names()
Groups <- Groups[c(1, 3, 2)]
plot_time_Adult_2 <- plot_report(obj = Adult_2_obj, 
                                 unsup = unsup_Adult_2, 
                                 Groups = Groups,
                                 obj_name = "Adult_2")

# Process Neonate_1
Groups <- Neonate_1_obj@meta.data$Group %>% table() %>% names()
Groups <- c("Neonate control",
            "Neonate RSV infected (NO IFN)",
            "Neonate (NO IFN) reinfected",
            "Neonate IFN and RSV infected",
            "Neonate IFN and RSV reinfected")

plot_time_Neonate_1 <- plot_report(obj = Neonate_1_obj, 
                                 unsup = unsup_Neonate_1,
                                 Groups = Groups,
                                 obj_name = "Neonate_1")

# Process Neonate_2
Groups <- Neonate_2_obj@meta.data$Group %>% table() %>% names()
Groups <- c("Neonate control",
            "Neonate RSV infected (NO IFN)",
            "Neonate (NO IFN) reinfected",
            "Neonate IFN and RSV infected",
            "Neonate IFN and RSV reinfected")
plot_time_Neonate_2 <- plot_report(obj = Neonate_2_obj, 
                                 unsup = unsup_Neonate_2, 
                                 Groups = Groups,
                                 obj_name = "Neonate_2")





# ## Adult_2 --------
# obj <- Adult_2_obj %>% NormalizeData()
# unsup <- unsup_Adult_2
# 
# pdf(file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Plots/Adult_2_report.pdf", width = 12, height = 6, )
# 
# ## Plot: Insisutype cluster quality 
# cols <-
#   c(
#     '#8DD3C7',
#     '#BEBADA',
#     '#FB8072',
#     '#80B1D3',
#     '#FDB462',
#     '#B3DE69',
#     '#FCCDE5',
#     '#D9D9D9',
#     '#BC80BD',
#     '#CCEBC5',
#     '#FFED6F',
#     '#E41A1C',
#     '#377EB8',
#     '#4DAF4A',
#     '#984EA3',
#     '#FF7F00',
#     '#FFFF33',
#     '#A65628',
#     '#F781BF',
#     '#999999'
#   )
# cols <- cols[seq_along(unique(unsup$clust))]
# names(cols) <- unique(unsup$clust)
# fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup, col = cols[unsup$clust])
# #class(fp)
# print(fp)
# 
# ## Plot: Cluster over spatial plot
# plt <- ImageDimPlot(obj, group.by = c('insitutype_cluster_10_clusters'), 
#                     cols = brewer.pal(10, 'Paired'),
#                     axes = TRUE) +
#   # coord_flip() + 
#   ggtitle("Adult_2")
# 
# print(plt)
# 
# ## Plot: Clusters proportion by TMA
# plt <- dittoBarPlot(obj, "insitutype_cluster_10_clusters", group.by = "TMA") 
# print(plt)
# 
# ## Plot: Markers Dotplot 
# library(presto)
# obj@assays$RNA <- obj@assays$Nanostring
# markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
# 
# top_markers <- markers %>%
#   group_by(group) %>%
#   dplyr::filter(logFC > 0 & padj <= 0.05) %>%
#   dplyr::top_n(wt = logFC, n = 20) %>%
#   ungroup() %>% data.frame()
# 
# # write_csv(top_markers, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult_2_top20_markers_per_cluster.csv")
# 
# 
# top_markers <- markers %>%
#   group_by(group) %>%
#   dplyr::filter(logFC > 0 & padj <= 0.05) %>%
#   dplyr::top_n(wt = logFC, n = 5) %>%
#   ungroup() %>%
#   pull(feature) %>%
#   unique()
# 
# insitutype_cluters_top_markers <- top_markers
# 
# # Ordering levels for plotting
# obj@meta.data$insitutype_cluster_10_clusters <- factor(obj@meta.data$insitutype_cluster_10_clusters, levels =c(names(table(obj@meta.data$insitutype_cluster_10_clusters)))) 
# 
# Idents(obj) <- "insitutype_cluster_10_clusters"
# plt <- DotPlot(obj, features = unique(c(insitutype_cluters_top_markers))) +
#   geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
#   scale_colour_viridis(option = "magma") +
#   guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle("Adult_2")
# 
# print(plt)
# 
# ## Plot: Markers Vlnplot
# library(presto)
# obj@assays$RNA <- obj@assays$Nanostring
# markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
# insitutype_cluters_top_auc_markers <- markers[markers$auc > 0.70, ]$feature
# 
# # VlnPlot(obj, features = insitutype_cluters_top_auc_markers[1:6], group.by = "TMA", pt.size = 0)
# features <- insitutype_cluters_top_auc_markers
# num_features <- length(features)
# for (i in seq(1, num_features, by = 6)) {
#   # Select the current batch of features (up to 6)
#   current_features <- features[i:min(i + 5, num_features)]
#   
#   # Create the violin plot for the current batch
#   plt <- VlnPlot(obj, features = current_features, group.by = "TMA", pt.size = 0)
#   
#   print(plt)
# }
# 
# 
# 
# ## GeoMx_DEG  %in% CosMx ---------
# # All 1k genes %in% GeoMx DEG
# rownames(obj) %in% deg_filt_adult$gene %>% table() 
# rownames(obj) %in% deg_filt_neonate$gene %>% table() 
# 
# # Insitutype cl. markers %in% GeoMx DEG
# insitutype_cluters_top_markers %in% deg_filt_neonate$gene %>% table() 
# insitutype_cluters_top_markers %in% deg_filt_adult$gene %>% table() 
# 
# CosMx_GeoMx_markers <- intersect(insitutype_cluters_top_markers, deg_filt_adult$gene)
# 
# 
# ## Plot: Markers Dotplot
# Idents(obj) <- "insitutype_cluster_10_clusters"
# plt <- DotPlot(obj, features = unique(c(CosMx_GeoMx_markers))) +
#   geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
#   scale_colour_viridis(option = "magma") +
#   guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle("GeoMx DEG in CosMx (Adult_2)")
# 
# print(plt)
# 
# ## Plot: Spatial plot 
# 
# 
# for(i in 1:length(CosMx_GeoMx_markers)) {
#   plt <- ImageFeaturePlot(#subset(obj, TMA == TMA_values[1]), 
#     obj,
#     border.color = 'white', 
#     size = 0.5,
#     #min.cutoff = 'q10',
#     # max.cutoff = 'q90',
#     features = c(CosMx_GeoMx_markers[i]))
#   
#   print(plt)
#   
# }
# 
# dev.off()
# 
# 
# 
# ## Neonate_1 --------
# obj <- Neonate_1_obj %>% NormalizeData()
# unsup <- unsup_Neonate_1
# 
# pdf(file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Plots/Neonate_1_report.pdf", width = 12, height = 6, )
# 
# ## Plot: Insisutype cluster quality 
# cols <-
#   c(
#     '#8DD3C7',
#     '#BEBADA',
#     '#FB8072',
#     '#80B1D3',
#     '#FDB462',
#     '#B3DE69',
#     '#FCCDE5',
#     '#D9D9D9',
#     '#BC80BD',
#     '#CCEBC5',
#     '#FFED6F',
#     '#E41A1C',
#     '#377EB8',
#     '#4DAF4A',
#     '#984EA3',
#     '#FF7F00',
#     '#FFFF33',
#     '#A65628',
#     '#F781BF',
#     '#999999'
#   )
# cols <- cols[seq_along(unique(unsup$clust))]
# names(cols) <- unique(unsup$clust)
# fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup, col = cols[unsup$clust])
# #class(fp)
# print(fp)
# 
# ## Plot: Cluster over spatial plot
# plt <- ImageDimPlot(obj, group.by = c('insitutype_cluster_10_clusters'), 
#                     cols = brewer.pal(10, 'Paired'),
#                     axes = TRUE) +
#   # coord_flip() + 
#   ggtitle("Neonate_1")
# 
# print(plt)
# 
# ## Plot: Clusters proportion by TMA
# plt <- dittoBarPlot(obj, "insitutype_cluster_10_clusters", group.by = "TMA") 
# print(plt)
# 
# ## Plot: Markers Dotplot 
# library(presto)
# obj@assays$RNA <- obj@assays$Nanostring
# markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
# 
# top_markers <- markers %>%
#   group_by(group) %>%
#   dplyr::filter(logFC > 0 & padj <= 0.05) %>%
#   dplyr::top_n(wt = logFC, n = 20) %>%
#   ungroup() %>% data.frame()
# 
# # write_csv(top_markers, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Neonate_1_top20_markers_per_cluster.csv")
# 
# top_markers <- markers %>%
#   group_by(group) %>%
#   dplyr::filter(logFC > 0 & padj <= 0.05) %>%
#   dplyr::top_n(wt = logFC, n = 5) %>%
#   ungroup() %>%
#   pull(feature) %>%
#   unique()
# 
# insitutype_cluters_top_markers <- top_markers
# 
# # Ordering levels for plotting
# obj@meta.data$insitutype_cluster_10_clusters <- factor(obj@meta.data$insitutype_cluster_10_clusters, levels =c(names(table(obj@meta.data$insitutype_cluster_10_clusters)))) 
# 
# Idents(obj) <- "insitutype_cluster_10_clusters"
# plt <- DotPlot(obj, features = unique(c(insitutype_cluters_top_markers))) +
#   geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
#   scale_colour_viridis(option = "magma") +
#   guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle("Neonate_1")
# 
# print(plt)
# 
# ## Plot: Markers Vlnplot
# library(presto)
# obj@assays$RNA <- obj@assays$Nanostring
# markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
# insitutype_cluters_top_auc_markers <- markers[markers$auc > 0.70, ]$feature
# 
# # VlnPlot(obj, features = insitutype_cluters_top_auc_markers[1:6], group.by = "TMA", pt.size = 0)
# features <- insitutype_cluters_top_auc_markers
# num_features <- length(features)
# for (i in seq(1, num_features, by = 6)) {
#   # Select the current batch of features (up to 6)
#   current_features <- features[i:min(i + 5, num_features)]
#   
#   # Create the violin plot for the current batch
#   plt <- VlnPlot(obj, features = current_features, group.by = "TMA", pt.size = 0)
#   
#   print(plt)
# }
# 
# 
# 
# ## GeoMx_DEG  %in% CosMx ---------
# # All 1k genes %in% GeoMx DEG
# rownames(obj) %in% deg_filt_adult$gene %>% table() 
# rownames(obj) %in% deg_filt_neonate$gene %>% table() 
# 
# # Insitutype cl. markers %in% GeoMx DEG
# insitutype_cluters_top_markers %in% deg_filt_neonate$gene %>% table() 
# insitutype_cluters_top_markers %in% deg_filt_adult$gene %>% table() 
# 
# CosMx_GeoMx_markers <- intersect(insitutype_cluters_top_markers, deg_filt_neonate$gene)
# 
# 
# ## Plot: Markers Dotplot
# Idents(obj) <- "insitutype_cluster_10_clusters"
# plt <- DotPlot(obj, features = unique(c(CosMx_GeoMx_markers))) +
#   geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
#   scale_colour_viridis(option = "magma") +
#   guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle("GeoMx DEG in CosMx (Neonate_1)")
# 
# print(plt)
# 
# ## Plot: Spatial plot 
# 
# 
# for(i in 1:length(CosMx_GeoMx_markers)) {
#   plt <- ImageFeaturePlot(#subset(obj, TMA == TMA_values[1]), 
#     obj,
#     border.color = 'white', 
#     size = 0.5,
#     #min.cutoff = 'q10',
#     # max.cutoff = 'q90',
#     features = c(CosMx_GeoMx_markers[i]))
#   
#   print(plt)
#   
# }
# 
# dev.off()
# 
# 
# 
# 
# 
# ## Neonate_2 --------
# obj <- Neonate_2_obj %>% NormalizeData()
# unsup <- unsup_Neonate_2
# 
# pdf(file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Plots/Neonate_2_report.pdf", width = 12, height = 6, )
# 
# ## Plot: Insisutype cluster quality 
# cols <-
#   c(
#     '#8DD3C7',
#     '#BEBADA',
#     '#FB8072',
#     '#80B1D3',
#     '#FDB462',
#     '#B3DE69',
#     '#FCCDE5',
#     '#D9D9D9',
#     '#BC80BD',
#     '#CCEBC5',
#     '#FFED6F',
#     '#E41A1C',
#     '#377EB8',
#     '#4DAF4A',
#     '#984EA3',
#     '#FF7F00',
#     '#FFFF33',
#     '#A65628',
#     '#F781BF',
#     '#999999'
#   )
# cols <- cols[seq_along(unique(unsup$clust))]
# names(cols) <- unique(unsup$clust)
# fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup, col = cols[unsup$clust])
# #class(fp)
# print(fp)
# 
# ## Plot: Cluster over spatial plot
# plt <- ImageDimPlot(obj, group.by = c('insitutype_cluster_10_clusters'), 
#                     cols = brewer.pal(10, 'Paired'),
#                     axes = TRUE) +
#   # coord_flip() + 
#   ggtitle("Neonate_2")
# 
# print(plt)
# 
# ## Plot: Clusters proportion by TMA
# plt <- dittoBarPlot(obj, "insitutype_cluster_10_clusters", group.by = "TMA") 
# print(plt)
# 
# ## Plot: Markers Dotplot 
# library(presto)
# obj@assays$RNA <- obj@assays$Nanostring
# markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
# 
# top_markers <- markers %>%
#   group_by(group) %>%
#   dplyr::filter(logFC > 0 & padj <= 0.05) %>%
#   dplyr::top_n(wt = logFC, n = 20) %>%
#   ungroup() %>% data.frame()
# 
# # write_csv(top_markers, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Neonate_2_top20_markers_per_cluster.csv")
# 
# top_markers <- markers %>%
#   group_by(group) %>%
#   dplyr::filter(logFC > 0 & padj <= 0.05) %>%
#   dplyr::top_n(wt = logFC, n = 5) %>%
#   ungroup() %>%
#   pull(feature) %>%
#   unique()
# 
# insitutype_cluters_top_markers <- top_markers
# 
# # Ordering levels for plotting
# obj@meta.data$insitutype_cluster_10_clusters <- factor(obj@meta.data$insitutype_cluster_10_clusters, levels =c(names(table(obj@meta.data$insitutype_cluster_10_clusters)))) 
# 
# Idents(obj) <- "insitutype_cluster_10_clusters"
# plt <- DotPlot(obj, features = unique(c(insitutype_cluters_top_markers))) +
#   geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
#   scale_colour_viridis(option = "magma") +
#   guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle("Neonate_2")
# 
# print(plt)
# 
# ## Plot: Markers Vlnplot
# library(presto)
# obj@assays$RNA <- obj@assays$Nanostring
# markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
# insitutype_cluters_top_auc_markers <- markers[markers$auc > 0.70, ]$feature
# 
# # VlnPlot(obj, features = insitutype_cluters_top_auc_markers[1:6], group.by = "TMA", pt.size = 0)
# features <- insitutype_cluters_top_auc_markers
# num_features <- length(features)
# for (i in seq(1, num_features, by = 6)) {
#   # Select the current batch of features (up to 6)
#   current_features <- features[i:min(i + 5, num_features)]
#   
#   # Create the violin plot for the current batch
#   plt <- VlnPlot(obj, features = current_features, group.by = "TMA", pt.size = 0)
#   
#   print(plt)
# }
# 
# 
# 
# ## GeoMx_DEG  %in% CosMx ---------
# # All 1k genes %in% GeoMx DEG
# rownames(obj) %in% deg_filt_adult$gene %>% table() 
# rownames(obj) %in% deg_filt_neonate$gene %>% table() 
# 
# # Insitutype cl. markers %in% GeoMx DEG
# insitutype_cluters_top_markers %in% deg_filt_neonate$gene %>% table() 
# insitutype_cluters_top_markers %in% deg_filt_adult$gene %>% table() 
# 
# CosMx_GeoMx_markers <- intersect(insitutype_cluters_top_markers, deg_filt_neonate$gene)
# 
# 
# ## Plot: Markers Dotplot
# Idents(obj) <- "insitutype_cluster_10_clusters"
# plt <- DotPlot(obj, features = unique(c(CosMx_GeoMx_markers))) +
#   geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
#   scale_colour_viridis(option = "magma") +
#   guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggtitle("GeoMx DEG in CosMx (Neonate_2)")
# 
# print(plt)
# 
# ## Plot: Spatial plot 
# 
# #DefaultBoundary(obj[["FOV"]]) <- 'segmentation'
# # DefaultBoundary(obj[["FOV"]]) <- 'centroids'
# 
# for(i in 1:length(CosMx_GeoMx_markers)) {
#   plt <- ImageFeaturePlot(#subset(obj, TMA == TMA_values[1]),
#     obj,
#     border.color = 'white',
#     size = 0.5,
#     #min.cutoff = 'q60',
#     #max.cutoff = 'q90',
#     features = c(CosMx_GeoMx_markers[i]))
# 
#   print(plt)
#   
#   
# }
#   
# dev.off()


## STILL IN TEST - gene exp visualization approach 
# for(i in 1:length(CosMx_GeoMx_markers)) {
#   gene_name <- CosMx_GeoMx_markers[i]  
#   cells <- colnames(obj)[GetAssayData(obj, slot = "data")[gene_name, ] > 4] # should change depending on the given gene 
#   obj_sub <- subset(x = obj, cells = cells)
#   
#   plt <- ImageFeaturePlot(
#     obj_sub,
#     border.color = 'white',
#     size = 0.5,
#     #min.cutoff = 'q60',
#     #max.cutoff = 'q90',
#     features = c(CosMx_GeoMx_markers[i]))
# 
#   print(plt)  
# 
# }
  







