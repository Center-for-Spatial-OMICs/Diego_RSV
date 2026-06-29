library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)

# Helper Functions 
Metadata_dir <- "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult2831439_08_11_2024_12_58_14_961/flatFiles/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1_metadata_file.csv.gz"

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


# # #   Clustering CosMx data - downstream analysis round 1 # # #

setwd("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1")

Adult_2_dir <- "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult2831439_08_11_2024_12_58_14_961/flatFiles/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1"

Neonate_2_dir <- "/mnt/scratch2/CosMX/Diego/DIegoRSVNeonate2831439_08_11_2024_13_26_50_655/flatFiles/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1"


### Adult_2 ----------------

## Create a list of Seurat objects by FOV/TMA - accordingly to Hannah ppt slides
Adult_2_obj <- LoadNanostring(Adult_2_dir, fov = 'FOV')

Adult_2_tx <- fread(paste0(Adult_2_dir, "/", "831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1_tx_file.csv.gz"))
Adult_2_tx <- as.data.frame(Adult_2_tx)

Adult_2_meta <- fread(paste0(Adult_2_dir, "/", "831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1_metadata_file.csv.gz"))

Adult_2_obj$fov <- Adult_2_meta$fov
Adult_2_obj$Area <- Adult_2_meta$Area
Adult_2_obj$Area.um2 <- Adult_2_meta$Area.um2

# match cell id in meta obj
Adult_2_meta$cell = paste0(as.character(Adult_2_meta$cell_ID), "_", Adult_2_meta$fov)

# remove control and system probes
which <- grep('Neg*|Syst*', rownames(Adult_2_obj))
Adult_2_obj <- Adult_2_obj[-c(which), ]

# Remove low quality cells
Adult_2_obj <- subset(Adult_2_obj, nCount_Nanostring > 10 & nFeature_Nanostring > 10 & Area.um2 < 500)

#Update meta
Adult_2_meta <- Adult_2_meta[Adult_2_meta$cell %in% colnames(Adult_2_obj), ]
dim(Adult_2_meta) #104893     58

# Checking duplicates on TMA names
TMA_names <- c(
  'Ap16-3781 A1',
  'Ap16-3780 A1',
  'Ap16-3780 A1',
  'Ap16-3894 A1',
  'Ap16-3779 A1',
  'Ap16-3779 A1',
  'Ap16-3781 A1',
  'Ap17-0331 A1',
  'Ap16-3896 A1',
  'Ap16-3896 A1',
  'Ap16-3894 A1',
  'Ap17-0333 A1',
  'Ap17-0332 A1',
  'Ap17-0332 A1',
  'Ap17-0331 A1',
  'Ap17-0333 A1'
)

table(TMA_names)

# Classify different embryos by FOV numbers
Adult_2_obj$TMA <- ""
Adult_2_obj[[]] <- Adult_2_obj[[]] %>%
  mutate(TMA = case_when(
    fov %in% 1:10 ~ 'Ap16-3781 A1 1',
    fov %in% 38:49 ~ 'Ap16-3780 A1 1',
    fov %in% 11:19 ~ 'Ap16-3780 A1 2',
    fov %in% 59:74 ~ 'Ap16-3894 A1 1',
    fov %in% 20:28 ~ 'Ap16-3779 A1 1',
    fov %in% 50:58 ~ 'Ap16-3779 A1 2',
    fov %in% 29:37 ~ 'Ap16-3781 A1 2',
    fov %in% 75:83 ~ 'Ap17-0331 A1 1',
    fov %in% 114:128 ~ 'Ap16-3896 A1 1',
    fov %in% 84:99 ~ 'Ap16-3896 A1 2',
    fov %in% 100:113 ~ 'Ap16-3894 A1 2',
    fov %in% 129:137 ~ 'Ap17-0333 A1 1',
    fov %in% 175:188 ~ 'Ap17-0332 A1 1',
    fov %in% 138:174 ~ 'Ap17-0332 A1 2',
    fov %in% 147:161 ~ 'Ap17-0331 A1 2',
    fov %in% 162:173 ~ 'Ap17-0333 A1 2',
    
    TRUE ~ TMA  # Retain existing value or specify a default if needed
  ))



# Create list of embryos per run (break down of ProcessObj())
obj_list <- list()
TMA_values <- names(table(Adult_2_obj$TMA))
TMA_values <- TMA_values[-1]
for (i in TMA_values) {
  # Subset the Seurat object based on the current TMA value
  subset_seurat <- subset(Adult_2_obj, TMA == i)
  
  # Add the subsetted Seurat object to the list
  obj_list[[i]] <- subset_seurat
}

# Set default boundary as segmentation for each object in the list
for (i in seq_along(obj_list)) {
  DefaultBoundary(obj_list[[i]][['FOV']]) <- 'segmentation'
}

Adult_2_Sobj_list <- obj_list
print(Adult_2_Sobj_list)

Adult_2_Sobj_flat <- Adult_2_obj



## Processing with Seurat Workflow and insitutype clustering
Out_list <- process_seurat_object(Adult_2_Sobj_flat,
                                  Metadata_dir = "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult2831439_08_11_2024_12_58_14_961/flatFiles/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1_metadata_file.csv.gz",
                                  # Downsize_to_test = TRUE,
                                  # Downsize_n_cells = 100,
                                  n_insitutype_clusters = 10, 
                                  run_insitutype = TRUE,
                                  run_seuratworkflow = TRUE)

Adult_2_Sobj_flat_proc <- Out_list$seurat_object
unsup <- Out_list$unsupervised_results

  
## Visualization 
obj <- Adult_2_Sobj_flat_proc
names(obj@meta.data)[names(obj@meta.data) == "seurat_clusters"] <- "seurat_clusters_all_FOVs"

## Plot: nCount nFeature
# TMA_values <- names(table(obj$TMA))
ImageFeaturePlot(#subset(obj, TMA == TMA_values[1]), 
  obj, 
  border.color = 'white', 
  size = 1,
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  features = c('nCount_Nanostring'))

hist(obj$nCount_Nanostring)
median(obj$nCount_Nanostring)


ImageFeaturePlot(#subset(obj, TMA == TMA_values[1]), 
  obj, 
  border.color = 'white', 
  size = 1,
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  features = c('nFeature_Nanostring'))

hist(obj$nFeature_Nanostring)
median(obj$nFeature_Nanostring)


## Plot: Clusters over spatial plot
ImageDimPlot(obj, group.by = c('insitutype_cluster_10_clusters'), 
             cols = brewer.pal(10, 'Paired'),
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("insitutype_cluster_10_clusters")

cols <-
  c(
    '#8DD3C7',
    '#BEBADA',
    '#FB8072',
    '#80B1D3',
    '#FDB462',
    '#B3DE69',
    '#FCCDE5',
    '#D9D9D9',
    '#BC80BD',
    '#CCEBC5',
    '#FFED6F',
    '#E41A1C',
    '#377EB8',
    '#4DAF4A',
    '#984EA3',
    '#FF7F00',
    '#FFFF33',
    '#A65628',
    '#F781BF',
    '#999999'
  )
cols <- cols[seq_along(unique(unsup$clust))]
names(cols) <- unique(unsup$clust)
fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup, col = cols[unsup$clust])
#class(fp)
print(fp)

## Re-run only insitutype 
Out_list_2 <- process_seurat_object(Adult_2_Sobj_flat_proc, #already have run_seuratworkflow ran on it
                                    Metadata_dir = "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult2831439_08_11_2024_12_58_14_961/flatFiles/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1_metadata_file.csv.gz",
                                    # Downsize_to_test = TRUE,
                                    # Downsize_n_cells = 100,
                                    n_insitutype_clusters = 5, 
                                    run_insitutype = TRUE,
                                    run_seuratworkflow = FALSE)


Adult_2_Sobj_flat_proc_2 <- Out_list_2$seurat_object
unsup_2 <- Out_list_2$unsupervised_results

# save(Adult_2_Sobj_flat_proc_2, unsup_2, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Clustering_Adult_2.rda")

obj <- Adult_2_Sobj_flat_proc_2
names(obj@meta.data)[names(obj@meta.data) == "seurat_clusters"] <- "seurat_clusters_all_FOVs"

cols <-
  c(
    '#8DD3C7',
    '#BEBADA',
    '#FB8072',
    '#80B1D3',
    '#FDB462',
    '#B3DE69',
    '#FCCDE5',
    '#D9D9D9',
    '#BC80BD',
    '#CCEBC5',
    '#FFED6F',
    '#E41A1C',
    '#377EB8',
    '#4DAF4A',
    '#984EA3',
    '#FF7F00',
    '#FFFF33',
    '#A65628',
    '#F781BF',
    '#999999'
  )
cols <- cols[seq_along(unique(unsup_2$clust))]
names(cols) <- unique(unsup_2$clust)
fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup_2, col = cols[unsup_2$clust])
#class(fp)
print(fp)

ImageDimPlot(obj, group.by = c('insitutype_cluster_5_clusters'), 
             cols = brewer.pal(5, 'Paired'),
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("insitutype_cluster_5_clusters")


ImageDimPlot(obj, group.by = c('seurat_clusters_all_FOVs'), 
             #cols = brewer.pal(10, 'Paired'), 
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("Seurat_clusters_all_FOVs")




## Plot: Clusters proportion by TMA
dittoBarPlot(obj, "insitutype_cluster_5_clusters", group.by = "TMA") 
dittoBarPlot(obj, "seurat_clusters_all_FOVs", group.by = "TMA")


## Plot: Markers Dotplot 

library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'insitutype_cluster_5_clusters')
top_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  dplyr::top_n(wt = logFC, n = 5) %>%
  ungroup() %>%
  pull(feature) %>%
  unique()

insitutype_cluters_top_markers <- top_markers

Idents(obj) <- "insitutype_cluster_5_clusters"
DotPlot(obj, features = unique(c(insitutype_cluters_top_markers))) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Top Markers from Insitutype Clusters")


library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'seurat_clusters_all_FOVs')
top_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  dplyr::top_n(wt = logFC, n = 5) %>%
  ungroup() %>%
  pull(feature) %>%
  unique()
seurat_cluters_top_markers <- top_markers

Idents(obj) <- "seurat_clusters_all_FOVs"
DotPlot(obj, features = unique(c(seurat_cluters_top_markers))) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Top Markers from Seurat Clusters")


## Plot: Markers Vlnplot
library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'insitutype_cluster_5_clusters')
insitutype_cluters_top_auc_markers <- markers[markers$auc > 0.70, ]$feature

VlnPlot(obj, features = insitutype_cluters_top_auc_markers, group.by = "TMA", pt.size = 0)

library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'seurat_clusters_all_FOVs')
seurat_cluters_top_auc_markers <- markers[markers$auc > 0.80, ]$feature

VlnPlot(obj, features = seurat_cluters_top_auc_markers, group.by = "TMA", pt.size = 0)


length(intersect(insitutype_cluters_top_auc_markers, 
                 seurat_cluters_top_auc_markers))
intersect(insitutype_cluters_top_auc_markers, 
          seurat_cluters_top_auc_markers)




### Neonate 2 ----------------

Neonate_2_dir <- "/mnt/scratch2/CosMX/Diego/DIegoRSVNeonate2831439_08_11_2024_13_26_50_655/flatFiles/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1"

## Create a list of Seurat objects by FOV/TMA - accordingly to Hannah ppt slides
Neonate_2_obj <- LoadNanostring(Neonate_2_dir, fov = 'FOV')

Neonate_2_tx <- fread(paste0(Neonate_2_dir, "/", "831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1_metadata_file.csv.gz"))
Neonate_2_tx <- as.data.frame(Neonate_2_tx)

Neonate_2_meta <- fread(paste0(Neonate_2_dir, "/", "831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1_metadata_file.csv.gz"))

Neonate_2_obj$fov <- Neonate_2_meta$fov
Neonate_2_obj$Area <- Neonate_2_meta$Area
Neonate_2_obj$Area.um2 <- Neonate_2_meta$Area.um2

# match cell id in meta obj
Neonate_2_meta$cell = paste0(as.character(Neonate_2_meta$cell_ID), "_", Neonate_2_meta$fov)

# remove control and system probes
which <- grep('Neg*|Syst*', rownames(Neonate_2_obj))
Neonate_2_obj <- Neonate_2_obj[-c(which), ]

# Remove low quality cells
Neonate_2_obj <- subset(Neonate_2_obj, nCount_Nanostring > 10 & nFeature_Nanostring > 10 & Area.um2 < 500)

#Update meta
Neonate_2_meta <- Neonate_2_meta[Neonate_2_meta$cell %in% colnames(Neonate_2_obj), ]
dim(Neonate_2_meta) #

# Neonate_2_obj_bk <- Neonate_2_obj


# Get the FOVs from Hannah's ppt slides later if needed
# Checking duplicates on TMA names

## Not working this way 
# TMA_names <- c(
#   'xxxxxxx',
#   (...)
# )
# 
# table(TMA_names)
# 
# # Classify different embryos by FOV numbers
Neonate_2_obj$TMA <- ""
Neonate_2_obj[[]] <- Neonate_2_obj[[]] %>%
  mutate(TMA = case_when(
    fov %in% 1:10 ~ 'Ap16-2653 A2 1',
    fov %in% 11:19 ~ 'Ap16-2653 A2 2',
    fov %in% 20:28 ~ 'Ap16-2653 A2 3',
    fov %in% 29:37 ~ 'Ap16-2656 A2 1',
    fov %in% 38:46 ~ 'Ap16-2655 A2 1',
    fov %in% 47:55 ~ 'Ap16-2655 A2 2',
    fov %in% 56:64 ~ 'Ap16-2654 A2 1',
    fov %in% 65:73 ~ 'Ap16-2657 A2 1',
    fov %in% 74:82 ~ 'Ap16-2657 A2 2',
    fov %in% 83:91 ~ 'Ap16-2656 A2 2',
    fov %in% 92:100 ~ 'Ap16-2656 A2 3',
    fov %in% 101:109 ~ 'Ap17-0325 A1 1',
    fov %in% 240:253 ~ 'Ap17-0325 A1 2',
    fov %in% 110:118 ~ 'Ap16-26587 A2 1',
    fov %in% 119:127 ~ 'Ap16-26587 A2 2',
    fov %in% 254:265 ~ 'Ap16-3768 A1 1',
    fov %in% 266:277 ~ 'Ap17-0327 A1 1',
    fov %in% 128:136 ~ 'Ap17-0327 A1 2',
    fov %in% 137:201 ~ 'Ap17-0326 A1 1',
    fov %in% 146:154 ~ 'Ap16-3770 A1 1',
    fov %in% 155:166 ~ 'Ap16-3769 A1 1',
    fov %in% 278:291 ~ 'Ap16-3769 A1 2',
    fov %in% 164:172 ~ 'Ap16-3772 A1 1',
    fov %in% 173:181 ~ 'Ap16-3771 A1 1',
    fov %in% 182:190 ~ 'Ap16-3771 A1 2',
    fov %in% 191:200 ~ 'Ap16-3768 A1 2',
    fov %in% 202:210 ~ 'Ap16-26587 A2 3',
    fov %in% 211:219 ~ 'Ap16-2657 A2 3',
    fov %in% 220:228 ~ 'Ap16-2656 A2 4',
    fov %in% 229:238 ~ 'Ap16-3770 A1 3',

    TRUE ~ TMA  # Retain existing value or specify a default if needed
  ))



# Create list of embryos per run (break down of ProcessObj())
obj_list <- list()
TMA_values <- names(table(Neonate_2_obj$TMA))
TMA_values <- TMA_values[-1]
for (i in TMA_values) {
  # Subset the Seurat object based on the current TMA value
  subset_seurat <- subset(Neonate_2_obj, TMA == i)

  # Add the subsetted Seurat object to the list
  obj_list[[i]] <- subset_seurat
}

# Set default boundary as segmentation for each object in the list
for (i in seq_along(obj_list)) {
  DefaultBoundary(obj_list[[i]][['FOV']]) <- 'segmentation'
}

Neonate_2_Sobj_list <- obj_list
print(Neonate_2_Sobj_list)

Neonate_2_Sobj_flat <- Neonate_2_obj


## Processing with Seurat Workflow and insitutype clustering
Out_list <- process_seurat_object(Neonate_2_Sobj_flat,
                                  Metadata_dir = "/mnt/scratch2/CosMX/Diego/DIegoRSVNeonate2831439_08_11_2024_13_26_50_655/flatFiles/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1_metadata_file.csv.gz",
                                  # Downsize_to_test = TRUE,
                                  # Downsize_n_cells = 100,
                                  n_insitutype_clusters = 10, 
                                  run_insitutype = TRUE,
                                  run_seuratworkflow = TRUE)

Neonate_2_Sobj_flat_proc <- Out_list$seurat_object
unsup <- Out_list$unsupervised_results


## Visualization 
obj <- Neonate_2_Sobj_flat_proc
names(obj@meta.data)[names(obj@meta.data) == "seurat_clusters"] <- "seurat_clusters_all_FOVs"

## Plot: nCount nFeature
# TMA_values <- names(table(obj$TMA))
ImageFeaturePlot(#subset(obj, TMA == TMA_values[1]), 
  obj, 
  border.color = 'white', 
  size = 1,
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  features = c('nCount_Nanostring'))

hist(obj$nCount_Nanostring)
median(obj$nCount_Nanostring)


ImageFeaturePlot(#subset(obj, TMA == TMA_values[1]), 
  obj, 
  border.color = 'white', 
  size = 1,
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  features = c('nFeature_Nanostring'))

hist(obj$nFeature_Nanostring)
median(obj$nFeature_Nanostring)


## Plot: Clusters over spatial plot
ImageDimPlot(obj, group.by = c('insitutype_cluster_10_clusters'), 
             cols = brewer.pal(10, 'Paired'),
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("insitutype_cluster_10_clusters")

cols <-
  c(
    '#8DD3C7',
    '#BEBADA',
    '#FB8072',
    '#80B1D3',
    '#FDB462',
    '#B3DE69',
    '#FCCDE5',
    '#D9D9D9',
    '#BC80BD',
    '#CCEBC5',
    '#FFED6F',
    '#E41A1C',
    '#377EB8',
    '#4DAF4A',
    '#984EA3',
    '#FF7F00',
    '#FFFF33',
    '#A65628',
    '#F781BF',
    '#999999'
  )
cols <- cols[seq_along(unique(unsup$clust))]
names(cols) <- unique(unsup$clust)
fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup, col = cols[unsup$clust])
#class(fp)
print(fp)

## Re-run only insitutype 
Out_list_2 <- process_seurat_object(Neonate_2_Sobj_flat_proc, #already have run_seuratworkflow ran on it
                                    Metadata_dir = "/mnt/scratch2/CosMX/Diego/DIegoRSVNeonate2831439_08_11_2024_13_26_50_655/flatFiles/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1_metadata_file.csv.gz",
                                    # Downsize_to_test = TRUE,
                                    # Downsize_n_cells = 100,
                                    n_insitutype_clusters = 5, 
                                    run_insitutype = TRUE,
                                    run_seuratworkflow = FALSE)


Neonate_2_Sobj_flat_proc_2 <- Out_list_2$seurat_object
unsup_2 <- Out_list_2$unsupervised_results

save(Neonate_2_Sobj_flat_proc_2, unsup_2, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Clustering_Neonate_2.rda")

obj <- Neonate_2_Sobj_flat_proc_2
names(obj@meta.data)[names(obj@meta.data) == "seurat_clusters"] <- "seurat_clusters_all_FOVs"

cols <-
  c(
    '#8DD3C7',
    '#BEBADA',
    '#FB8072',
    '#80B1D3',
    '#FDB462',
    '#B3DE69',
    '#FCCDE5',
    '#D9D9D9',
    '#BC80BD',
    '#CCEBC5',
    '#FFED6F',
    '#E41A1C',
    '#377EB8',
    '#4DAF4A',
    '#984EA3',
    '#FF7F00',
    '#FFFF33',
    '#A65628',
    '#F781BF',
    '#999999'
  )
cols <- cols[seq_along(unique(unsup_2$clust))]
names(cols) <- unique(unsup_2$clust)
fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup_2, col = cols[unsup_2$clust])
#class(fp)
print(fp)

ImageDimPlot(obj, group.by = c('insitutype_cluster_5_clusters'), 
             cols = brewer.pal(5, 'Paired'),
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("insitutype_cluster_5_clusters")


ImageDimPlot(obj, group.by = c('seurat_clusters_all_FOVs'), 
             #cols = brewer.pal(10, 'Paired'), 
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("Seurat_clusters_all_FOVs")




## Plot: Clusters proportion by TMA
dittoBarPlot(obj, "insitutype_cluster_5_clusters", group.by = "TMA") 
dittoBarPlot(obj, "seurat_clusters_all_FOVs", group.by = "TMA")


## Plot: Markers Dotplot 

library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'insitutype_cluster_5_clusters')
top_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  dplyr::top_n(wt = logFC, n = 5) %>%
  ungroup() %>%
  pull(feature) %>%
  unique()

insitutype_cluters_top_markers <- top_markers

Idents(obj) <- "insitutype_cluster_5_clusters"
DotPlot(obj, features = unique(c(insitutype_cluters_top_markers))) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Top Markers from Insitutype Clusters")


library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'seurat_clusters_all_FOVs')
top_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  dplyr::top_n(wt = logFC, n = 5) %>%
  ungroup() %>%
  pull(feature) %>%
  unique()
seurat_cluters_top_markers <- top_markers

Idents(obj) <- "seurat_clusters_all_FOVs"
DotPlot(obj, features = unique(c(seurat_cluters_top_markers))) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Top Markers from Seurat Clusters")


## Plot: Markers Vlnplot
library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'insitutype_cluster_5_clusters')
insitutype_cluters_top_auc_markers <- markers[markers$auc > 0.70, ]$feature

VlnPlot(obj, features = insitutype_cluters_top_auc_markers[], group.by = "TMA", pt.size = 0)

library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'seurat_clusters_all_FOVs')
seurat_cluters_top_auc_markers <- markers[markers$auc > 0.80, ]$feature

VlnPlot(obj, features = seurat_cluters_top_auc_markers[13:24], group.by = "TMA", pt.size = 0)


length(intersect(insitutype_cluters_top_auc_markers, 
                 seurat_cluters_top_auc_markers))
intersect(insitutype_cluters_top_auc_markers, 
          seurat_cluters_top_auc_markers)




### Using 10 clusters -------------
load("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Clustering_Adult_2.rda")

obj <- Adult_2_Sobj_flat_proc_2
names(obj@meta.data)[names(obj@meta.data) == "seurat_clusters"] <- "seurat_clusters_all_FOVs"


ImageDimPlot(obj, group.by = c('insitutype_cluster_10_clusters'), 
             cols = brewer.pal(10, 'Paired'),
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("insitutype_cluster_10_clusters")



## Plot: Clusters proportion by TMA
dittoBarPlot(obj, "insitutype_cluster_10_clusters", group.by = "TMA") 

## Plot: Markers Dotplot 

library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
top_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  dplyr::top_n(wt = logFC, n = 5) %>%
  ungroup() %>%
  pull(feature) %>%
  unique()

insitutype_cluters_top_markers <- top_markers
# Ordering levels for plotting
obj@meta.data$insitutype_cluster_10_clusters <- factor(obj@meta.data$insitutype_cluster_10_clusters, levels =c(names(table(obj@meta.data$insitutype_cluster_10_clusters)))) 

Idents(obj) <- "insitutype_cluster_10_clusters"
DotPlot(obj, features = unique(c(insitutype_cluters_top_markers))) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1)) + 
  ggtitle("Top Markers from Insitutype Clusters")


## DEG from GeoMx 
library(readr)
Adu1_Neo1_GeoMx_DEG_YL <- read_csv("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adu1_Neo1_GeoMx_DEG_YL.csv") %>% data.frame()

names(Adu1_Neo1_GeoMx_DEG_YL)[names(Adu1_Neo1_GeoMx_DEG_YL) == "...1"] <- "gene"

dim(Adu1_Neo1_GeoMx_DEG_YL)
head(Adu1_Neo1_GeoMx_DEG_YL)

deg_filt_up <- Adu1_Neo1_GeoMx_DEG_YL %>% filter(adj.P.Val <= 0.05 &
                                                   logFC > 0)
deg_filt_down <- Adu1_Neo1_GeoMx_DEG_YL %>% filter(adj.P.Val <= 0.05 &
                                                     logFC < 0)

library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'insitutype_cluster_10_clusters')
insitutype_cluters_top_auc_markers <- markers[markers$auc > 0.70, ]$feature

# Insitutype cl. markers 
insitutype_cluters_top_auc_markers %in% deg_filt_up$gene %>% table() # 7/33
insitutype_cluters_top_auc_markers %in% deg_filt_down$gene %>% table() # 22/33

Adu2_CosMx_Adu1GeoMx_markers <- intersect(insitutype_cluters_top_auc_markers, deg_filt_down$gene)

## Dotplot 
Idents(obj) <- "insitutype_cluster_10_clusters"
DotPlot(obj, features = unique(c(Adu2_CosMx_Adu1GeoMx_markers))) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1)) + 
  ggtitle("GeoMx DEG in CosMx (Adult2)")


## Spatial plot 
ImageDimPlot(subset(obj, insitutype_cluster_10_clusters %in% c( "e")), 
  #obj,
             group.by = c('insitutype_cluster_10_clusters'), 
             cols = brewer.pal(10, 'Paired'),
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("GeoMx DEG in CosMx - best clusters (Adult2)")



VlnPlot(obj, features = Adu2_CosMx_Adu1GeoMx_markers[1:10], 
        pt.size = 0,
        group.by = "TMA")

VlnPlot(obj, features = Adu2_CosMx_Adu1GeoMx_markers[11:20], 
        pt.size = 0, 
        group.by = "TMA")


ImageFeaturePlot(#subset(obj, TMA == TMA_values[1]), 
  obj,
  #subset(obj, insitutype_cluster_5_clusters %in% c("c")),
  border.color = 'white', 
  size = 1,
  # min.cutoff = 'q10',
  # max.cutoff = 'q90',
  features = c('Cd74'))



## Saving markers to share with Diego
Adult_2_CosMx_Cluster_Markers <- markers[markers$auc > 0.70, ]
Adult_2_CosMx_Cluster_Markers$inGeoMx <- "No"
Adult_2_CosMx_Cluster_Markers[Adult_2_CosMx_Cluster_Markers$feature %in% Adu2_CosMx_Adu1GeoMx_markers, ]$inGeoMx <- "Yes"

# write.csv(Adult_2_CosMx_Cluster_Markers, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult_2_CosMx_Cluster_Markers.csv")




        # # #   GeoMx DEG -----> CosMx  # # #
        # # #   GeoMx DEG -----> CosMx  # # #
        # # #   GeoMx DEG -----> CosMx  # # #

# # #   Find DEG genes from GeoMx into CosMx data # # #

library(readr)
Adu1_Neo1_GeoMx_DEG_YL <- read_csv("Adu1_Neo1_GeoMx_DEG_YL.csv") %>% data.frame()

names(Adu1_Neo1_GeoMx_DEG_YL)[names(Adu1_Neo1_GeoMx_DEG_YL) == "...1"] <- "gene"

dim(Adu1_Neo1_GeoMx_DEG_YL)
head(Adu1_Neo1_GeoMx_DEG_YL)

deg_filt_up <- Adu1_Neo1_GeoMx_DEG_YL %>% filter(adj.P.Val <= 0.05 &
                                    logFC > 0)
deg_filt_down <- Adu1_Neo1_GeoMx_DEG_YL %>% filter(adj.P.Val <= 0.05 &
                                                   logFC < 0)



## GeoMx_DEG %in% CosMx_Adu2 ---------
load("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Clustering_Adult_2.rda")
obj <- Adult_2_Sobj_flat_proc_2
names(obj@meta.data)[names(obj@meta.data) == "seurat_clusters"] <- "seurat_clusters_all_FOVs"
library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'insitutype_cluster_5_clusters')
insitutype_cluters_top_auc_markers_Adult2 <- markers[markers$auc > 0.70, ]$feature


# All 1k genes
rownames(obj) %in% deg_filt_up$gene %>% table() # 312/997 intersected genes (from up)
rownames(obj) %in% deg_filt_down$gene %>% table() # 340/997 intersected genes (from down)

# Insitutype cl. markers 
insitutype_cluters_top_auc_markers_Adult2 %in% deg_filt_up$gene %>% table() # 0/7
insitutype_cluters_top_auc_markers_Adult2 %in% deg_filt_down$gene %>% table() # 6/7

Adu2_CosMx_Neo1GeoMx_markers <- intersect(insitutype_cluters_top_auc_markers_Adult2, deg_filt_down$gene)


# Visualizations 
## Dotplot 
Idents(obj) <- "insitutype_cluster_5_clusters"
DotPlot(obj, features = unique(c(Adu2_CosMx_Neo1GeoMx_markers))) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("GeoMx DEG in CosMx (Adult2)")


## Spatial plot 
ImageDimPlot(subset(obj, insitutype_cluster_5_clusters %in% c( "c", "d")), 
             group.by = c('insitutype_cluster_5_clusters'), 
             cols = brewer.pal(5, 'Paired'),
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("GeoMx DEG in CosMx - best clusters (Adult2)")

ImageDimPlot(subset(obj, insitutype_cluster_5_clusters %in% c( "c")),
             group.by = c('insitutype_cluster_5_clusters'), 
             cols = brewer.pal(5, 'Paired'),
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("GeoMx DEG in CosMx - best clusters (Adult2)")

ImageDimPlot(subset(obj, insitutype_cluster_5_clusters %in% c( "d")), 
             group.by = c('insitutype_cluster_5_clusters'), 
             cols = brewer.pal(5, 'Paired'),
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("GeoMx DEG in CosMx - best clusters (Adult2)")


ImageFeaturePlot(#subset(obj, TMA == TMA_values[1]), 
  obj,
  #subset(obj, insitutype_cluster_5_clusters %in% c("c")),
  border.color = 'white', 
  size = 1,
  # min.cutoff = 'q10',
  # max.cutoff = 'q90',
  features = c('Gsn'))

ImageFeaturePlot(#subset(obj, TMA == TMA_values[1]), 
  obj,
  #subset(obj, insitutype_cluster_5_clusters %in% c("c")),
  border.color = 'white', 
  size = 1,
  # min.cutoff = 'q10',
  # max.cutoff = 'q90',
  features = c('Col1a1'))

ImageFeaturePlot(#subset(obj, TMA == TMA_values[1]), 
  obj,
  #subset(obj, insitutype_cluster_5_clusters %in% c("c")),
  border.color = 'white', 
  size = 1,
  # min.cutoff = 'q10',
  # max.cutoff = 'q90',
  features = c('Bgn'))




## GeoMx_DEG %in% CosMx_Neo2 ----------
load("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Clustering_Neonate_2.rda")
obj <- Neonate_2_Sobj_flat_proc_2
names(obj@meta.data)[names(obj@meta.data) == "seurat_clusters"] <- "seurat_clusters_all_FOVs"
library(presto)
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, 'insitutype_cluster_5_clusters')
insitutype_cluters_top_auc_markers_Neonate2 <- markers[markers$auc > 0.70, ]$feature

# All 1k genes
rownames(obj) %in% deg_filt_up$gene %>% table() # 312/997 intersected genes (from up)
rownames(obj) %in% deg_filt_down$gene %>% table() # 340/997 intersected genes (from down)

# Insitutype cl. markers 
insitutype_cluters_top_auc_markers_Neonate2 %in% deg_filt_up$gene %>% table() # 0/25
insitutype_cluters_top_auc_markers_Neonate2 %in% deg_filt_down$gene %>% table() # 1/4

Neo2_CosMx_Neo1GeoMx_markers <- intersect(insitutype_cluters_top_auc_markers_Neonate2, deg_filt_down$gene)


# Visualizations 

## Dotplot 
Idents(obj) <- "insitutype_cluster_5_clusters"
DotPlot(obj, features = unique(c(Neo2_CosMx_Neo1GeoMx_markers))) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("GeoMx DEG in CosMx (Neonate2)")


## Spatial plot 
ImageDimPlot(subset(obj, insitutype_cluster_5_clusters %in% c("e")), 
             group.by = c('insitutype_cluster_5_clusters'), 
             cols = brewer.pal(5, 'Paired'),
             axes = TRUE) +
  # coord_flip() + 
  ggtitle("GeoMx DEG in CosMx - best clusters (Neonate2)")











# ## GeoMx gene exp - from Yutian --------
# 
# ## Prepare exp mtx and metadata 
# exp_mtx <- data.table::fread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/ex_df.csv") %>% data.frame()
# 
# head(exp_mtx)
# exp_mtx$V1
# rownames(exp_mtx) <- exp_mtx$V1
# exp_mtx$V1 <- NULL
# 
# meta <- readxl::read_excel("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/annotation_Diego_Geomx_rsv.xlsx") %>% data.frame()
# head(meta)
# 
# meta$Sample_ID %>% length() #94
# colnames(exp_mtx) %>% length() #92
# 
# meta$Sample_ID <- paste0(meta$Sample_ID, ".dcc")
# 
# setdiff(meta$Sample_ID, colnames(exp_mtx))
# setdiff(colnames(exp_mtx), meta$Sample_ID)
# 
# meta$Sample_ID <- gsub("-", ".", meta$Sample_ID)
# 
# setdiff(meta$Sample_ID, colnames(exp_mtx))
# setdiff(colnames(exp_mtx), meta$Sample_ID)
# 
# length(intersect(meta$Sample_ID, colnames(exp_mtx)))
# 
# 
# ## Check metadata columns
# head(meta)
# table(meta$Slide.Name)
# meta$Group <- NA
# meta[meta$Slide.Name %in% "Adult.1_8.30.24", ]$Group <- "Adult_1"
# meta[meta$Slide.Name %in% "Adult.2_8.30.24", ]$Group <- "Adult_2"
# meta[meta$Slide.Name %in% "Neonate.1_8.30.24", ]$Group <- "Neonate_1"
# meta[meta$Slide.Name %in% "Neonate.2_8.30.24", ]$Group <- "Neonate_2"
# meta[meta$Slide.Name %in% "No Template Control", ]$Group <- "Control"
# 
# 
# 
# ## Yutian said that in GeoMx, Adult_1 and Neonate_1 are the best data. 
# ## Then I should get the DEG from Adult_1 vs Neonate_1 then try to find those genes in CosMx (Adult_2 and Neonate_2)
# 
# Ad1_Ne1_meta <- meta[meta$Group %in% c("Neonate_1", "Adult_1"), ]
# Ad1_Ne1_exp <- exp_mtx[, colnames(exp_mtx) %in% Ad1_Ne1$Sample_ID]
# 
# 
# 
# ## Run Differential Expression Analysis on GeoMx -------------
# library(standR)
# library(SpatialExperiment)
# library(limma)
# library(ExperimentHub)
# library(dplyr)
# library(GeomxTools)
# library(SpatialDecon)
# library(Seurat)
# library(lme4)
# library(parallel)
# library(EnhancedVolcano)
# library(stringr)
# 
# 
# 
# 
# ### 1. DEG directly from the exp mtx and metada Yutian handled to me
# 
# ## Filter low expressed genes 
# # cout == 5 in more than 90% of ROIs 
# # Doing it separately for Adult_1 and Neonate_1
# 
# ## logCPM it
# library(edgeR)
# # Create DGEList object
# y <- DGEList(counts=Ad1_Ne1_exp)
# # Normalize library sizes
# y <- calcNormFactors(y)
# # Calculate logCPM
# Ad1_Ne1_exp_logCPM <- cpm(y, log=TRUE, prior.count=1)
# 
# # Perform PCA - counts
# pca1 <- prcomp(t(Ad1_Ne1_exp)) 
# aux <- as.data.frame(pca1$x[, 1:3]) 
# # Merge data
# scores1 <- merge(Ad1_Ne1_meta, aux, by.y=0, by.x="Sample_ID", all.x=T)
# # Load necessary library and set theme
# library(ggplot2); theme_set(theme_classic())
# # Create PCA plot for all samples separated by articles
# ggplot(scores1, aes(x=PC1, y=PC2, colour=factor(Group), )) +
#   geom_point(size = 4) +
#   scale_color_manual(values=c("gold", "darkgreen"), name="Group") +
#   xlab(paste0("PC1 (", prettyNum(summary(pca1)$importance[2,1]*100, digits = 2), "%)")) +
#   ylab(paste0("PC2 (", prettyNum(summary(pca1)$importance[2,2]*100, digits = 2), "%)")) +
#   scale_x_continuous(labels = scales::scientific_format()) +
#   scale_y_continuous(labels = scales::scientific_format()) +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14, face = "bold"),
#         legend.title = element_text(size = 16),  
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, color = "black", face = "bold")) +
#   ggtitle("PCA counts") 
# 
# 
# # Perform PCA - logCPM 
# pca1 <- prcomp(t(Ad1_Ne1_exp_logCPM)) 
# aux <- as.data.frame(pca1$x[, 1:3]) 
# # Merge data
# scores1 <- merge(Ad1_Ne1_meta, aux, by.y=0, by.x="Sample_ID", all.x=T)
# # Load necessary library and set theme
# library(ggplot2); theme_set(theme_classic())
# # Create PCA plot for all samples separated by articles
# ggplot(scores1, aes(x=PC1, y=PC2, colour=factor(Roi), )) +
#   geom_point(size = 4) +
#   #scale_color_manual(values=c("gold", "darkgreen"), name="Group") +
#   xlab(paste0("PC1 (", prettyNum(summary(pca1)$importance[2,1]*100, digits = 2), "%)")) +
#   ylab(paste0("PC2 (", prettyNum(summary(pca1)$importance[2,2]*100, digits = 2), "%)")) +
#   scale_x_continuous(labels = scales::scientific_format()) +
#   scale_y_continuous(labels = scales::scientific_format()) +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14, face = "bold"),
#         legend.title = element_text(size = 16),  
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16, color = "black", face = "bold")) +
#   ggtitle("PCA logCPM") 
# 
# 
# # It’s missing the patient ID in the metadata
# # So I don’t know if each ROI comes from a different patient or what 
# 
# # Yutian now told me she should have de DEG list. Waiting on her to share it with me ...
# 
# 
# 
# 
# 
# 
# 
# 
# Ne1_smp_ID <- Ad1_Ne1_meta[Ad1_Ne1_meta$Group %in% "Neonate_1", ]$Sample_ID
# Ne1_logCPM <- Ad1_Ne1_exp_logCPM[, colnames(Ad1_Ne1_exp_logCPM) %in% Ne1_smp_ID]
# 
# 
# 
# 
# ### 2. DEG when we have a "NanoStringGeoMxSet" - i don't have it cause i dont have the "bundle files" out of the machine (DCCFiles, PKCFiles, SampleAnnotationFile to read through "readNanoStringGeoMxSet()")
# 
# load("/mnt/scratch1/maycon/Ayo_GeoMx/Previous_analysis/Rds/20240319_114_sDAS_Miami_Omotoso_WTA.RData")
# print(target_data) # input data has to be like that 
# # learn how to do it or get a new DEG tutorial from Nanostring tutorials 
# 
# 
# 
# Get_RNA_de <- function(target_data = NULL, 
#                        Country_of_origin_chr = NULL, 
#                        Diagnosis_chr = NULL, 
#                        segment_chr = NULL,
#                        nCores_num = NULL) {
#   
#   contrast_samples <- rownames(filter(pData(target_data),
#                                       
#                                       Country.of.origin %in% Country_of_origin_chr,
#                                       Diagnosis %in% Diagnosis_chr,
#                                       segment %in% segment_chr))
#   contrast <- NULL
#   contrast <- target_data[,contrast_samples] 
#   model_formula <- ~ Country.of.origin + (1|PatientID) #variable of interest + "random effect" variable. It should "regress out" every variation due to patients  
#   group_variable <- "Country.of.origin"
#   
#   # parallel::detectCores()-1 # 191 cores 
#   start.time <- Sys.time()
#   de_contrast <- GeomxTools::mixedModelDE(contrast,
#                                           elt = "log_q",
#                                           modelFormula = model_formula, # Model for contrast in Index
#                                           groupVar = group_variable, # Group variable for contrast is Index
#                                           nCores = nCores_num,
#                                           multiCore = FALSE)
#   
#   de_contrast <- formatLMMResults(de_contrast, p_adjust_method = "none")
#   
#   end.time <- Sys.time() 
#   print(end.time - start.time) # Time difference of 4.42347 mins
#   return(de_contrast)
# }
# 
# 
# 
# 
# Plot_Volcano <- function(DEGs = NULL, 
#                          Group_1 = NULL,
#                          Group_2 = NULL,
#                          Segment = NULL,
#                          ylim = NULL) {
#   
#   
#   # Read DEGs
#   names(DEGs)[names(DEGs) == "Estimate"] <- "logFC"
#   names(DEGs)[names(DEGs) == "Feature"] <- "feature"
#   names(DEGs)[names(DEGs) == "P"] <- "pval"
#   names(DEGs)[names(DEGs) == "NONE"] <- "padj"
#   
#   # Filter by sig. DEGs 
#   DEGs_sig <- DEGs[DEGs$pval <= 0.05 & (DEGs$logFC >= 0.5 | DEGs$logFC <= -0.5), ]
#   
#   # Define top and bottom DEGs
#   DEGs_sig[order(DEGs_sig$logFC,  decreasing = TRUE), c("feature", "Comparison", "logFC")][1:20, c(1)] -> top20;       print(paste0("Top Up 20 DEGs_sig: ", top20))
#   print("-----")
#   DEGs_sig[order(DEGs_sig$logFC,  decreasing = FALSE), c("feature", "Comparison", "logFC")][1:20, c(1)] -> bottom_20;
#   print(paste0("Top Down 20 DEGs_sig: ", bottom_20))
#   
#   Total_DEGs_sig <- dim(DEGs_sig[DEGs_sig$logFC >= 0.5 | DEGs_sig$logFC <= -0.5, ])[1]
#   Up_DEGs_sig <- dim(DEGs_sig[DEGs_sig$logFC >= 0.5, ])[1]
#   Down_DEGs_sig <- dim(DEGs_sig[DEGs_sig$logFC <= -0.5, ])[1]
#   
#   # Plot volcano plot 
#   library(EnhancedVolcano)
#   rownames(DEGs) <- DEGs$feature 
#   p <- EnhancedVolcano(
#     DEGs,
#     lab = DEGs$feature,
#     x = 'logFC',
#     y = 'padj',
#     ylim = c(0, ylim),
#     selectLab = c(top20[1:6], bottom_20[1:6]),
#     pCutoff = 0.05, #pvalue cutoff line
#     FCcutoff = 0.5, #foldChange cutoff line
#     xlab = paste0("<---- ",Group_1,"   " ,"Log2 FoldChange", "   ",Group_2," ---->" ),
#     pointSize = 4.0,
#     labSize = 6.0,
#     labCol = 'black',
#     labFace = 'bold',
#     boxedLabels = TRUE,
#     #colAlpha = 4/5,
#     legendPosition = 'right',
#     legendLabSize = 14,
#     legendIconSize = 4.0,
#     drawConnectors = TRUE,
#     widthConnectors = 1.0,
#     colConnectors = 'black',
#     max.overlaps = Inf,
#     caption = paste0("Total DEGs: ", Total_DEGs_sig, "\n",
#                      "Up genes: ", Up_DEGs_sig, " towards ", Group_2, "\n",
#                      "Down genes: ", Down_DEGs_sig, " towards ", Group_1, "\n",
#                      "p <= 0.05, FC >= |0.5|"),
#     title = paste0("DEGs ", Group_1, " vs ", Group_2, " [", Segment, "]"),
#     subtitle = "" )
#   
#   return(list(p, DEGs_sig))
# }








