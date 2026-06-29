library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(shiny)
library(plotly)
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(qs)

## Tutorials 
# insitutype with cell type ref: file:///Users/mmarcao/Downloads/rstudio-export%20(41)/NSCLC-semi-supervised-cell-typing-vignette.html


## Cell type reference - profile mtx from CosMx 
# at https://github.com/Nanostring-Biostats/CellProfileLibrary/blob/master/Mouse/Adult/Lung_MCA.RData

load("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/Lung_MCA.RData")
metadata
cellGroups
profile_matrix[1:10, 1:10] # row: mouse genes; col: cell type; content: cell type "fraction" by gene
dim(profile_matrix)

rownames(profile_matrix)


## Load seurat obj 
Sobj <- qs::qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Objects/all_4_slides_processed.qs")
rownames(Sobj)

 
intersect(rownames(profile_matrix), rownames(Sobj)) # matching genes (863)
setdiff(rownames(Sobj), rownames(profile_matrix)) # genes not matched from CosMx 1k panel (144)

# profile_matrix <- profile_matrix[rownames(profile_matrix) %in% rownames(Sobj), ]
# dim(profile_matrix)
# 
# Sobj_sub <- subset(Sobj, features = rownames(profile_matrix))
# 
# Sobj_100 <- subset(Sobj_sub, cells = colnames(Sobj_sub)[1:100])


print(head(profile_matrix))
# Monocyte.progenitor Neutrophil    NK.cell     Nuocyte
# Aatk            0.0000000 0.00000000 0.00000000 0.000000000
# Abl1            0.0000000 0.02857143 0.01090909 0.006622517
# Abl2            0.0000000 0.00000000 0.01090909 0.013245033
# Acacb           0.0000000 0.00000000 0.00000000 0.000000000
# Ace             0.1538462 0.02857143 0.02909091 0.026490066
# Ackr3           0.0000000 0.00000000 0.00000000 0.006622517

# ## Run insitutype-unsupervised 
# unsup <- insitutype(
#   x = as.matrix(t(obj@assays$Nanostring$counts)),
#   neg = obj@meta.data$nCount_negprobes,
#   reference_profiles = NULL,
#   bg = NULL,
#   #n_clusts = 10,
#   n_clusts = 15,
#   n_phase1 = 200,
#   n_phase2 = 500,
#   n_phase3 = 2000,
#   n_starts = 1,
#   max_iters = 5)
# 
# # unsup_cl10 <- unsup
# unsup_cl15 <- unsup
# 
# obj$insitutype_cluster <- 'other'
# obj$insitutype_cluster[names(unsup$clust)] <- unsup$clust
# # names(obj@meta.data)[names(obj@meta.data) == "insitutype_cluster"] <- paste0("insitutype_cluster_", "10", "_clusters")
# names(obj@meta.data)[names(obj@meta.data) == "insitutype_cluster"] <- paste0("insitutype_cluster_", "15", "_clusters")


# ## Run insitutype-celltype_ref (semi-supervised)
# semisup <- insitutype(
#   x = t(as.matrix(Sobj@assays$Nanostring$counts)),
#   neg = Sobj@meta.data$nCount_negprobes,
#   # cohort = cohort,
#   # Enter your own per-cell background estimates here if you
#   # have them; otherwise insitutype will use the negprobes to
#   # estimate background for you.
#   bg = NULL,
#   # condensed to save time. n_clusts = 5:15 would be more optimal
#   n_clusts = c(13),
#   reference_profiles = as.matrix(profile_matrix),
#   update_reference_profiles = FALSE,
#   # choosing inadvisably low numbers to speed the vignette; using the defaults
#   # in recommended.
#   n_phase1 = 200,
#   n_phase2 = 500,
#   n_phase3 = 2000,
#   n_starts = 1,
#   max_iters = 5
# ) 
# 
# Sobj$insitutype_cluster <- 'other'
# Sobj$insitutype_cluster[names(semisup$clust)] <- semisup$clust
# names(Sobj@meta.data)[names(Sobj@meta.data) == "insitutype_cluster"] <- paste0("insitutype_cluster_", "SemiSuperv","_clusters")


## Run (supervised)
sup <- insitutypeML(x = t(as.matrix(Sobj@assays$Nanostring$counts)),
                    neg = Sobj@meta.data$nCount_negprobes,
                    #cohort = cohort,
                    reference_profiles = as.matrix(profile_matrix)) 


Sobj$insitutype_cluster <- 'other'
Sobj$insitutype_cluster[names(sup$clust)] <- sup$clust
names(Sobj@meta.data)[names(Sobj@meta.data) == "insitutype_cluster"] <- paste0("insitutype_cluster_", "Superv","_celltypes")




## Visualizations 


library(dittoSeq)
dittoBarPlot(Sobj, "insitutype_cluster_Superv_celltypes", 
             group.by = "TMA_2",
             split.by = "Slide",
             #color.panel = cols, 
             main = "Clusters by TMAs (ROIs)")


## Renaming cell types to major cell types 
# Alveolar bipotent progenitor
Sobj@meta.data$major_celltype <- "Other"
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Alveolar.bipotent.progenitor", ]$major_celltype <- "Alveolar_bipotent_progenitor"

# Alveolar macrophage
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Alveolar.macrophage.Ear2.high", "Alveolar.macrophage.Pclaf.high"), ]$major_celltype <- "Alveolar_macrophage"

# AT1 and AT2 cells
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "AT1.cell", ]$major_celltype <- "AT1_cell"
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "AT2.cell", ]$major_celltype <- "AT2_cell"

# Lymphocytes
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("B.cell", "IgA.producing.B.cell", "NK.cell", "T.cell.Cd8b1.high"), ]$major_celltype <- "Lymphocyte"

# Granulocytes
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Basophil", "Eosinophil", "Neutrophil"), ]$major_celltype <- "Granulocyte"

# Epithelial cells
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Ciliated.cell", "Clara.cell"), ]$major_celltype <- "Epithelial_cell"

# Dendritic cells
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Dendritic.cell.Gngt2.high", "Dendritic.cell.H2.M2.high", "Dendritic.cell.Mgl2.high", "Dendritic.cell.Naaa.high", "Dendritic.cell.Tubb5.high", "Plasmacytoid.dendritic.cell"), ]$major_celltype <- "Dendritic_cell"

# Endothelial cells
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Endothelial.cell.Kdr.high", "Endothelial.cell.Tmem100.high", "Endothelial.cell.Vwf.high"), ]$major_celltype <- "Endothelial_cell"

# Interstitial macrophage
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Interstitial.macrophage", ]$major_celltype <- "Interstitial_macrophage"

# Monocyte progenitor
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Monocyte.progenitor", ]$major_celltype <- "Monocyte_progenitor"

# Stromal cells
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Stromal.cell.Acta2.high", "Stromal.cell.Dcn.high", "Stromal.cell.Inmt.high"), ]$major_celltype <- "Stromal_cell"

# Nuocyte
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Nuocyte", ]$major_celltype <- "Nuocyte"


## Set color map 
cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- Sobj$major_celltype %>% table() %>% names()
cols <- setNames(cols[1:length(factor_levels)], factor_levels)



library(dittoSeq)
dittoBarPlot(Sobj, "major_celltype", 
             group.by = "TMA_2",
             split.by = "Slide",
             color.panel = cols, 
             main = "Clusters by TMAs (ROIs)")

celltype_freq <- data.frame(table(Sobj@meta.data$major_celltype))
names(celltype_freq) <- c("Cell_type", "N_of_cells")
library(ggplot2); theme_set(theme_classic())
ggplot(celltype_freq, aes(reorder(x = factor(Cell_type), N_of_cells), y = N_of_cells, fill = cols)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_identity() +
  labs(#title = "Frequency of Cells by Sample ID",
    x = "Cell_type",
    y = "N_of_cells") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16),  
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, color = "black", face = "bold")) +
  labs(fill = "")
ggtitle("Frequency of Cells by Cell type")



## Plotting celltype barplots for each Comparison 
library(dittoSeq)
library(ggplot2)
library(Seurat)

# Define group comparisons
comparisons <- list(
  list(
    group_1 = "Neonate IFN and RSV infected", 
    group_2 = "Adult RSV infected"
  ),
  list(
    group_1 = "Neonate IFN and RSV reinfected", 
    group_2 = "Adult reinfected"
  ),
  list(
    group_1 = "Neonate RSV infected (NO IFN)", 
    group_2 = "Neonate IFN and RSV infected"
  )
)

# Create plots for each comparison
for (comp in comparisons) {
  # Subset Seurat object
  Sobj_subset <- subset(Sobj, subset = Group %in% c(comp$group_1, comp$group_2))
  
  ## dittoSeq Plot ##
p1 <-  dittoBarPlot(
    Sobj_subset, 
    var = "major_celltype",
    group.by = "Group",  # Compare between groups
    #split.by = "Slide",
    color.panel = cols,
    main = paste("Cell Types:", comp$group_1, "vs", comp$group_2)
  )

print(p1)
  
  ## ggplot2 Analysis ##
  # Calculate frequencies
  celltype_freq <- data.frame(table(
    Sobj_subset@meta.data$major_celltype
  ))
  
  colnames(celltype_freq) <- c("Cell_type", "N_of_cells")
  
  # Create comparative plot
  library(ggplot2); theme_set(theme_classic())
p2 <-  ggplot(celltype_freq, aes(reorder(x = factor(Cell_type), N_of_cells), y = N_of_cells, fill = cols)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_identity() +
    labs(#title = "Frequency of Cells by Sample ID",
      x = "Cell_type",
      y = "N_of_cells") +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 16),  
          legend.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 16, color = "black", face = "bold")) +
    labs(fill = "") +
  ggtitle(paste("Cell Types:", comp$group_1, "vs", comp$group_2))
  
  print(p2)
}


## The dittoseq barplots are already really informative about the differences in cell types across Groups 

## To be added to this block of analysis (do it in less then 1h)
## 1. check cell markers for each cell type (upseplot, dotplot top markers, heatmap with top markers and insitutype clusters)
## 2. describe the difference in cell type proportion found within each comparison
## 3. keep each dittoseq barplot for each Group along with each "N of cells" plots for each comparison



## Checking cell markers --------

# Get celltype markers 
library(presto)
Sobj@assays$RNA <- Sobj@assays$Nanostring
markers <- wilcoxauc(Sobj, 'major_celltype')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 10) %>% 
  ungroup() -> top_markers

# Check for overlapped celltype markers 
cluster_gene_list <- list()
list_names <- table(top_markers$group) %>% names()
for(i in 1:length(list_names)) {
  group <- top_markers[top_markers$group %in% list_names[i], ]$feature %>% list()
  cluster_gene_list[i] <- group
  names(cluster_gene_list)[i] <- list_names[i]
  
}
# Plot upsetplot 
library(UpSetR)
upset(fromList(cluster_gene_list), #order.by = "freq", 
      nsets = length(cluster_gene_list)) 

# Plot dotplot 
Idents(Sobj) <- "major_celltype"
DotPlot(Sobj, 
        features =  unique(c(top_markers$feature))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16)) + ggtitle("Top Cell Type Markers")


# Plot heatmap
library(dittoSeq)
Idents(Sobj) <- "major_celltype"
SeuObj_downsmp <- subset(x = Sobj, downsample = 500)
SeuObj_downsmp <- subset(SeuObj_downsmp, insitutype_cluster_15_clusters %in% c("e", "f"), invert = TRUE)
dittoHeatmap(SeuObj_downsmp, 
             unique(top_markers$feature),
             annot.by = c("insitutype_cluster_15_clusters", "major_celltype"),
             order.by = "major_celltype", #to make it work, complex and use_raster shoudl set to FALSE
             complex = F,
             use_raster = F, main = "Top Celltype Markers - downsampled data") #+ theme(legend.position = "none")

# qsave(Sobj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj.qs")







#

#

#

#

#

#

#

#

#

#

#

#




### Re-raning cell type and niche analysis - now keeping FOV for the whole data ### ----------------
### Apr 14, 2025 

# ### Keeping objects in a list - sure about it?? we want a niche analys across them all. So we need to find a way to keep the spatial info while using Reduce() to merge all this objects into one Seurat obj --------------
# Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj.qs")
# 
# # Set color map 
# cols <-
#   c(
#     '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
#     '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
#     '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
#   )
# 
# factor_levels <- Sobj$major_celltype %>% table() %>% names()
# cols <- setNames(cols[1:length(factor_levels)], factor_levels)
# 
# # Define group comparisons
# comparisons <- list(
#   list(
#     group_1 = "Neonate IFN and RSV infected", 
#     group_2 = "Adult RSV infected"
#   ),
#   list(
#     group_1 = "Neonate IFN and RSV reinfected", 
#     group_2 = "Adult reinfected"
#   ),
#   list(
#     group_1 = "Neonate RSV infected (NO IFN)", 
#     group_2 = "Neonate IFN and RSV infected"
#   )
# )
# 
# 
# 
# sobj_list <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult12_Neonate12_seuobj_list.rds")
# 
# # Define metadata directories
# metadata_dirs <- list(
#   Adult_1 = "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult1831439_08_11_2024_9_49_46_243/flatFiles/831439_IND_CR_FFPE_TMA_Adult1_STJ_N_R1",
#   Neonate_1 = "/mnt/scratch2/CosMX/Diego/DiegoRSVNeonate1831439_08_11_2024_11_18_07_116/flatFiles/831439_IND_CR_FFPE_TMA_Neonate1_STJ_N_R1",
#   Adult_2 = "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult2831439_08_11_2024_12_58_14_961/flatFiles/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1",
#   Neonate_2 = "/mnt/scratch2/CosMX/Diego/DIegoRSVNeonate2831439_08_11_2024_13_26_50_655/flatFiles/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1"
# )
# 
# Sobj@meta.data$cell
# metadata_file <- list.files(path = metadata_dirs[['Adult_1']], pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)
# metadata <- fread(metadata_file) %>% as.data.frame()
# head(metadata)
# 
# 
# # Adding flat files metadata, clusters, and cell type back to the listed seurat objects
# for (slide_name in names(metadata_dirs)) {
#   # Load metadata for the current slide
#   metadata_file <- list.files(path = metadata_dirs[[slide_name]], pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)
#   metadata <- fread(metadata_file) %>% as.data.frame()
#   
#   # Add Slide and cell columns to the metadata
#   metadata$Slide <- slide_name
#   metadata$cell <- paste0(as.character(metadata$cell_ID), "_", metadata$fov, "_", metadata$Slide)
#   
#   # Set rownames of the metadata to match cell identifiers
#   rownames(metadata) <- metadata$cell
#   
#   # Filter to seurat obj
#   seurat_obj <- sobj_list[[slide_name]]
#   seurat_obj@meta.data$Slide <- slide_name
#   seurat_obj@meta.data$cell <- paste0(as.character(seurat_obj@meta.data$cell_ID), #"_", seurat_obj@meta.data$fov,
#                                       "_", seurat_obj@meta.data$Slide)
#   #rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$cell
#   
#   # Filtering by flat-files metadata
#   metadata <- metadata[metadata$cell %in% seurat_obj@meta.data$cell, ]
#   # Filtering by processed files metadata
#   metadata <- merge(metadata, Sobj@meta.data, by="cell")
#   # Now actually filtering the metadata from the list based on the metadata to be merged 
#   outdated_meta <- seurat_obj@meta.data
#   metadata <- metadata[rownames(outdated_meta), ]
#   outdated_meta <- outdated_meta[rownames(metadata), ]
#   
#   if (!identical(outdated_meta$cell, metadata$cell)) {
#     stop("The 'cell' columns in 'outdated_meta' and 'metadata' are not identical. Stopping the script.")
#   }
#   
#   # Merge the metadata into the corresponding Seurat object
#   sobj_list[[slide_name]] <- AddMetaData(object = seurat_obj,
#                                          metadata = metadata)
#   
#   
# }


# Sobj <- BuildNicheAssay(object = Sobj, fov = "FOV", group.by = "major_celltype", niches.k = 10, neighbors.k = 20, cluster.name = "niches_10")


Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj.qs")

library(dittoSeq)
dittoBarPlot(Sobj, 
             "niches_5", 
             group.by = "Group", 
             #color.panel = cols,
             main = "")

library(dittoSeq)
dittoBarPlot(Sobj, 
             "niches_10", 
             group.by = "Group", 
             #color.panel = cols,
             main = "")

## Why we are getting niches only for a part of the data ??? 

# Checking coordinates (fine)
coord <- GetTissueCoordinates(Sobj)
table(is.na(coord$x))
table(is.na(coord$y))
table(is.na(coord$cell))


# Checking 
table(Sobj$niches_5) %>% sum()
table(Sobj$niches_5, useNA = "always") #403391 NAs ... 78% of cells getting niche NAs 

Sobj@images$FOV$centroids # only 112945 cells 
Sobj@images$FOV.2$centroids # only 201988 cells 

# The "Reduce()" merging problem is not working to keep FOVs
# Let's find a way to merge all the objects withou loosing the FOVs 
# Load the Seurat object list

## Adding flat files metadata info 
sobj_list <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult12_Neonate12_seuobj_list.rds")

# Define metadata directories
metadata_dirs <- list(
  Adult_1 = "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult1831439_08_11_2024_9_49_46_243/flatFiles/831439_IND_CR_FFPE_TMA_Adult1_STJ_N_R1",
  Neonate_1 = "/mnt/scratch2/CosMX/Diego/DiegoRSVNeonate1831439_08_11_2024_11_18_07_116/flatFiles/831439_IND_CR_FFPE_TMA_Neonate1_STJ_N_R1",
  Adult_2 = "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult2831439_08_11_2024_12_58_14_961/flatFiles/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1",
  Neonate_2 = "/mnt/scratch2/CosMX/Diego/DIegoRSVNeonate2831439_08_11_2024_13_26_50_655/flatFiles/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1"
)

# Load metadata files and merge them into the respective Seurat objects
for (slide_name in names(metadata_dirs)) {
  # Load metadata for the current slide
  metadata_file <- list.files(path = metadata_dirs[[slide_name]], pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)
  metadata <- fread(metadata_file) %>% as.data.frame()
  
  # Add Slide and cell columns to the metadata
  metadata$Slide <- slide_name
  metadata$cell <- paste0(as.character(metadata$cell_ID), "_", metadata$fov, "_", metadata$Slide)
  
  # Set rownames of the metadata to match cell identifiers
  rownames(metadata) <- metadata$cell
  
  # Filter to seurat obj
  seurat_obj <- sobj_list[[slide_name]]
  seurat_obj@meta.data$Slide <- slide_name
  seurat_obj@meta.data$cell <- paste0(as.character(seurat_obj@meta.data$cell_ID), #"_", seurat_obj@meta.data$fov,
                                      "_", seurat_obj@meta.data$Slide)
  #rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$cell
  
  metadata <- metadata[metadata$cell %in% seurat_obj@meta.data$cell, ]
  
  # Merge the metadata into the corresponding Seurat object
  sobj_list[[slide_name]] <- AddMetaData(seurat_obj, metadata)
}


## Creating a unified FOV for centroids 
# Initialize empty containers for centroids, segmentation, and molecules
centroid_list <- list()
# all_segmentation <- list()
# all_molecules <- data.frame()

# Loop through each object in the list
for (obj_name in names(sobj_list)) {
  fov <- sobj_list[[obj_name]]@images$FOV
  obj <- sobj_list[[obj_name]]
  
  # Extract and validate centroids
  centroids <- fov$centroids@coords
  if (!identical(rownames(obj@meta.data), fov$centroids@cells)) {
    stop(paste("Mismatch between meta.data rownames and FOV cells for", obj_name))
  }
  centroids <- data.frame(centroids)
  centroids$cell <- obj@meta.data$cell
  centroid_list[[obj_name]] <- centroids
  
  # # Extract segmentation boundaries
  # segmentation <- fov$segmentation@polygons$ ..... ???
  # all_segmentation <- c(all_segmentation, segmentation)
  # 
  # # Extract molecules
  # molecules <- fov$molecules@.Data ..... ???
  # molecules$cell_id <- paste0(obj_name, "_", molecules$cell_id) # Ensure unique cell IDs
  # all_molecules <- rbind(all_molecules, molecules)
}

# 2. Combine components across all objects
combined_centroids <- do.call(rbind, centroid_list)
# combined_segmentation <- do.call(c, all_segmentation)
# combined_molecules <- do.call(rbind, all_molecules)
rownames(combined_centroids) <- NULL

# 3. Create unified FOV
combined_fov <- CreateFOV(
  type = c("centroids"),
  coords = combined_centroids, 
  assay = "Nanostring",
  key = "Nanostring_",
  name = "FOV"
)

# 4. Merge Seurat objects and add combined FOV

merged_Sobj <- qs::qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Objects/all_4_slides_one_Seuobj_and_Metadata.qs")

# cell_id <- merged_Sobj@meta.data$cell
# merged_Sobj <- RenameCells(merged_Sobj, new.names = cell_id)

rownames(merged_Sobj@meta.data) <- merged_Sobj@meta.data$cell
  
merged_Sobj@images$FOV <- combined_fov
merged_Sobj@images$FOV.2 <- NULL

Sobj <- merged_Sobj


## Run Cell Type again 
load("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/Lung_MCA.RData")
metadata
cellGroups
profile_matrix[1:10, 1:10] # row: mouse genes; col: cell type; content: cell type "fraction" by gene
dim(profile_matrix)

rownames(profile_matrix)

intersect(rownames(profile_matrix), rownames(Sobj@assays$Nanostring$counts)) # matching genes (863) CosMx 1k

## Run (supervised)



counts <- data.frame(Sobj@assays$Nanostring$counts)
colnames(counts) <- merged_Sobj@meta.data$cell
sup <- insitutypeML(x = t(as.matrix(counts)),
                    neg = Sobj@meta.data$nCount_negprobes,
                    #cohort = cohort,
                    reference_profiles = as.matrix(profile_matrix)) 

Sobj@meta.data$insitutype_cluster <- NA
Sobj@meta.data$insitutype_cluster <- 'other'
Sobj$insitutype_cluster[names(sup$clust)] <- sup$clust
names(Sobj@meta.data)[names(Sobj@meta.data) == "insitutype_cluster"] <- paste0("insitutype_cluster_", "Superv","_celltypes")


table(Sobj@meta.data$insitutype_cluster_Superv_celltypes, useNA = "always")

## Renaming cell types to major cell types 
# Alveolar bipotent progenitor
Sobj@meta.data$major_celltype <- "Other"
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Alveolar.bipotent.progenitor", ]$major_celltype <- "Alveolar_bipotent_progenitor"

# Alveolar macrophage
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Alveolar.macrophage.Ear2.high", "Alveolar.macrophage.Pclaf.high"), ]$major_celltype <- "Alveolar_macrophage"

# AT1 and AT2 cells
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "AT1.cell", ]$major_celltype <- "AT1_cell"
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "AT2.cell", ]$major_celltype <- "AT2_cell"

# Lymphocytes
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("B.cell", "IgA.producing.B.cell", "NK.cell", "T.cell.Cd8b1.high"), ]$major_celltype <- "Lymphocyte"

# Granulocytes
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Basophil", "Eosinophil", "Neutrophil"), ]$major_celltype <- "Granulocyte"

# Epithelial cells
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Ciliated.cell", "Clara.cell"), ]$major_celltype <- "Epithelial_cell"

# Dendritic cells
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Dendritic.cell.Gngt2.high", "Dendritic.cell.H2.M2.high", "Dendritic.cell.Mgl2.high", "Dendritic.cell.Naaa.high", "Dendritic.cell.Tubb5.high", "Plasmacytoid.dendritic.cell"), ]$major_celltype <- "Dendritic_cell"

# Endothelial cells
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Endothelial.cell.Kdr.high", "Endothelial.cell.Tmem100.high", "Endothelial.cell.Vwf.high"), ]$major_celltype <- "Endothelial_cell"

# Interstitial macrophage
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Interstitial.macrophage", ]$major_celltype <- "Interstitial_macrophage"

# Monocyte progenitor
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Monocyte.progenitor", ]$major_celltype <- "Monocyte_progenitor"

# Stromal cells
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% 
                 c("Stromal.cell.Acta2.high", "Stromal.cell.Dcn.high", "Stromal.cell.Inmt.high"), ]$major_celltype <- "Stromal_cell"

# Nuocyte
Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Nuocyte", ]$major_celltype <- "Nuocyte"


table(Sobj@meta.data$major_celltype, useNA = "always")



## Set color map 
cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- Sobj$major_celltype %>% table() %>% names()
cols <- setNames(cols[1:length(factor_levels)], factor_levels)



## Run Niche analysis again 
# It's taking too long ...
DefaultBoundary(Sobj@images$FOV) <- "segmentation"
DefaultFOV(Sobj) <- "FOV"

Sobj <- BuildNicheAssay(object = Sobj, fov = "FOV", group.by = "major_celltype", niches.k = 5, neighbors.k = 20, cluster.name = "niches_5",
                            )
# RESUME FROM HERE 



# Let's try a downsampled version of it 
Idents(Sobj) <- "Group"
SeuObj_downsmp <- subset(x = Sobj, downsample = 500)

DefaultAssay(SeuObj_downsmp) <- "Nanostring"
SeuObj_downsmp <- BuildNicheAssay(object = SeuObj_downsmp, fov = "FOV", group.by = "major_celltype", niches.k = 5, neighbors.k = 20, cluster.name = "niches_5")

library(dittoSeq)
dittoBarPlot(SeuObj_downsmp, 
             "niches_5", 
             group.by = "Group", 
             #color.panel = cols,
             main = "")

Sobj@images$FOV$segmentation







library(dittoSeq)
comp_1_prop_niche_celltype <- dittoBarPlot(subset(Sobj, Group %in% c("Neonate IFN and RSV infected", "Adult RSV infected")), "major_celltype", 
             group.by = "niches_5",
             color.panel = cols,
             main = "")
comp_1_prop_niche_celltype




library(dittoSeq)
dittoBarPlot(Sobj, "major_celltype", 
             group.by = "niches_10",
             color.panel = cols,
             main = "")
















# Sobj = AddMetaData(
# object = Sobj_Valid_Merged,
# metadata = pData_sub)





## Test insitutype - semisupervised  example (it's working) ----------
library(InSituType)
data("ioprofiles")
data("iocolors")
data("mini_nsclc")
str(mini_nsclc)
set.seed(0)

counts <- mini_nsclc$counts
str(counts)
counts[25:30, 9:14]

negmean <- Matrix::rowMeans(mini_nsclc$neg)
head(negmean)

data(ioprofiles)
str(ioprofiles)
ioprofiles[1:5, 1:10]

astats <- get_anchor_stats(counts = mini_nsclc$counts,
                           neg = Matrix::rowMeans(mini_nsclc$neg),
                           profiles = ioprofiles)
#> The following genes in the count data are missing from fixed_profiles and will be omitted from anchor selection: CCL15,CCL18,CCL23,CCL3,CCL3L3,CCL4,CCL4L2,CCL5,CD24,CYTOR,FCGBP,FYB1,GDF10,H2AZ1,H4C3,LINC02446,MMP12,MRC1,MSMB,PECAM1,PRSS2,SELENOP,VEGFD,VSIR

# estimate per-cell bg as a fraction of total counts:
negmean.per.totcount <- mean(rowMeans(mini_nsclc$neg)) / mean(rowSums(counts))
per.cell.bg <- rowSums(counts) * negmean.per.totcount

# now choose anchors:
anchors <- choose_anchors_from_stats(counts = counts, 
                                     neg = mini_nsclc$negmean, 
                                     bg = per.cell.bg,
                                     anchorstats = astats, 
                                     # a very low value chosen for the mini
                                     # dataset. Typically hundreds of cells
                                     # would be better.
                                     n_cells = 50, 
                                     min_cosine = 0.4, 
                                     min_scaled_llr = 0.03, 
                                     insufficient_anchors_thresh = 5)
#> The following cell types had too few anchors and so are being removed from consideration: mast, mDC, monocyte, neutrophil, NK, pDC, T CD4 memory, T CD4 naive, T CD8 memory, T CD8 naive, Treg

# plot the anchors atop the UMAP:
par(mfrow = c(1, 1))
plot(mini_nsclc$umap, pch = 16, cex = 0.1, col = "peachpuff1", xaxt = "n",  yaxt = "n", xlab = "", ylab = "",
     main = "Selected anchor cells")
points(mini_nsclc$umap[!is.na(anchors), ], col = iocolors[anchors[!is.na(anchors)]], pch = 16, cex = 0.6)
legend("topright", pch = 16, col = iocolors[setdiff(unique(anchors), NA)], legend = setdiff(unique(anchors), NA), cex = 0.65)

updatedprofiles <- updateReferenceProfiles(reference_profiles = ioprofiles, 
                                           counts = mini_nsclc$counts, 
                                           neg = mini_nsclc$neg, 
                                           bg = per.cell.bg,
                                           anchors = NULL) # one may want to work on anchors argument 

immunofluordata <- matrix(rpois(n = nrow(counts) * 4, lambda = 100), 
                          nrow(counts))
cohort <- fastCohorting(immunofluordata,
                        gaussian_transform = TRUE) 
table(cohort)


updatedprofiles$updated_profiles


semisup <- insitutype(
  x = counts,
  neg = negmean,
  cohort = cohort,
  # Enter your own per-cell background estimates here if you
  # have them; otherwise insitutype will use the negprobes to
  # estimate background for you.
  bg = NULL,
  # condensed to save time. n_clusts = 5:15 would be more optimal
  n_clusts = c(5, 6),
  reference_profiles = updatedprofiles$updated_profiles,
  update_reference_profiles = FALSE,
  # choosing inadvisably low numbers to speed the vignette; using the defaults
  # in recommended.
  n_phase1 = 200,
  n_phase2 = 500,
  n_phase3 = 2000,
  n_starts = 1,
  max_iters = 5
) 
#> 6.74 cells per geometric bin.
#> Selecting optimal number of clusters from a range of 5 - 6
#> Clustering with n_clust = 5
#> iter 1
#> iter 2
#> iter 3
#> iter 4
#> iter 5
#> iter 6
#> Converged: <= 0.5% of cell type assignments changed in the last iteration.
#> ==========================================================================
#> Clustering with n_clust = 6
#> iter 1
#> iter 2
#> iter 3
#> iter 4
#> iter 5
#> iter 6
#> iter 7
#> iter 8
#> Converged: <= 0.5% of cell type assignments changed in the last iteration.
#> ==========================================================================
#> phase 1: random starts in 200 cell subsets
#> iter 1
#> iter 2
#> iter 3
#> Converged: <= 0.2% of cell type assignments changed in the last iteration.
#> ==========================================================================
#> phase 2: refining best random start in a 500 cell subset
#> iter 1
#> iter 2
#> iter 3
#> iter 4
#> iter 5
#> phase 3: finalizing clusters in a 2000 cell subset
#> iter 1
#> iter 2
#> iter 3
#> iter 4
#> iter 5
#> phase 4: classifying all 2198 cells





