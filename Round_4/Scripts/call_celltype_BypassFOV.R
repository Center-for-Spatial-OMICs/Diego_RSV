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

### Celltype from the merged objects by Reduce()


### Try to use Reduce specifying "cell_ID" - use a different script


### Fixing FOV -----------

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

merged_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj.qs")

# cell_id <- merged_Sobj@meta.data$cell
# merged_Sobj <- RenameCells(merged_Sobj, new.names = cell_id)

rownames(merged_Sobj@meta.data) <- merged_Sobj@meta.data$cell

merged_Sobj@images$FOV <- combined_fov
merged_Sobj@images$FOV.2 <- NULL

Sobj <- merged_Sobj




### Running Niche  -----------
rownames(Sobj)
colnames(Sobj) 
Sobj@assays$Nanostring$counts %>% rownames()
Sobj@assays$Nanostring$counts %>% colnames()

# cell ID is different across the whole seurat object
# Not working
Sobj <- RenameCells(Sobj, new.names = Sobj@meta.data$cell)
# Not working
colnames(Sobj@assays$Nanostring$counts) <- Sobj@meta.data$cell 

DefaultAssay(Sobj) <- "Nanostring"
Idents(Sobj) <- "Group"
SeuObj_downsmp <- subset(x = Sobj, downsample = 500)


## ERRORR .......

DefaultAssay(SeuObj_downsmp) <- "Nanostring"
SeuObj_downsmp <- BuildNicheAssay(object = SeuObj_downsmp, fov = "FOV", group.by = "major_celltype", niches.k = 5, neighbors.k = 20, cluster.name = "niches_5")

library(dittoSeq)
dittoBarPlot(SeuObj_downsmp, 
             "niches_5", 
             group.by = "Group", 
             #color.panel = cols,
             main = "")







