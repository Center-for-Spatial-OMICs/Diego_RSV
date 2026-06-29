

## Current using -------
library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(qs)

options(future.globals.maxSize = 99999 * 1024^2)

# Load the Seurat object list
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

# Start with the first object
All_Sobj <- sobj_list[[1]]

# Loop through the remaining objects and merge them one by one
for (i in 2:length(sobj_list)) {
  All_Sobj <- merge(All_Sobj, y = sobj_list[[i]])
}

All_Sobj <- JoinLayers(All_Sobj)

qsave(All_Sobj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Objects/all_4_slides_one_Seuobj_and_Metadata.qs")









## Out of usage ------- 

# library(Seurat); library(tidyverse); library(data.table); library(dplyr)
# library(RColorBrewer); library(InSituType)
# library(ggplot2)
# library(dplyr)
# library(viridis)
# library(dittoSeq)
# library(qs)
# 
# options(future.globals.maxSize = 99999 * 1024^2)
# 
# ## Combine seurat obj
# sobj_list <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult12_Neonate12_seuobj_list.rds")
# 
# sobj_list$Adult_1$Slide <- "Adult_1"
# sobj_list$Adult_1$cell <- paste0(as.character(sobj_list$Adult_1$cell_ID), "_", sobj_list$Adult_1$fov, "_",sobj_list$Adult_1$Slide)
# #rownames(sobj_list$Adult_1@meta.data) <- sobj_list$Adult_1$cell
# 
# sobj_list$Adult_2$Slide <- "Adult_2"
# sobj_list$Adult_2$cell <- paste0(as.character(sobj_list$Adult_2$cell_ID), "_", sobj_list$Adult_2$fov, "_", sobj_list$Adult_2$Slide)
# #rownames(sobj_list$Adult_2@meta.data) <- sobj_list$Adult_2$cell
# 
# sobj_list$Neonate_1$Slide <- "Neonate_1"
# sobj_list$Neonate_1$cell <- paste0(as.character(sobj_list$Neonate_1$cell_ID), "_", sobj_list$Neonate_1$fov, "_", sobj_list$Neonate_1$Slide)
# #rownames(sobj_list$Neonate_1@meta.data) <- sobj_list$Neonate_1$cell
# 
# sobj_list$Neonate_2$Slide <- "Neonate_2"
# sobj_list$Neonate_2$cell <- paste0(as.character(sobj_list$Neonate_2$cell_ID), "_", sobj_list$Neonate_2$fov, "_",sobj_list$Neonate_2$Slide)
# #rownames(sobj_list$Neonate_2@meta.data) <- sobj_list$Neonate_2$cell
# 
# 
# All_Sobj <- Reduce(function(x, y) merge(x, y), sobj_list)
# All_Sobj <- JoinLayers(All_Sobj)
# 
# qsave(All_Sobj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Objects/all_4_slides_one_Seuobj.qs")
# 
# 
# 
# ## Combine metadata
# # Path dir to flat files
# Adult_1_dir <- "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult1831439_08_11_2024_9_49_46_243/flatFiles/831439_IND_CR_FFPE_TMA_Adult1_STJ_N_R1"
# 
# Neonate_1_dir <- "/mnt/scratch2/CosMX/Diego/DiegoRSVNeonate1831439_08_11_2024_11_18_07_116/flatFiles/831439_IND_CR_FFPE_TMA_Neonate1_STJ_N_R1"
# 
# Adult_2_dir <- "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult2831439_08_11_2024_12_58_14_961/flatFiles/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1"
# 
# Neonate_2_dir <- "/mnt/scratch2/CosMX/Diego/DIegoRSVNeonate2831439_08_11_2024_13_26_50_655/flatFiles/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1"
# 
# # Load and combine metadata
# metadata_list <- list(
#   fread(list.files(path = Adult_1_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)),
#   fread(list.files(path = Neonate_1_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)),
#   fread(list.files(path = Adult_2_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)),
#   fread(list.files(path = Neonate_2_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE))
# )
# 
# # Define the names for each element in the list
# names(metadata_list) <- c("Adult_1", "Neonate_1", "Adult_2", "Neonate_2")
# 
# # Apply the operation to create the 'cell' column for all metadata files
# metadata_list <- lapply(names(metadata_list), function(name) {
#   metadata <- metadata_list[[name]]
#   metadata$Slide <- name
#   metadata$cell <- paste0(as.character(metadata$cell_ID), "_", metadata$fov, "_", metadata$Slide) #(cell_ID + fov + source slide) == cell
#   return(metadata)
# })
# 
# # Combining all metadata into one
# metadata <- plyr::rbind.fill(metadata_list)
# head(metadata)
# 
# qsave(metadata, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Objects/all_4_slides_one_Metadata.qs")


