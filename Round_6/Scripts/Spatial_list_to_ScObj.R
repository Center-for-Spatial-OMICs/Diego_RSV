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


# ## Load objects -----------
# obj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/Sobj_list_TMA_ID_fixed.qs")
# 
# 
# meta_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/metadata_list.qs")
# 
# if (!identical(names(obj_list), names(meta_list))) {
#   stop("Names of obj_list and meta_list do not match!")
# }
# obj_list <- lapply(names(obj_list), function(name) {
#   AddMetaData(obj_list[[name]], metadata = meta_list[[name]])
# })
# # Putting list element's names back
# names(obj_list) <- names(meta_list)
# 
# 
# 
# 
# All_Sobj <- obj_list[[1]]
# 
# # Loop through the remaining objects and merge them one by one
# for (i in 2:length(obj_list)) {
#   obj_name <- names(obj_list)[i]
#   
#   if (!is.na(obj_name) && obj_name != "NA") {
#     All_Sobj <- merge(All_Sobj, y = obj_list[[i]])
#     print(paste0(obj_name, " has been merged (", round(i/length(obj_list), 2), ")"))
#     qsave(All_Sobj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/ScObj_merging.qs")
#   } else {
#     print(paste0("Skipping index ", i, " (name is NA)"))
#   }
# }
# 
# # All_Sobj <- Reduce(function(x, y) merge(x, y), obj_list)
# 
# print("Starting JoinLayers")
# All_Sobj <- JoinLayers(All_Sobj, assay = "Nanostring")
# 
# 
# print("Saving All_Sobj")
# qsave(All_Sobj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/ScObj.qs")


All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/ScObj.qs")


# ##### Troubleshooting -----
# Idents(All_Sobj) <- "major_celltype"
# obj <- subset(All_Sobj, downsample = 100)
# 
# obj <- obj %>%
#   NormalizeData() %>% 
#   FindVariableFeatures(nfeatures = nrow(obj@assays$Nanostring$counts)) %>%
#   ScaleData() %>%
#   RunPCA() %>%
#   FindNeighbors(dims = 1:10) %>%
#   #FindClusters() %>% #already have insitutype 
#   RunUMAP(dims = 1:10)
# 
# #####


All_Sobj <- All_Sobj %>%
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = nrow(All_Sobj@assays$Nanostring$counts)) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:10) %>%
  #FindClusters() %>% #already have insitutype 
  RunUMAP(dims = 1:10)

print("Saving All_Sobj + UMAP")
qsave(All_Sobj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/ScObj_umap.qs")


