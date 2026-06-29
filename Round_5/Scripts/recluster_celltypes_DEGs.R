## ====== Date: Jun 02, 2025  ===== 

# ========================================================
## Directions: 
## DEG Analysis in Re-clustered Cell Types: Neonate NO IFN infected vs Neonate IFN infected vs Adult infected
## Re-cluster separated cell types 
# Dendritic
# Epithelial
# Lymphocytes
## Get Volcano plots for their DEGs
# Eg: Group1_dendritic vs Group2_dendritic vs Group3_dendritic
# ========================================================


### Load and define functions ---------
library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(qs)

### Load and define variables  ---------
Sobj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj_list_celltype_niche.qs")

# All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_merged_Sobj.qs")
All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_processed.qs")


comparisons <- list(
  "Comparison" = list(
    group_1 = "Neonate RSV infected (NO IFN)",
    group_2 = "Neonate IFN and RSV infected",
    group_3 = "Adult RSV infected"
  )
)

### Re-clustering based on cell type  ---------
Sobj_den <- subset(All_Sobj, major_celltype %in% "Dendritic_cell")
Sobj_epi <- subset(All_Sobj, major_celltype %in% "Epithelial_cell")
Sobj_lyn <- subset(All_Sobj, major_celltype %in% "Lymphocyte")


## Epithelial ----
Sobj_epi@assays$RNA <- Sobj_epi@assays$Nanostring
obj <- Sobj_epi
DefaultAssay(obj) <- "RNA"
obj <- obj %>%
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = length(rownames(obj))) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = c(0.2, 0.5, 1, 1.5, 2.0)) %>% 
  RunUMAP(dims = 1:15)

Idents(obj) <- "RNA_snn_res.0.5"
p1 <- DimPlot(obj)

Idents(obj) <- "RNA_snn_res.1"
p2 <- DimPlot(obj)

Idents(obj) <- "RNA_snn_res.0.2"
p3 <- DimPlot(obj)

p3 |p1 | p2 


library(data.table)
library(dplyr)
library(clustree)

cluster_df_sub <- obj@meta.data[, c("RNA_snn_res.0.2", 
                                    "RNA_snn_res.0.5", 
                                    "RNA_snn_res.1",
                                    "RNA_snn_res.1.5",
                                    "RNA_snn_res.2")]

# Change column names
names(cluster_df_sub)
names(cluster_df_sub) <- c("seurat_row_1", 
                           "seurat_row_2", 
                           "seurat_row_3",
                           "seurat_row_4",
                           "seurat_row_5") 

library(clustree)
clustree_out <- clustree(cluster_df_sub[, ], prefix = "seurat_row_", node_colour = "sc3_stability")

print(clustree_out)
# Decreases cluster stability 
# <---------------
  
## Manual annotation based on clustree_out
obj@meta.data$Resolved_clusters <- NA
obj@meta.data[obj@meta.data$RNA_snn_res.1.5 %in% c("11"), ]$Resolved_clusters <- "11_res_1.5"
obj@meta.data[obj@meta.data$RNA_snn_res.1.5 %in% c("6"), ]$Resolved_clusters <- "6_res_1.5"

obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("0"), ]$Resolved_clusters <- "0_res_0.5"
obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("2"), ]$Resolved_clusters <- "2_res_0.5"
obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("6"), ]$Resolved_clusters <- "6_res_0.5"

obj@meta.data[obj@meta.data$RNA_snn_res.1 %in% c("3"), ]$Resolved_clusters <- "3_res_1"
obj@meta.data[obj@meta.data$RNA_snn_res.1 %in% c("7"), ]$Resolved_clusters <- "7_res_1"
obj@meta.data[obj@meta.data$RNA_snn_res.1 %in% c("9"), ]$Resolved_clusters <- "9_res_1"
obj@meta.data[obj@meta.data$RNA_snn_res.1 %in% c("10"), ]$Resolved_clusters <- "10_res_1"

obj@meta.data[obj@meta.data$RNA_snn_res.1 %in% c("1"), ]$Resolved_clusters <- "1_res_1"
obj@meta.data[obj@meta.data$RNA_snn_res.1 %in% c("6"), ]$Resolved_clusters <- "6_res_1"

obj@meta.data[obj@meta.data$RNA_snn_res.0.2 %in% c("5"), ]$Resolved_clusters <- "5_res_0.2"

obj@meta.data[obj@meta.data$Resolved_clusters %in% NA, ]$Resolved_clusters <- "non_resolved_clusters"

table(obj@meta.data$Resolved_clusters, useNA = "always")



Idents(obj) <- "Resolved_clusters"
p4 <- DimPlot(obj)

p3 |p1 | p2 | p4

p4

# Remove "non_resolved_clusters"
obj_bk <- obj
obj <- subset(obj, Resolved_clusters %in% "non_resolved_clusters", invert = TRUE)

## Resolved_clusters cluster markers 
library(presto)
markers <- wilcoxauc(obj, 'Resolved_clusters')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 5) %>% 
  ungroup() -> top_markers

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "Resolved_clusters"
DotPlot(obj, 
        #features = unique(top_markers$feature)
        features =  unique(c(top_markers$feature))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) + ggtitle("")


Idents(obj) <- "Group"
p5 <- DimPlot(obj)
p5 

obj_test <- subset(obj, Group %in% c(comparisons$Comparison))
Idents(obj_test) <- "Group"
p6 <- DimPlot(obj_test)
p6 

Idents(obj_test) <- "Resolved_clusters"
p7 <- DimPlot(obj_test)
p7 

obj_test@meta.data$Resolved_clusters_good_markers <- "other"
obj_test@meta.data[obj_test@meta.data$Resolved_clusters %in% 
                     c("6_res_1",
                       "5_res_0.2",
                       "2_res_0.5", 
                       "10_res_1"), ]

# qsave(obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_epithelial.qs")

## Summary:
## 2_res_0.5 sounds like a real phenotype and it's across the three groups of interest. Fewer cells on Neonate IFN and RSV infected compared to Neonate RSV infected (NO IFN) and Adult RSV infected.



## ====== Date: Jun 03, 2025  ===== 

## Lymphocyte ----
Sobj_lyn@assays$RNA <- Sobj_lyn@assays$Nanostring
obj <- Sobj_lyn
DefaultAssay(obj) <- "RNA"
obj <- obj %>%
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = length(rownames(obj))) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = c(0.2, 0.5, 1, 1.5, 2.0)) %>% 
  RunUMAP(dims = 1:15)

Idents(obj) <- "RNA_snn_res.0.5"
p1 <- DimPlot(obj)

Idents(obj) <- "RNA_snn_res.1"
p2 <- DimPlot(obj)

Idents(obj) <- "RNA_snn_res.0.2"
p3 <- DimPlot(obj)

p3 |p1 | p2 


library(data.table)
library(dplyr)
library(clustree)

cluster_df_sub <- obj@meta.data[, c("RNA_snn_res.0.2", 
                                    "RNA_snn_res.0.5", 
                                    "RNA_snn_res.1",
                                    "RNA_snn_res.1.5",
                                    "RNA_snn_res.2")]

# Change column names
names(cluster_df_sub)
names(cluster_df_sub) <- c("seurat_row_1", 
                           "seurat_row_2", 
                           "seurat_row_3",
                           "seurat_row_4",
                           "seurat_row_5") 

library(clustree)
clustree_out <- clustree(cluster_df_sub[, ], prefix = "seurat_row_", node_colour = "sc3_stability")

print(clustree_out)
# Decreases cluster stability 
# <---------------

## Manual annotation based on clustree_out
obj@meta.data$Resolved_clusters <- NA

obj@meta.data[obj@meta.data$RNA_snn_res.1 %in% c("0"), ]$Resolved_clusters <- "0_res_1"

obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("4"), ]$Resolved_clusters <- "4_res_0.5"
obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("5"), ]$Resolved_clusters <- "5_res_0.5"
obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("6"), ]$Resolved_clusters <- "6_res_0.5"

obj@meta.data[obj@meta.data$RNA_snn_res.0.2 %in% c("1"), ]$Resolved_clusters <- "1_res_0.2_TooComplex"
obj@meta.data[obj@meta.data$RNA_snn_res.0.2 %in% c("2"), ]$Resolved_clusters <- "2_res_0.2_TooComplex"

obj@meta.data[obj@meta.data$RNA_snn_res.0.2 %in% c("4"), ]$Resolved_clusters <- "4_res_0.2"
obj@meta.data[obj@meta.data$RNA_snn_res.0.2 %in% c("5"), ]$Resolved_clusters <- "5_res_0.2"
obj@meta.data[obj@meta.data$RNA_snn_res.0.2 %in% c("6"), ]$Resolved_clusters <- "6_res_0.2"
obj@meta.data[obj@meta.data$RNA_snn_res.0.2 %in% c("7"), ]$Resolved_clusters <- "7_res_0.2"

obj@meta.data[obj@meta.data$Resolved_clusters %in% NA, ]$Resolved_clusters <- "non_resolved_clusters"

table(obj@meta.data$Resolved_clusters, useNA = "always")


Idents(obj) <- "Resolved_clusters"
p4 <- DimPlot(obj)

p3 |p1 | p2 | p4

p4

# Remove "non_resolved_clusters"
#obj_bk <- obj
obj <- subset(obj, Resolved_clusters %in% "non_resolved_clusters", invert = TRUE)

## Resolved_clusters cluster markers 
library(presto)
markers <- wilcoxauc(obj, 'Resolved_clusters')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 5) %>% 
  ungroup() -> top_markers

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "Resolved_clusters"
DotPlot(obj, 
        #features = unique(top_markers$feature)
        features =  unique(c(top_markers$feature))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) + ggtitle("")


Idents(obj) <- "Group"
p5 <- DimPlot(obj)
p5 

obj_test <- subset(obj, Group %in% c(comparisons$Comparison))
Idents(obj_test) <- "Group"
p6 <- DimPlot(obj_test)
p6 

Idents(obj_test) <- "Resolved_clusters"
p7 <- DimPlot(obj_test)
p7 

obj_test@meta.data$Resolved_clusters_good_markers <- obj_test@meta.data$Resolved_clusters
obj_test@meta.data[!obj_test@meta.data$Resolved_clusters %in% 
                     c("4_res_0.2",
                       "5_res_0.2",
                       "6_res_0.2", 
                       "6_res_0.2",
                       "4_res_0.5",
                       "5_res_0.5"), ]$Resolved_clusters_good_markers <- "other"

Idents(obj_test) <- "Resolved_clusters_good_markers"
p7 <- DimPlot(obj_test)
p7 

table(obj_test@meta.data$Resolved_clusters_good_markers,
      obj_test@meta.data$Group)


# qsave(obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_lymphocytes.qs")


## Summary:
## 4_res_0.5 is the most analyzble sub-cluster among those 3 groups 


## Dendritic ----
Sobj_den@assays$RNA <- Sobj_den@assays$Nanostring
obj <- Sobj_den
DefaultAssay(obj) <- "RNA"
obj <- obj %>%
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = length(rownames(obj))) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = c(0.2, 0.5, 1, 1.5, 2.0)) %>% 
  RunUMAP(dims = 1:15)

Idents(obj) <- "RNA_snn_res.0.5"
p1 <- DimPlot(obj)

Idents(obj) <- "RNA_snn_res.1"
p2 <- DimPlot(obj)

Idents(obj) <- "RNA_snn_res.0.2"
p3 <- DimPlot(obj)

p3 |p1 | p2 


library(data.table)
library(dplyr)
library(clustree)

cluster_df_sub <- obj@meta.data[, c("RNA_snn_res.0.2", 
                                    "RNA_snn_res.0.5", 
                                    "RNA_snn_res.1",
                                    "RNA_snn_res.1.5",
                                    "RNA_snn_res.2")]

# Change column names
names(cluster_df_sub)
names(cluster_df_sub) <- c("seurat_row_1", 
                           "seurat_row_2", 
                           "seurat_row_3",
                           "seurat_row_4",
                           "seurat_row_5") 

library(clustree)
clustree_out <- clustree(cluster_df_sub[, ], prefix = "seurat_row_", node_colour = "sc3_stability")

print(clustree_out)
# Decreases cluster stability 
# <---------------

## Manual annotation based on clustree_out
obj@meta.data$Resolved_clusters <- NA

obj@meta.data[obj@meta.data$RNA_snn_res.1 %in% c("0"), ]$Resolved_clusters <- "0_res_1"

obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("2"), ]$Resolved_clusters <- "2_res_0.5_TooComplex"
obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("3"), ]$Resolved_clusters <- "3_res_0.5_TooComplex"
obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("6"), ]$Resolved_clusters <- "6_res_0.5_TooComplex"

obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("1"), ]$Resolved_clusters <- "1_res_0.5"
obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("4"), ]$Resolved_clusters <- "4_res_0.5"
obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("7"), ]$Resolved_clusters <- "7_res_0.5"
obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("8"), ]$Resolved_clusters <- "8_res_0.5"
obj@meta.data[obj@meta.data$RNA_snn_res.0.5 %in% c("9"), ]$Resolved_clusters <- "9_res_0.5"

obj@meta.data[obj@meta.data$RNA_snn_res.0.2 %in% c("3"), ]$Resolved_clusters <- "3_res_0.2"

obj@meta.data[obj@meta.data$Resolved_clusters %in% NA, ]$Resolved_clusters <- "non_resolved_clusters"

table(obj@meta.data$Resolved_clusters, useNA = "always")


Idents(obj) <- "Resolved_clusters"
p4 <- DimPlot(obj)

p3 |p1 | p2 | p4

p4

# Remove "non_resolved_clusters"
#obj_bk <- obj
obj <- subset(obj, Resolved_clusters %in% "non_resolved_clusters", invert = TRUE)

## Resolved_clusters cluster markers 
library(presto)
markers <- wilcoxauc(obj, 'Resolved_clusters')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 5) %>% 
  ungroup() -> top_markers

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "Resolved_clusters"
DotPlot(obj, 
        #features = unique(top_markers$feature)
        features =  unique(c(top_markers$feature))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) + ggtitle("")


Idents(obj) <- "Group"
p5 <- DimPlot(obj)
p5 

obj_test <- subset(obj, Group %in% c(comparisons$Comparison))
Idents(obj_test) <- "Group"
p6 <- DimPlot(obj_test)
p6 

Idents(obj_test) <- "Resolved_clusters"
p7 <- DimPlot(obj_test)
p7 

obj_test@meta.data$Resolved_clusters_good_markers <- obj_test@meta.data$Resolved_clusters
obj_test@meta.data[!obj_test@meta.data$Resolved_clusters %in% 
                     c("9_res_0.5",
                       "7_res_0.5",
                       "8_res_0.5",
                       "3_res_0.5"), ]$Resolved_clusters_good_markers <- "other"

Idents(obj_test) <- "Resolved_clusters_good_markers"
p7 <- DimPlot(obj_test)
p7 

table(obj_test@meta.data$Resolved_clusters_good_markers,
      obj_test@meta.data$Group)


#qsave(obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_dendrocytes.qs")



## Summary:
## 8_res_0.5 is the most analyzble sub-cluster among those 3 groups 







### Visualizing the highlighted sub-clusters on spatial -------

## Naming back sub-clusters
# sc_obj Dendrocytes
obj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_dendrocytes.qs")


obj_test <- subset(obj, Group %in% c(comparisons$Comparison))

dittoBarPlot(obj_test,
             var = "Resolved_clusters",
             group.by = "TMA",
             split.by = "Group")

obj_test@meta.data$Resolved_clusters_good_markers <- obj_test@meta.data$Resolved_clusters
obj_test@meta.data[!obj_test@meta.data$Resolved_clusters %in% 
                     c("9_res_0.5",
                       "7_res_0.5",
                       "8_res_0.5",
                       "3_res_0.5"), ]$Resolved_clusters_good_markers <- "other"


dittoBarPlot(obj_test,
             var = "Resolved_clusters_good_markers",
             group.by = "TMA",
             split.by = "Group")




## Looking for on sub-cluster
obj_test <- subset(obj, Group %in% c(comparisons$Comparison) &
                     Resolved_clusters %in% "8_res_0.5")

# Checking if there's any duplicated cell ID within each slide 
duplicated(obj_test@meta.data[obj_test@meta.data$Slide %in% "Adult_1", ]$cell_id,
           Sobj_list$Adult_1@meta.data$cell_id) %>% table()
duplicated(
  Sobj_list$Adult_1@meta.data$cell_id,
  obj_test@meta.data[obj_test@meta.data$Slide %in% "Adult_1", ]$cell_id) %>% table()



label_back_to_spatial <- function(obj_meta, Sobj_list, slide_names, column_name,label_name) {
  
  for(slide in slide_names) {
    # Get cell IDs from current slide in obj_test
    cells_to_anno <- obj_meta@meta.data[obj_meta@meta.data$Slide %in% slide, ]$cell_id
    
    # Get cell IDs from corresponding Sobj_list element
    Sobj_list[[slide]]@meta.data[[column_name]] <- "other"
    Sobj_list[[slide]]@meta.data[Sobj_list[[slide]]@meta.data$cell_id %in% cells_to_anno, ] %>% dim() %>% print()
    Sobj_list[[slide]]@meta.data[Sobj_list[[slide]]@meta.data$cell_id %in% cells_to_anno, ][[column_name]] <- label_name
  }
  
  # Convert back to data.frame format
  return(Sobj_list)
}


# Label dendritic cells
Sobj_list <- label_back_to_spatial(
  obj_meta = obj_test,
  Sobj_list = Sobj_list,
  slide_names = c("Adult_1", "Adult_2", "Neonate_1", "Neonate_2"),
  columns_name = "Cell_type_subcluster",
  label_name = "Dendritic_8_res_0.5"
)

table(Sobj_list$Adult_1$TMA,
      Sobj_list$Adult_1$Cell_type_subcluster)

table(Sobj_list$Adult_1$TMA_2,
      Sobj_list$Adult_1$Cell_type_subcluster)

ImageDimPlot(
  #subset(Sobj_list$Adult_1, Cell_type_subcluster %in% "Dendritic_8_res_0.5"),
  # subset(Sobj_list$Adult_1, TMA_2 %in% c("Ap16-3779", "Ap16-3780", "Ap16-3781")),
  subset(Sobj_list$Adult_1, TMA_2 %in% c("Ap16-3780")),
  #Sobj_list$Adult_1,
  group.by = 'Cell_type_subcluster',
  fov = "FOV",
  size = 0.5,
  #cols = cluster_colors, # I don't want it to match now cause the clusters were calculated separated 
  border.size = 0.01
) + ggtitle("")



# sc_obj Epithelial 
obj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_epithelial.qs")

obj_test <- subset(obj, Group %in% c(comparisons$Comparison))

dittoBarPlot(obj_test,
             var = "Resolved_clusters",
             group.by = "TMA",
             split.by = "Group")

obj_test@meta.data$Resolved_clusters_good_markers <- obj_test@meta.data$Resolved_clusters
obj_test@meta.data[!obj_test@meta.data$Resolved_clusters %in% 
                     c("6_res_1",
                       "5_res_0.2",
                       "2_res_0.5", 
                       "10_res_1"), ]$Resolved_clusters_good_markers <- "other"

dittoBarPlot(obj_test,
             var = "Resolved_clusters_good_markers",
             group.by = "TMA",
             split.by = "Group")






# sc_obj Lymphocytes
obj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_lymphocytes.qs")


obj_test <- subset(obj, Group %in% c(comparisons$Comparison))

dittoBarPlot(obj_test,
             var = "Resolved_clusters",
             group.by = "TMA",
             split.by = "Group")

obj_test@meta.data$Resolved_clusters_good_markers <- obj_test@meta.data$Resolved_clusters
obj_test@meta.data[!obj_test@meta.data$Resolved_clusters %in% 
                     c("4_res_0.2",
                       "5_res_0.2",
                       "6_res_0.2", 
                       "6_res_0.2",
                       "4_res_0.5",
                       "5_res_0.5"), ]$Resolved_clusters_good_markers <- "other"


dittoBarPlot(obj_test,
             var = "Resolved_clusters_good_markers",
             group.by = "TMA",
             split.by = "Group")


obj_test <- subset(obj, Group %in% c(comparisons$Comparison) &
                     Resolved_clusters %in% "4_res_0.5")

Sobj_list <- label_back_to_spatial(
  obj_meta = obj_test,
  Sobj_list = Sobj_list,
  slide_names = c("Adult_1", "Adult_2", "Neonate_1", "Neonate_2"),
  column_name = "Cell_type_subcluster",
  label_name = "Lymphocyte_4_res_0.5"
)


 
ImageDimPlot(
  #subset(Sobj_list$Adult_1, Cell_type_subcluster %in% "Dendritic_8_res_0.5"),
  # subset(Sobj_list$Adult_1, TMA_2 %in% c("Ap16-3779", "Ap16-3780", "Ap16-3781")),
  subset(Sobj_list$Adult_1, TMA_2 %in% c("Ap16-3780")),
  #Sobj_list$Adult_1,
  group.by = 'Cell_type_subcluster',
  fov = "FOV",
  size = 0.5,
  #cols = cluster_colors, # I don't want it to match now cause the clusters were calculated separated 
  border.size = 0.01
) + ggtitle("")



Idents(obj) <- "Resolved_clusters"
DotPlot(obj, features =  unique(c("Ptprc", "Cd3e", "Cd4", #TCD4
                                  "Ptprc", "Cd3e", "Cd8a" #TCD8
                                  ))[]) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  ggtitle('TCD4 and TCD8 cell markers within lymphocytes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))









### Volcano DEGs (do we really need it for analysis?) -------
# Wanna force comparisons? 
# Dendritic Group 1 vs Dendritic Group 2 vs Dendrict Group 3 ?? 

# what if our cell type ref is not that tight? Let's annotate/check the current cell anno with Diego's markers before forcing differences 



All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_merged_Sobj.qs")

### Processing --------
obj <- All_Sobj
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- SketchData(
  object = obj,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = length(rownames(obj)))
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = c(0.2, 0.5, 1, 1.5, 2))
obj <- RunUMAP(obj, dims = 1:30, return.model = T)

obj@assays$RNA <- obj@assays$Nanostring
DefaultAssay(obj) <- "sketch"
obj <- ProjectData(
  object = obj,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:30,
  refdata = list(cluster_0.2_full = "sketch_snn_res.0.2",
                 cluster_0.5_full = "sketch_snn_res.0.5",
                 cluster_1.0_full = "sketch_snn_res.1", 
                 cluster_1.5_full = "sketch_snn_res.1.5",
                 cluster_2.0_full = "sketch_snn_res.2")
)
# switch back to analyzing all cells
DefaultAssay(obj) <- "RNA"

Idents(obj) <- "cluster_0.5_full"
DimPlot(obj)

# qsave(obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_processed.qs")


Idents(obj) <- "cluster_0.5_full"
DotPlot(obj, features =  unique(c("Ptprc", "Cd3e", "Cd4", #TCD4
                                  "Ptprc", "Cd3e", "Cd8a" #TCD8
))[]) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  ggtitle('TCD4 and TCD8 cell markers within all cells') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Idents(obj) <- "cluster_0.2_full"
DotPlot(obj, features =  unique(c("Ptprc", "Cd3e", "Cd4", #TCD4
                                  "Ptprc", "Cd3e", "Cd8a" #TCD8
))[]) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  ggtitle('TCD4 and TCD8 cell markers within all cells') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


Idents(obj) <- "cluster_1.5_full"
DotPlot(obj, features =  unique(c("Ptprc", "Cd3e", "Cd4", #TCD4
                                  "Ptprc", "Cd3e", "Cd8a" #TCD8
))[]) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  ggtitle('TCD4 and TCD8 cell markers within all cells') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


Idents(obj) <- "cluster_2.0_full"
DotPlot(obj, features =  unique(c("Ptprc", "Cd3e", "Cd4", #TCD4
                                  "Ptprc", "Cd3e", "Cd8a" #TCD8
))[]) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  ggtitle('TCD4 and TCD8 cell markers within all cells') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


Idents(obj) <- "major_celltype"
DotPlot(obj, features =  unique(c("Ptprc", "Cd3e", "Cd4", #TCD4
                                  "Ptprc", "Cd3e", "Cd8a" #TCD8
))[]) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  ggtitle('TCD4 and TCD8 cell markers within all cells') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Summary: 
# It seems the markers Diego has givenus can’t resolve T CD4 any better  than the cell anno. we already have Therefore, I believe we’re good to go forward (DEG, spatial cor. etc) with the cell anno wehave used so far.




## Differential Expression Analysis ------

Plot_Volcano_DEG_sc <- function(SeuratObj = NULL, Identity = NULL, Group_1 = NULL, Group_2 = NULL, markers_Diego) {
  # Subset to comparison 
  Idents(SeuratObj) <- Identity
  SeuratObj_sub <- subset(x = SeuratObj, idents = c(Group_1, Group_2), invert = FALSE)
  
  
  # Do Differential Expression Analysis - find all markers 
  library(presto)
  markers <- presto::wilcoxauc(SeuratObj_sub, group_by = Identity)
  # markers <- markers[markers$padj < 0.05 & markers$auc > 0.80, ]
  markers <- markers[order(markers$logFC, decreasing = TRUE), ]
  markers <- markers[markers$padj <= 0.05, ]
  markers_Group_2 <- markers[markers$group %in% Group_2, ]
  
  # Define top and bottom DEGs
  DEGs <- markers_Group_2
  DEGs[order(DEGs$logFC,  decreasing = TRUE), c("feature", "group", "logFC")][1:20, c(1)] -> top20;       print(paste0(Group_1," Top 20 DEGs: ", top20))
  print("-----")
  DEGs[order(DEGs$logFC,  decreasing = FALSE), c("feature", "group", "logFC")][1:20, c(1)] -> bottom_20;
  print(paste0(Group_2, " Top 20 DEGs: ", bottom_20))
  
  Total_DEGs <- dim(DEGs[DEGs$logFC >= 0.5 | DEGs$logFC < -0.5, ])[1]
  Up_DEGs <- dim(DEGs[DEGs$logFC >= 0.5, ])[1]
  Down_DEGs <- dim(DEGs[DEGs$logFC < -0.5, ])[1]
  
  # Plot volcano plot 
  library(EnhancedVolcano)
  DEGs <- DEGs[!grepl("MT-", DEGs$feature), ] #removing mitochondrial genes from the plot, but they will still be on the DEGs table
  rownames(DEGs) <- DEGs$feature 
  p <- EnhancedVolcano(
    DEGs,
    lab = DEGs$feature,
    x = 'logFC',
    y = 'padj',
    #selectLab = c(top20[1:6], bottom_20[1:6]),
    selectLab = markers_Diego,
    pCutoff = 10e-2, #pvalue cutoff line
    FCcutoff = 0.5, #foldChange cutoff line
    xlab = paste0("<---- ",Group_1,"   " ,"Log2 FoldChange", "   ",Group_2," ---->"),
    pointSize = 4.0,
    labSize = 6.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    #colAlpha = 4/5,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    max.overlaps = Inf,
    # caption = paste0("Total DEGs: ", Total_DEGs, "\n",
    #                  "Up genes: ", Up_DEGs, " towards ", Group_1, "\n",
    #                  "Up genes: ", Down_DEGs, " towards ", Group_2, "\n",
    #                  "p < 0.05, FC > |0.5|"),
    title = paste0("DEGs ", Group_1, " vs ", Group_2),
    subtitle = "" )
  
  return(list(p, markers_Group_2))
}

All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_processed.qs")
DefaultAssay(All_Sobj) <- "RNA"

comparisons <- list(
  "Comparison" = list(
    group_1 = "Neonate RSV infected (NO IFN)",
    group_2 = "Neonate IFN and RSV infected",
    group_3 = "Adult RSV infected"
  )
)

library(readr)
markers_Diego <- read_csv("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/markers_Diego_on_panel.csv") %>% data.frame()

dendri_genes <- unique(c(trimws(unlist(strsplit(markers_Diego$all_markers[1], ",")))))
epi_genes <- unique(c(trimws(unlist(strsplit(markers_Diego$all_markers[2], ",")))))
lynph_genes <- unique(c(trimws(unlist(strsplit(markers_Diego$all_markers[3], ",")))))


## Dendritic_cell DEGs ----
plt_dend_1_2 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Dendritic_cell")),
                                Identity = "Group", 
                                Group_1 = comparisons$Comparison$group_1, 
                                Group_2 = comparisons$Comparison$group_2, 
                                markers_Diego = dendri_genes)

plt_dend_1_2[[2]][plt_dend_1_2[[2]]$feature %in% dendri_genes, ] #no Diego's dendri_genes 
plt_dend_1_2[[1]]

DEGs_markers <- plt_dend_1_2[[2]]
top_gene_G1 <- DEGs_markers[order(DEGs_markers$logFC, decreasing = FALSE), ]$feature[1:3]
plt_dend_1_2 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Dendritic_cell")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_1, 
                                    Group_2 = comparisons$Comparison$group_2, 
                                    markers_Diego = c(top_gene_G1, dendri_genes))

# if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
plt_dend_1_2[[1]]





plt_dend_1_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Dendritic_cell")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_1, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = dendri_genes)

plt_dend_1_3[[2]][plt_dend_1_3[[2]]$feature %in% dendri_genes, ] #no Diego's dendri_genes 

plt_dend_1_3[[1]]


DEGs_markers <- plt_dend_1_3[[2]]
top_gene_G1 <- DEGs_markers[order(DEGs_markers$logFC, decreasing = FALSE), ]$feature[1:4]
plt_dend_1_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Dendritic_cell")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_1, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = c(top_gene_G1, dendri_genes))
                                    
plt_dend_1_3[[1]]
                                    
                                    
        




plt_dend_2_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Dendritic_cell")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_2, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = dendri_genes)

plt_dend_2_3[[2]][plt_dend_2_3[[2]]$feature %in% dendri_genes, ] #no Diego's dendri_genes 

plt_dend_2_3[[1]]

DEGs_markers <- plt_dend_2_3[[2]]
top_gene_G1 <- DEGs_markers[order(DEGs_markers$logFC, decreasing = FALSE), ]$feature[1:4]
plt_dend_2_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Dendritic_cell")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_2, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = c(top_gene_G1, dendri_genes))
                                    
plt_dend_2_3[[1]]                            
                                    
# all padj <= 0.05      
DEG_1_2 <- plt_dend_1_2[[2]]    
DEG_1_3 <- plt_dend_1_3[[2]]   
DEG_2_3 <- plt_dend_2_3[[2]]

comparisons$Comparison

DEG_list <- list(
  Neo_RSV__Neo_RSV_IFN = DEG_1_2[DEG_1_2$logFC <= -0.5 | DEG_1_2$logFC >= 0.05, ]$feature,
  Neo_RSV__Adu_RSV = DEG_1_3[DEG_1_3$logFC <= -0.5 | DEG_1_3$logFC >= 0.05, ]$feature,
  Neo_RSV_IFN__Adu_RSV = DEG_2_3[DEG_2_3$logFC <= -0.5 | DEG_2_3$logFC >= 0.05, ]$feature
)

library(ggplot2)
library(UpSetR)
upset(fromList(DEG_list), order.by = "freq", nsets = length(DEG_list)) 

x <- upset(fromList(DEG_list), nsets = length(DEG_list))
x$New_data[1:3, 1:3]
dim(x$New_data) #

x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
             x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
             x$New_data$Neo_RSV__Adu_RSV %in% 0, ]

x1 <- unique(unlist(DEG_list, use.names = FALSE))
gene_1 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
                x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

gene_2 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                x$New_data$Neo_RSV_IFN__Adu_RSV %in% 0 &
                x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

RSV_IFN_gene_effects_dendritic <- c(gene_1, gene_2)

RSV_IFN_genes_df <- data.frame(genes = RSV_IFN_gene_effects_dendritic,
                               CellType = "Dendritic",
           description = "genes most likely related to the RSV_IFN effect")

comparisons$Comparison

DEG_1_2 <- plt_dend_1_2[[2]]    
DEG_1_3 <- plt_dend_1_3[[2]]   
DEG_2_3 <- plt_dend_2_3[[2]]


# Filtering top genes to share with Diego
DEG_1_2_top <- rbind(DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% head(n = 20),
      DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_2_top$CellType <- "Dendritic"


DEG_1_3_top <- rbind(DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_3_top$CellType <- "Dendritic"


DEG_2_3_top <- rbind(DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_2_3_top$CellType <- "Dendritic"

library(openxlsx)
sheets_list <- list("RSV_IFN_genes_hit" = RSV_IFN_genes_df, 
                    "NeoRSV_NeoRSVIFN" = DEG_1_2_top, 
                    "NeoRSV_AduRSV" = DEG_1_3_top,
                    "NeoRSVIFN_AduRSV" = DEG_2_3_top)

write.xlsx(sheets_list, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/Dendritic_DEG_AduRSV_NeoRSV_NeoRSVIFN.xlsx")


## Epithelial_cell DEGs ----
plt_dend_1_2 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Epithelial_cell")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_1, 
                                    Group_2 = comparisons$Comparison$group_2, 
                                    markers_Diego = epi_genes)

plt_dend_1_2[[2]][plt_dend_1_2[[2]]$feature %in% epi_genes, ] #no Diego's epi_genes 
plt_dend_1_2[[1]]

DEGs_markers <- plt_dend_1_2[[2]]
top_gene_G1 <- DEGs_markers[order(DEGs_markers$logFC, decreasing = FALSE), ]$feature[1:3]
plt_dend_1_2 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Epithelial_cell")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_1, 
                                    Group_2 = comparisons$Comparison$group_2, 
                                    markers_Diego = c(top_gene_G1, epi_genes))

# if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
plt_dend_1_2[[1]]





plt_dend_1_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Epithelial_cell")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_1, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = epi_genes)

plt_dend_1_3[[2]][plt_dend_1_3[[2]]$feature %in% epi_genes, ] #no Diego's epi_genes 

plt_dend_1_3[[1]]


DEGs_markers <- plt_dend_1_3[[2]]
top_gene_G1 <- DEGs_markers[order(DEGs_markers$logFC, decreasing = FALSE), ]$feature[1:4]
plt_dend_1_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Epithelial_cell")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_1, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = c(top_gene_G1, epi_genes))

plt_dend_1_3[[1]]







plt_dend_2_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Epithelial_cell")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_2, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = epi_genes)

plt_dend_2_3[[2]][plt_dend_2_3[[2]]$feature %in% epi_genes, ] #no Diego's epi_genes 

plt_dend_2_3[[1]]

DEGs_markers <- plt_dend_2_3[[2]]
top_gene_G1 <- DEGs_markers[order(DEGs_markers$logFC, decreasing = FALSE), ]$feature[1:4]
plt_dend_2_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Epithelial_cell")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_2, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = c(top_gene_G1, epi_genes))

plt_dend_2_3[[1]]                            

# all padj <= 0.05      
DEG_1_2 <- plt_dend_1_2[[2]]    
DEG_1_3 <- plt_dend_1_3[[2]]   
DEG_2_3 <- plt_dend_2_3[[2]]

comparisons$Comparison

DEG_list <- list(
  Neo_RSV__Neo_RSV_IFN = DEG_1_2[DEG_1_2$logFC <= -0.5 | DEG_1_2$logFC >= 0.05, ]$feature,
  Neo_RSV__Adu_RSV = DEG_1_3[DEG_1_3$logFC <= -0.5 | DEG_1_3$logFC >= 0.05, ]$feature,
  Neo_RSV_IFN__Adu_RSV = DEG_2_3[DEG_2_3$logFC <= -0.5 | DEG_2_3$logFC >= 0.05, ]$feature
)

library(ggplot2)
library(UpSetR)
upset(fromList(DEG_list), order.by = "freq", nsets = length(DEG_list)) 

x <- upset(fromList(DEG_list), nsets = length(DEG_list))
x$New_data[1:3, 1:3]
dim(x$New_data) #

x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
             x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
             x$New_data$Neo_RSV__Adu_RSV %in% 0, ]

x1 <- unique(unlist(DEG_list, use.names = FALSE))
gene_1 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

gene_2 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 0 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

RSV_IFN_gene_effects_Epithelial <- c(gene_1, gene_2)

RSV_IFN_genes_df <- data.frame(genes = RSV_IFN_gene_effects_Epithelial,
                               CellType = "Epithelial",
                               description = "genes most likely related to the RSV_IFN effect")

comparisons$Comparison

DEG_1_2 <- plt_dend_1_2[[2]]    
DEG_1_3 <- plt_dend_1_3[[2]]   
DEG_2_3 <- plt_dend_2_3[[2]]


# Filtering top genes to share with Diego
DEG_1_2_top <- rbind(DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_2_top$CellType <- "Epithelial"


DEG_1_3_top <- rbind(DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_3_top$CellType <- "Epithelial"


DEG_2_3_top <- rbind(DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_2_3_top$CellType <- "Epithelial"

library(openxlsx)
sheets_list <- list("RSV_IFN_genes_hit" = RSV_IFN_genes_df, 
                    "NeoRSV_NeoRSVIFN" = DEG_1_2_top, 
                    "NeoRSV_AduRSV" = DEG_1_3_top,
                    "NeoRSVIFN_AduRSV" = DEG_2_3_top)

write.xlsx(sheets_list, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/Epithelial_DEG_AduRSV_NeoRSV_NeoRSVIFN.xlsx")



## Lymphocyte DEGs ----
plt_dend_1_2 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Lymphocyte")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_1, 
                                    Group_2 = comparisons$Comparison$group_2, 
                                    markers_Diego = lynph_genes)

plt_dend_1_2[[2]][plt_dend_1_2[[2]]$feature %in% lynph_genes, ] #no Diego's lynph_genes 
plt_dend_1_2[[1]]

DEGs_markers <- plt_dend_1_2[[2]]
top_gene_G1 <- DEGs_markers[order(DEGs_markers$logFC, decreasing = FALSE), ]$feature[1:3]
plt_dend_1_2 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Lymphocyte")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_1, 
                                    Group_2 = comparisons$Comparison$group_2, 
                                    markers_Diego = c(top_gene_G1, lynph_genes))

# if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
plt_dend_1_2[[1]]





plt_dend_1_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Lymphocyte")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_1, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = lynph_genes)

plt_dend_1_3[[2]][plt_dend_1_3[[2]]$feature %in% lynph_genes, ] #no Diego's lynph_genes 

plt_dend_1_3[[1]]


DEGs_markers <- plt_dend_1_3[[2]]
top_gene_G1 <- DEGs_markers[order(DEGs_markers$logFC, decreasing = FALSE), ]$feature[1:4]
plt_dend_1_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Lymphocyte")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_1, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = c(top_gene_G1, lynph_genes))

plt_dend_1_3[[1]]







plt_dend_2_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Lymphocyte")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_2, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = lynph_genes)

plt_dend_2_3[[2]][plt_dend_2_3[[2]]$feature %in% lynph_genes, ] #no Diego's lynph_genes 

plt_dend_2_3[[1]]

DEGs_markers <- plt_dend_2_3[[2]]
top_gene_G1 <- DEGs_markers[order(DEGs_markers$logFC, decreasing = FALSE), ]$feature[1:4]
plt_dend_2_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Lymphocyte")),
                                    Identity = "Group", 
                                    Group_1 = comparisons$Comparison$group_2, 
                                    Group_2 = comparisons$Comparison$group_3, 
                                    markers_Diego = c(top_gene_G1, lynph_genes))

plt_dend_2_3[[1]]                            

# all padj <= 0.05      
DEG_1_2 <- plt_dend_1_2[[2]]    
DEG_1_3 <- plt_dend_1_3[[2]]   
DEG_2_3 <- plt_dend_2_3[[2]]

comparisons$Comparison

DEG_list <- list(
  Neo_RSV__Neo_RSV_IFN = DEG_1_2[DEG_1_2$logFC <= -0.5 | DEG_1_2$logFC >= 0.05, ]$feature,
  Neo_RSV__Adu_RSV = DEG_1_3[DEG_1_3$logFC <= -0.5 | DEG_1_3$logFC >= 0.05, ]$feature,
  Neo_RSV_IFN__Adu_RSV = DEG_2_3[DEG_2_3$logFC <= -0.5 | DEG_2_3$logFC >= 0.05, ]$feature
)

library(ggplot2)
library(UpSetR)
upset(fromList(DEG_list), order.by = "freq", nsets = length(DEG_list)) 

x <- upset(fromList(DEG_list), nsets = length(DEG_list))
x$New_data[1:3, 1:3]
dim(x$New_data) #

x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
             x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
             x$New_data$Neo_RSV__Adu_RSV %in% 0, ]

x1 <- unique(unlist(DEG_list, use.names = FALSE))
gene_1 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

gene_2 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 0 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

RSV_IFN_gene_effects_Lymphocyte <- c(gene_1, gene_2)

RSV_IFN_genes_df <- data.frame(genes = RSV_IFN_gene_effects_Lymphocyte,
                               CellType = "Lymphocyte",
                               description = "genes most likely related to the RSV_IFN effect")

comparisons$Comparison

DEG_1_2 <- plt_dend_1_2[[2]]    
DEG_1_3 <- plt_dend_1_3[[2]]   
DEG_2_3 <- plt_dend_2_3[[2]]


# Filtering top genes to share with Diego
DEG_1_2_top <- rbind(DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_2_top$CellType <- "Lymphocyte"


DEG_1_3_top <- rbind(DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_3_top$CellType <- "Lymphocyte"


DEG_2_3_top <- rbind(DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_2_3_top$CellType <- "Lymphocyte"

library(openxlsx)
sheets_list <- list("RSV_IFN_genes_hit" = RSV_IFN_genes_df, 
                    "NeoRSV_NeoRSVIFN" = DEG_1_2_top, 
                    "NeoRSV_AduRSV" = DEG_1_3_top,
                    "NeoRSVIFN_AduRSV" = DEG_2_3_top)

write.xlsx(sheets_list, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/Lymphocyte_DEG_AduRSV_NeoRSV_NeoRSVIFN.xlsx")













