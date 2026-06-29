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

## Load data -------
obj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/Sobj_list_TMA_ID_fixed.qs")


## Load cell type reference -------
## Cell type reference - profile mtx from CosMx 
# at https://github.com/Nanostring-Biostats/CellProfileLibrary/blob/master/Mouse/Adult/Lung_MCA.RData
load("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/Lung_MCA.RData")
metadata
cellGroups
profile_matrix[1:10, 1:10] # row: mouse genes; col: cell type; content: cell type "fraction" by gene
dim(profile_matrix)

rownames(profile_matrix)


## Run cell type -------
# Define the function to process each Sobj
call_celltype <- function(Sobj, reference_profiles = profile_matrix) {
  # Automatically detect which assay to use
  preferred_assays <- c("Nanostring", "RNA")
  found_assay <- NULL
  for (a in preferred_assays) {
    if (a %in% names(Sobj@assays)) {
      found_assay <- a
      break
    }
  }
  if (is.null(found_assay)) {
    warning("Neither 'Nanostring' nor 'RNA' assay found in Sobj.")
    return(Sobj)
  }
  message("Using assay: ", found_assay)
  
  # Use GetAssayData for Seurat v5+ compatibility
  x_mat <- t(as.matrix(GetAssayData(Sobj, assay = found_assay, layer = "counts")))
  
  # Run insitutypeML
  sup <- insitutypeML(
    x = x_mat,
    neg = Sobj@meta.data$nCount_negprobes,
    # cohort = cohort, # Uncomment if you have a cohort variable
    reference_profiles = as.matrix(reference_profiles)
  )
  
  # Assign clusters
  Sobj$insitutype_cluster <- "other"
  Sobj$insitutype_cluster[names(sup$clust)] <- sup$clust
  names(Sobj@meta.data)[names(Sobj@meta.data) == "insitutype_cluster"] <- "Inst_Celltype"
  
  # # Plot flightpath_plot
  # cols <- brewer.pal(10, 'Paired')
  # cols <- cols[seq_along(unique(sup$clust))]
  # names(cols) <- unique(sup$clust)
  # fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = sup, col = cols[sup$clust])
  # print(fp)
  
  message(Sobj@meta.data$Samples[1], " has been cell typed")
  
  
  return(Sobj)
}


# Initialize an empty list to store metadata
meta_list <- list()

for (i in seq_along(obj_list)) {
  Sobj <- obj_list[[i]]
  names(Sobj@meta.data)[names(Sobj@meta.data) == "Inst_Celltype"] <- "OLD_Inst_Celltype"
  names(Sobj@meta.data)[names(Sobj@meta.data) == "major_celltype"] <- "OLD_major_celltype"
  
  Sobj <- call_celltype(Sobj, reference_profiles = profile_matrix)
  
  # Store metadata in the list, using a sample-specific name if desired
  sample_name <- Sobj@meta.data$TMA_4[1]
  meta_list[[sample_name]] <- Sobj@meta.data
  
  # Save the list after each iteration
  qsave(meta_list, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/metadata_list.qs")
  
  message("Saved metadata for ", sample_name)
}





# # # Resume Processing After a Crash
# meta_list <- if (file.exists("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/metadata_list.qs")) {
#   qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/metadata_list.qs")
# } else {
#   list()
# }
# 
# for (i in seq_along(obj_list)) {
#   Sobj <- obj_list[[i]]
#   sample_name <- Sobj@meta.data$TMA_3[1]
# 
#   # Skip if already processed
#   if (sample_name %in% names(meta_list)) {
#     message("Skipping already processed sample: ", sample_name)
#     next
#   }
# 
#   Sobj <- call_celltype(Sobj, reference_profiles = profile_matrix)
#   meta_list[[sample_name]] <- Sobj@meta.data
#   qsave(meta_list, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/metadata_list.qs")
#   message("Saved metadata for ", sample_name)
# }






## Reading metadata back to the obj list
library(qs)
meta_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/metadata_list.qs")
for (i in seq_along(obj_list)) {
  sample_name <- obj_list[[i]]@meta.data$TMA_3[1]
  if (sample_name %in% names(meta_list)) {
    obj_list[[i]]@meta.data <- meta_list[[sample_name]]
  }
}


all_meta <- lapply(obj_list, function(obj) {
  meta <- obj@meta.data
  return(meta)
})
all_meta <- do.call(plyr::rbind.fill, all_meta)


table(all_meta$Inst_Celltype) %>% names()
table(all_meta$OLD_major_celltype) %>% names()

major_celltype_dict <- c(
  "Alveolar.bipotent.progenitor" = "Alveolar_bipotent_progenitor",
  "Alveolar.macrophage.Ear2.high" = "Alveolar_macrophage",
  "Alveolar.macrophage.Pclaf.high" = "Alveolar_macrophage",
  "AT1.cell" = "AT1_cell",
  "AT2.cell" = "AT2_cell",
  "B.cell" = "Lymphocyte",
  "IgA.producing.B.cell" = "Lymphocyte",
  "Basophil" = "Granulocyte",
  "Ciliated.cell" = "Epithelial_cell",
  "Clara.cell" = "Epithelial_cell",
  "Dendritic.cell.Gngt2.high" = "Dendritic_cell",
  "Dendritic.cell.H2.M2.high" = "Dendritic_cell",
  "Dendritic.cell.Mgl2.high" = "Dendritic_cell",
  "Dendritic.cell.Naaa.high" = "Dendritic_cell",
  "Dendritic.cell.Tubb5.high" = "Dendritic_cell",
  "Endothelial.cell.Kdr.high" = "Endothelial_cell",
  "Endothelial.cell.Tmem100.high" = "Endothelial_cell",
  "Endothelial.cell.Vwf.high" = "Endothelial_cell",
  "Eosinophil" = "Granulocyte",
  "Interstitial.macrophage" = "Interstitial_macrophage",
  "Monocyte.progenitor" = "Monocyte_progenitor",
  "Neutrophil" = "Granulocyte",
  "NK.cell" = "Lymphocyte",
  "Nuocyte" = "Nuocyte",
  "Plasmacytoid.dendritic.cell" = "Dendritic_cell",
  "Stromal.cell.Acta2.high" = "Stromal_cell",
  "Stromal.cell.Dcn.high" = "Stromal_cell",
  "Stromal.cell.Inmt.high" = "Stromal_cell",
  "T.cell.Cd8b1.high" = "Lymphocyte"
)

## Get major_celltype
for (i in seq_along(obj_list)) {
  obj <- obj_list[[i]]
  # Map minor to major cell type using your dictionary 
  obj@meta.data$major_celltype <- major_celltype_dict[obj@meta.data$Inst_Celltype]
  obj_list[[i]] <- obj
}


all_meta <- lapply(obj_list, function(obj) {
  meta <- obj@meta.data
  return(meta)
})
all_meta <- do.call(plyr::rbind.fill, all_meta)

test <- all_meta[all_meta$TMA_3 %in% "Ap17-0333 A1 1", ]
pheatmap(table(test$OLD_major_celltype,
               test$major_celltype))

pheatmap(table(all_meta$OLD_major_celltype,
               all_meta$major_celltype))

names(table(all_meta$OLD_major_celltype))
names(table(all_meta$major_celltype))

## Comparing the old and the current cell anno: 
# Main diff - old cell anno was done by Slide; current cell anno was done by TMA
library(pheatmap)
# create contingency table of major cell types (after cleaning)
ct_table <- table(all_meta$OLD_major_celltype, all_meta$major_celltype)

# Optional: aggregate rare categories or filter top frequent ones before plotting

# Scale rows (cell types) to show relative proportions within each OLD_major_celltype
mat_scaled <- t(apply(ct_table, 1, function(x) x / sum(x)))

# plot with better colors and clustering
pheatmap(mat_scaled,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("white", "blue"))(50),
         fontsize_row = 10,
         fontsize_col = 10,
         angle_col = 45)


## Set color map 
cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- all_meta$major_celltype %>% table() %>% names()
cols <- setNames(cols[1:length(factor_levels)], factor_levels)


library(ggplot2)
ggplot(all_meta, aes(x = TMA_3, fill = major_celltype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols) +
  labs(
    x = "",
    y = "",
    fill = "",
    title = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      face = "bold", size = 16, #hjust = 0, 
      color = "black"
    ),
    axis.title.x = element_text(
      face = "bold", size = 12, #color = "black"
    ),
    axis.title.y = element_text(
      face = "bold", size = 14, color = "black"
    ),
    axis.text.x = element_text(
      face = "plain", size = 10, angle = 0, #hjust = 0, 
      color = "black"
    ),
    axis.text.y = element_text(
      face = "plain", size = 12, color = "black"
    ),
    legend.title = element_text(
      face = "bold", size = 14, color = "black"
    ),
    legend.text = element_text(
      face = "plain", size = 12, color = "black"
    )
  )



library(ggplot2)
ggplot(all_meta, aes(x = Group, fill = major_celltype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols) +
  labs(
    x = "",
    y = "",
    fill = "",
    title = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      face = "bold", size = 16, #hjust = 0, 
      color = "black"
    ),
    axis.title.x = element_text(
      face = "bold", size = 12, #color = "black"
    ),
    axis.title.y = element_text(
      face = "bold", size = 14, color = "black"
    ),
    axis.text.x = element_text(
      face = "plain", size = 10, angle = 0, #hjust = 0, 
      color = "black"
    ),
    axis.text.y = element_text(
      face = "plain", size = 12, color = "black"
    ),
    legend.title = element_text(
      face = "bold", size = 14, color = "black"
    ),
    legend.text = element_text(
      face = "plain", size = 12, color = "black"
    )
  )









