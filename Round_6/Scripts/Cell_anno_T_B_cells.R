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



## Load data ---------
obj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/Sobj_list_TMA_ID_fixed.qs")

# All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/ScObj.qs")

ScObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/ScObj_umap.qs")

# plot markers dotplot
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


All_Sobj <- NormalizeData(All_Sobj, assay = "Nanostring")
DefaultAssay(All_Sobj) <- "Nanostring"

cd4_markers <- intersect(cell_types$Lymphocytes$T_CD4_Cells$Markers, All_Sobj@assays$Nanostring$data %>% rownames())
cd8_markers <- intersect(cell_types$Lymphocytes$T_CD8_Cells$Markers, All_Sobj@assays$Nanostring$data %>% rownames())
bcell_markers <- intersect(cell_types$Lymphocytes$Activated_B_Cells$Markers, All_Sobj@assays$Nanostring$data %>% rownames())

Idents(All_Sobj) <- "major_celltype"
DotPlot(All_Sobj,
        features =  unique(c(cd4_markers, cd8_markers, bcell_markers))
        #features =  unique(c(top_markers_auc))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("")



## We need to get all the anno. from the individuals TMAs 
obj_list

lyn_markers <- list(cd4_markers = cd4_markers,
                    cd8_markers = cd8_markers,
                    bcell_markers = bcell_markers)

# Run sctype
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
# install.packages('openxlsx')
library(openxlsx)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

  
gs_list <- list(
  gs_positive = lyn_markers
)



seurat_obj <- obj_list[[1]] 
seurat_obj <- NormalizeData(seurat_obj, assay = "Nanostring")
# # get cell-type by cell matrix
# es.max = sctype_score(scRNAseqData = seurat_obj@assays$Nanostring$data %>% as.matrix(),
#                       scaled = TRUE,
#                       gs = gs_list$gs_positive,
#                       gs2 = NULL)
# 
# 
# # get cell-type by cell matrix
# expression_matrix <- GetAssayData(object = seurat_obj, assay = "RNA", slot = "data") %>% as.matrix()


# # I had to run this chunk by chunk to get to the "es.max" 
# 
# es.max = sctype_score(scRNAseqData = seurat_obj@assays$RNA$data %>% as.matrix(), scaled = FALSE, 
#                       gs = gs_list$gs_positive, gs2 = NULL) 

scRNAseqData <- seurat_obj@assays$Nanostring$data %>% as.matrix()
scaled = TRUE
# gs = gs_list$gs_positive
gs = gs_list$gs_positive[2:3]
gs2 = NULL

# function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){

# check input matrix
if(!is.matrix(scRNAseqData)){
  warning("scRNAseqData doesn't seem to be a matrix")
} else {
  if(sum(dim(scRNAseqData))==0){
    warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
  }
}

# marker sensitivity
marker_stat = sort(table(unlist(gs)), decreasing = T); 
marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                gene_ = names(marker_stat), stringsAsFactors = !1)

# convert gene names to Uppercase
gene_names_to_uppercase <- F
if(gene_names_to_uppercase){
  rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
}

# subselect genes only found in data
names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
gs = lapply(1:length(gs), function(d_){ 
  GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
gs2 = lapply(1:length(gs2), function(d_){ 
  GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]

# z-scale if not
if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData

# multiple by marker sensitivity
for(jj in 1:nrow(cell_markers_genes_score)){
  Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
}

# subselect only with marker genes
Z = Z[unique(c(unlist(gs),unlist(gs2))), ]

# combine scores
es = do.call("rbind", lapply(names(gs), function(gss_){ 
  sapply(1:ncol(Z), function(j) {
    gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
    sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
    if(is.na(sum_t2)){
      sum_t2 = 0;
    }
    sum_t1 + sum_t2
  })
})) 

dimnames(es) = list(names(gs), colnames(Z))
es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows

es.max




# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seurat_obj@meta.data$major_celltype), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$major_celltype==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$major_celltype==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  


# set low-confident (low ScType score) clusters to "unknown"
## gs2 = NULL when there is no negative marker
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < 
                     sctype_scores$ncells / 5] <- "Unsigned"

table(sctype_scores$type)

print(sctype_scores[,1:3])

seurat_obj@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_obj@meta.data$customclassif[seurat_obj@meta.data$major_celltype == j] = as.character(cl_type$type[1])
}

DefaultAssay(seurat_obj) <- "Nanostring"
Idents(seurat_obj) <- "major_celltype"
DotPlot(seurat_obj,
        features =  unique(c(cd4_markers, cd8_markers, bcell_markers))
        #features =  unique(c(top_markers_auc))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("")

Idents(seurat_obj) <- "customclassif"
DotPlot(seurat_obj,
        features =  unique(c(cd4_markers, cd8_markers, bcell_markers))
        #features =  unique(c(top_markers_auc))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("")

table(seurat_obj$customclassif)






sctype_score_custom <- function(
    seurat_obj,
    assay = "Nanostring",
    gs, 
    gs2 = NULL,
    scaled = TRUE,
    cluster_anno = "major_celltype",
    confidence_factor = 5,
    new_anno_name = "customclassif"
){
  # extract expression matrix
  scRNAseqData <- seurat_obj@assays[[assay]]$data %>% as.matrix()
  
  # check input matrix
  if (!is.matrix(scRNAseqData)) {
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if (sum(dim(scRNAseqData)) == 0) {
      warning("The dimension of input scRNAseqData matrix equals 0, is it empty?")
    }
  }
  
  # marker sensitivity
  marker_stat <- sort(table(unlist(gs)), decreasing = TRUE)
  marker_sensitivity <- data.frame(
    score_marker_sensitivity = scales::rescale(
      as.numeric(marker_stat), to = c(0, 1), from = c(length(gs), 1)
    ),
    gene_ = names(marker_stat),
    stringsAsFactors = FALSE
  )
  
  # subselect genes only found in data
  names_gs_cp <- names(gs)
  names_gs_2_cp <- names(gs2)
  gs <- lapply(1:length(gs), function(d_) {
    GeneIndToKeep <- rownames(scRNAseqData) %in% as.character(gs[[d_]])
    rownames(scRNAseqData)[GeneIndToKeep]
  })
  gs2 <- lapply(1:length(gs2), function(d_) {
    GeneIndToKeep <- rownames(scRNAseqData) %in% as.character(gs2[[d_]])
    rownames(scRNAseqData)[GeneIndToKeep]
  })
  names(gs) <- names_gs_cp
  names(gs2) <- names_gs_2_cp
  
  cell_markers_genes_score <- marker_sensitivity[
    marker_sensitivity$gene_ %in% unique(unlist(gs)), ]
  
  # scale if needed
  if (!scaled) {
    Z <- t(scale(t(scRNAseqData)))
  } else {
    Z <- scRNAseqData
  }
  
  # multiply scores by sensitivity
  for (jj in 1:nrow(cell_markers_genes_score)) {
    Z[cell_markers_genes_score[jj, "gene_"], ] <- 
      Z[cell_markers_genes_score[jj, "gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect marker genes
  Z <- Z[unique(c(unlist(gs), unlist(gs2))), ]
  
  # combine scores
  es <- do.call("rbind", lapply(names(gs), function(gss_) {
    sapply(1:ncol(Z), function(j) {
      gs_z <- Z[gs[[gss_]], j]
      gz_2 <- Z[gs2[[gss_]], j] * -1
      sum_t1 <- (sum(gs_z) / sqrt(length(gs_z)))
      sum_t2 <- sum(gz_2) / sqrt(length(gz_2))
      if (is.na(sum_t2)) {
        sum_t2 <- 0
      }
      sum_t1 + sum_t2
    })
  }))
  
  dimnames(es) <- list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all), ]
  
  # merge by cluster
  library(dplyr)
  if (length(gs) < 2) {
    cL_results <- unique(seurat_obj@meta.data[[cluster_anno]]) %>%
      lapply(function(cl) {
        cells_in_cluster <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[cluster_anno]] == cl, ])
        
        # sum scores only for these cells, not entire matrix
        es.max <- data.frame(es.max)
        cluster_scores <- es.max$es.max %>% sum()
        
        cluster_scores_sorted <- sort(cluster_scores, decreasing = TRUE)
        
        data.frame(
          cluster = cl,
          type = names(gs),
          scores = cluster_scores_sorted,
          ncells = length(cells_in_cluster)
        ) %>% head(10)
      }) %>%
      bind_rows()
  
  } else {
    cL_results <- unique(seurat_obj@meta.data[[cluster_anno]]) %>%
      lapply(function(cl) {
        cells_in_cluster <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[cluster_anno]] == cl, ])
        cluster_scores <- rowSums(es.max[, cells_in_cluster])
        cluster_scores_sorted <- sort(cluster_scores, decreasing = TRUE)
        data.frame(
          cluster = cl,
          type = names(cluster_scores_sorted),
          scores = cluster_scores_sorted,
          ncells = length(cells_in_cluster)
        ) %>% head(10)
      }) %>%
      bind_rows()
  }
  
  sctype_scores <- cL_results %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = scores)
  
  # low-confidence clusters to "Unsigned"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) <
                       sctype_scores$ncells / confidence_factor] <- "Unsigned"
  
  # assign back to Seurat object with custom annotation column name
  seurat_obj@meta.data[[new_anno_name]] <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, ]
    seurat_obj@meta.data[[new_anno_name]][seurat_obj@meta.data[[cluster_anno]] == j] <- as.character(cl_type$type[1])
  }
  
  # return list of seurat object and scores dataframe
  return(list(
    seurat_obj = seurat_obj,
    sctype_scores = sctype_scores
  ))
}



seurat_obj <- obj_list[[1]] 
seurat_obj <- NormalizeData(seurat_obj, assay = "Nanostring")
out <- sctype_score_custom(
    seurat_obj,
    assay = "Nanostring",
    gs = gs_list$gs_positive, 
    gs2 = NULL,
    scaled = TRUE,
    cluster_anno = "major_celltype",
    confidence_factor = 4, # 2 conservative   
    new_anno_name = "customclassif"
)
  
# out$sctype_scores
# seurat_obj <- out$seurat_obj
 
Idents(seurat_obj) <- "customclassif"
DotPlot(seurat_obj,
        features =  unique(c(cd4_markers, cd8_markers, bcell_markers))
        #features =  unique(c(top_markers_auc))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("")



### IT IS NOT WORKING FOR ALL THE TMAS FOR SOME REASON .....
library(Seurat)
library(ggplot2)

pdf("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Plots/tcell_bcell_markers_dotplot_by_TMA.pdf", width=10, height=7)  # Open one PDF for output

for (i in seq_along(obj_list)) {
  seurat_obj <- obj_list[[i]]
  seurat_obj <- NormalizeData(seurat_obj, assay = "Nanostring")
  out <- sctype_score_custom(
    seurat_obj,
    assay = "Nanostring",
    gs = gs_list$gs_positive, 
    gs2 = NULL,
    scaled = TRUE,
    cluster_anno = "major_celltype",
    confidence_factor = 4,
    new_anno_name = "customclassif"
  )
  seurat_obj <- out$seurat_obj
  Idents(seurat_obj) <- "customclassif"
  
  plot <- DotPlot(
    seurat_obj,
    features = unique(c(cd4_markers, cd8_markers, bcell_markers))
  ) +
    geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
    scale_colour_viridis(option = "magma") +
    guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste("DotPlot:", names(obj_list)[i]))
  
  print(plot)  # Print each plot to the PDF
}

dev.off()  # Close the PDF device



### Let's work off the single merged object then 

# Make sure Seurat and dplyr are loaded
library(Seurat)
library(dplyr)

# Extract raw counts matrix from Nanostring assay
counts <- ScObj@assays$Nanostring$counts

# Initialize annotation vector with default "other"
annot <- rep("other", ncol(counts))

# Hierarchical annotation based on transcript counts
# Check Ptprc > 0
has_Ptprc <- counts["Ptprc", ] > 0
# Among Ptprc positive, check Cd3e > 0
has_Cd3e <- counts["Cd3e", ] > 0
# Among Cd3e positive, check Cd4 and Cd8a
has_Cd4 <- counts["Cd4", ] > 0
has_Cd8a <- counts["Cd8a", ] > 0
# Check B cell markers: Ptprc and Cd19 > 0
has_Cd19 <- counts["Cd19", ] > 0

# Assign Tcd4
annot[has_Ptprc & has_Cd3e & has_Cd4] <- "Tcd4"
# Assign Tcd8
annot[has_Ptprc & has_Cd3e & has_Cd8a] <- "Tcd8"
# Assign B cell
annot[has_Ptprc & !has_Cd3e & has_Cd19] <- "Bcell"

# Remaining cells retain "other"

# Add annotation to Seurat meta data
ScObj@meta.data$lyn_minor_celltype <- annot
table(ScObj@meta.data$lyn_minor_celltype)

# now plot it by TMAs 
ScObj <- NormalizeData(ScObj, assay = "Nanostring")
DefaultAssay(ScObj) <- "Nanostring"
Idents(ScObj) <- "lyn_minor_celltype"
DotPlot(ScObj,
        features =  unique(c(cd4_markers, cd8_markers, bcell_markers))
        #features =  unique(c(top_markers_auc))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("")

dittoBarPlot(
  subset(ScObj, lyn_minor_celltype %in% "other", invert = TRUE),
  var = "lyn_minor_celltype",
  group.by = "TMA_4",
  #color.panel = minor_celltype_concate_colors,
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)


dittoBarPlot(
  subset(ScObj, lyn_minor_celltype %in% "other", invert = TRUE),
  var = "lyn_minor_celltype",
  group.by = "Group",
  #color.panel = minor_celltype_concate_colors,
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)

dittoBarPlot(
  ScObj,
  #subset(ScObj, lyn_minor_celltype %in% "other", invert = TRUE),
  var = "lyn_minor_celltype",
  group.by = "major_celltype",
  #color.panel = minor_celltype_concate_colors,
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)

Objsub <- subset(ScObj, lyn_minor_celltype %in% "other", invert = TRUE)
pheatmap(table(Objsub$lyn_minor_celltype,Objsub$major_celltype) )


ScObj_sub <- subset(ScObj, lyn_minor_celltype %in% "other", invert = TRUE)
dittoBarPlot(
  subset(ScObj_sub, Group %in% c("Adult RSV infected", "Adult Control")),
  var = "lyn_minor_celltype",
  group.by = "Group",
  color.panel = c("#F0E442", "#E69F00", "#D55E00"),
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)


dittoBarPlot(
  subset(ScObj_sub, Group %in% c("Neonate control", "Neonate RSV infected (NO IFN)", "Neonate IFN and RSV infected")),
  var = "lyn_minor_celltype",
  group.by = "Group",
  color.panel = c("#F0E442", "#E69F00", "#D55E00"),
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)

Idents(ScObj_sub) <- "lyn_minor_celltype"
DimPlot(ScObj_sub)

# NOOo so crap UMAP =[
Idents(ScObj) <- "Inst_Celltype"
DimPlot(ScObj)


table(ScObj$Group)


### Using the single merged obj we've been using -----------
Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_lym_Neocy_Granu_to_minor.qs") 
Sobj_bcells <- subset(Sobj, lym_nuo_granu_recluster %in% "Bcell")

Sobj@meta.data$cell
obj <- obj_list$`Ap17-0333 A1 1`
obj$cell

intersect(obj$cell,
          Sobj@meta.data$cell)


obj@meta.data$Bcell <- "other"
obj@meta.data[obj@meta.data$cell %in% Sobj_bcells$cell, ]$Bcell <- "Bcell"

ImageDimPlot(obj,
             group.by = 'Bcell',
             border.color = "white",
             fov = "FOV",
             border.size = 0.1,
             cols = c("gray", "gold"),
             # molecules = rownames(obj_comp1_1[1:4]),
             # mols.size = 0.5,
             # mols.cols = c("red", "blue", "yellow", "green"),
             # alpha = 0.5,
             size = 1.5,
             axes = F) +
  ggtitle(paste0(names(obj_list)[1]))



obj_list <- obj_list[!grepl("NA", names(obj_list))]

pdf("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Plots/bcell_over_fix_TMA_anno.pdf", onefile = TRUE)


for (i in seq_along(obj_list)) {
  obj <- obj_list[[i]]
  # Check if object is empty (e.g. has no cells or relevant slots)
  if (length(obj@meta.data$cell) == 0) next
  
  obj@meta.data$Bcell <- "other"
  obj@meta.data[obj@meta.data$cell %in% Sobj_bcells$cell, ]$Bcell <- "Bcell"
  
  print(
    ImageDimPlot(obj,
                 group.by = 'Bcell',
                 border.color = "white",
                 fov = "FOV",
                 border.size = 0.1,
                 cols = c("black", "gold"),
                 size = 1.5,
                 axes = FALSE) +
      ggtitle(paste0(names(obj_list)[i], " : ",as.vector(obj_list[[i]]$Group[1])))
  )
}

dev.off()

grep("0332 A1 1.1", names(obj_list))
names(obj_list)[40]
grep("Ap17-0332 A1 1.1", names(obj_list))







# collapse major cell types if necessary
# subset T and B cells pop.
# run seurat workflow 
# keep cell_ID and cell type 
# put it back into ScObj and obj_list






