## ====== Date: Jun 11, 2025  ===== 
## ====== Date: Jun 18, 2025  ===== 
## ====== Date: Jul 03, 2025  ===== 

# ========================================================
## Directions: 
## DEG Analysis in Re-clustered Cell Types: Neonate NO IFN infected vs Neonate IFN infected vs Adult infected
## Re-cluster separated MINOR cell types 
# Dendritic (D1, D2 ...)
# Epithelial (etc ...)
# Lymphocytes (Bcells, Tcells)
# Based on Diego's markers 
## Get Volcano plots for their DEGs
# Eg: Group1_dendritic_01 vs Group2_dendritic_01 vs Group3_dendritic_01
## Workflow 
# 1. Filter major cell types 
# 2. Run through institutype clustering (2, 5, 8, 11, 15)
# 3. Run clustree (not necessary now - I already know the number of cell types/clusters I need to find)
# 4. Find main minor cell types (eg Lymp should have Bcell and Tcells separated before any major dig in)
# Note: No need to go for "more granular cell types" than Diego's markers can tell
# 5. Once we have found those main minor cell types, go for DEG the same way we did for the major cell types
# 6. Export DEGs in a spreadsheet. Share plots and DEGs with Diego's team.
# ========================================================


### Load and define functions ---------
library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(qs)


call_insitu_cluster <- function(Sobj, n_of_cluster, assay, return_flightpath = FALSE) {
  # Check if assay exists
  if (!assay %in% names(Sobj@assays)) {
    message(sprintf("Assay '%s' not found in the Seurat object. Returning input object unchanged.", assay))
    if (return_flightpath) {
      return(list(Sobj = Sobj, flightpath = NULL))
    } else {
      return(Sobj)
    }
  }
  
  # Store flightpaths if requested
  flightpaths <- list()
  
  for (n in n_of_cluster) {
    # Run insitutypeML
    unsup <- insitutype(
      x = as.matrix(t(as.matrix(Sobj@assays[[assay]]$counts))),
      neg = Sobj@meta.data$nCount_negprobes,
      reference_profiles = NULL,
      bg = NULL,
      n_clusts = n,
      n_phase1 = 200,
      n_phase2 = 500,
      n_phase3 = 2000,
      n_starts = 1,
      max_iters = 5
    )
    
    # Create new column for this cluster number, default 'other'
    cluster_col <- paste0("Inst_Cluster_n", n)
    Sobj@meta.data[[cluster_col]] <- "other"
    if(identical(rownames(Sobj@meta.data), names(unsup$clust))) {
      Sobj@meta.data[[cluster_col]] <- unsup$clust
    } else {
      print("We can't add cluster label into Sobj@meta.data that way, bro...")
    }
    
    # Plot flightpath if requested
    if (return_flightpath) {
      cols = RColorBrewer::brewer.pal(10, 'Paired')
      cols <- cols[seq_along(unique(unsup$clust))]
      names(cols) <- unique(unsup$clust)
      fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup, col = cols[unsup$clust])
      print(fp)
      flightpaths[[cluster_col]] <- fp
    }
  }
  
  if (return_flightpath) {
    return(list(Sobj = Sobj, flightpaths = flightpaths))
  } else {
    return(Sobj)
  }
}

Plot_Volcano_DEG_sc <- function(
    SeuratObj = NULL, # seurat object
    Identity = NULL, # group column name
    Group_1 = NULL, # group to compare
    Group_2 = NULL, # group to compare
    logFC_direction = NULL,    # base line comparison (eg: control)
    markers_to_label = NULL,   # gene names on vulcan plot
    filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
    filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
    filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
    remove_mt = FALSE,
    genes_to_remove_from_vis = NULL# TRUE to remove mitochondrial genes
) {
  # Subset to comparison 
  Idents(SeuratObj) <- Identity
  SeuratObj_sub <- subset(x = SeuratObj, idents = c(Group_1, Group_2), invert = FALSE)
  
  # Differential Expression Analysis
  library(presto)
  markers <- presto::wilcoxauc(SeuratObj_sub, group_by = Identity)
  
  # --- Filtering Section ---
  if (!is.null(filter_logFC)) {
    markers <- markers[abs(markers$logFC) >= filter_logFC, ]
  }
  if (!is.null(filter_padj)) {
    markers <- markers[markers$padj < filter_padj, ]
  }
  if (!is.null(filter_auc)) {
    markers <- markers[markers$auc > filter_auc, ]
  }
  if (remove_mt) {
    markers <- markers[!grepl("^MT-", markers$feature), ]
  }
  # -------------------------
  
  if (!is.null(genes_to_remove_from_vis)) {
    pattern <- paste(genes_to_remove_from_vis, collapse = "|")
    markers <- markers[!grepl(pattern, markers$feature), ]
  }
  
  # Only keep markers for specified group(s)
  if (is.null(logFC_direction)) {
    logFC_direction <- Group_2   # default to Group_2 if not specified
  }
  markers <- markers[!markers$group %in% logFC_direction, ]
  markers <- markers[order(markers$logFC, decreasing = TRUE), ]
  
  # Define top and bottom DEGs
  DEGs <- markers
  top20 <- DEGs[order(DEGs$logFC, decreasing = TRUE), "feature"][1:20]
  bottom_20 <- DEGs[order(DEGs$logFC, decreasing = FALSE), "feature"][1:20]
  
  if (!is.null(markers_to_label)) {
    markers_to_label <- c(top20[1:5], bottom_20[1:5], markers_to_label)
  } else {
    markers_to_label <- c(top20[1:5], bottom_20[1:5])
  }
  
  # Plot volcano plot 
  library(EnhancedVolcano)
  rownames(DEGs) <- DEGs$feature 
  p <- EnhancedVolcano(
    DEGs,
    lab = DEGs$feature,
    x = 'logFC',
    y = 'padj',
    selectLab = markers_to_label,
    pCutoff = 10e-2,
    FCcutoff = 0.5,
    xlab = paste0("<---- ",Group_1,"   Log2 FoldChange   ",Group_2," ---->"),
    pointSize = 4.0,
    labSize = 6.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    max.overlaps = Inf,
    caption = "",
    title = paste0("DEGs ", Group_1, " vs ", Group_2),
    subtitle = "" )
  
  return(list(p, DEGs))
}



### Load and define variables  ---------
Sobj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj_list_celltype_niche.qs")

All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_processed.qs")

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


### Lymphocytes re-clustering -------

# 1. Define minor cell types
## Keep target cell type 
Sobj <- subset(All_Sobj, major_celltype %in% c("Granulocyte",
                                                   "Nuocyte",
                                                   "Lymphocyte") )
Sobj@assays$RNA$data <- NULL

## insitutype clustering by Slide
out <- call_insitu_cluster(Sobj = Sobj, 
                                n_of_cluster = 7, # there're 7 genesets for 7 minor cell types
                                assay = "RNA", 
                                return_flightpath = TRUE)
Sobj <- out$Sobj

## Check n of real clusters
Idents(Sobj) <- "Inst_Cluster_n8"
DotPlot(Sobj, features =  unique(c("Ptprc", "Cd3e", "Cd4", #TCD4
                                  "Ptprc", "Cd3e", "Cd8a", #TCD8
                                  "Ptprc", "Cd19", "Cd69", #B cell
                                  "Tnfsf13b", "Tnfsf13", #B cell
                                  "Tnfrsf13c", "Tnfrsf13b", "Tnfrsf17",#B cell
                                  "Cd27",#B cell
                                  "Sdc1"#B cell
))[]) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Analyst: it looks like we should 
# Collapse clusters d and c into Tcd4
# Collapse e and b into Bcells 
# h is Tcd8
# now f, a, and g we need to look it up on their own markers
# I added this labels in $lym_nuo_granu_recluster. For some reason I'm missing this part of the code ...


## Look up clusters f, a, and g 
library(presto)
markers <- wilcoxauc(Sobj, 'Inst_Cluster_n8')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 10) %>% 
  ungroup() %>% 
  data.frame() -> top_markers

Idents(Sobj) <- "Inst_Cluster_n8"
DotPlot(Sobj, features =  unique(top_markers$feature
))[]) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



## Get Seurat UMAP
DefaultAssay(Sobj) <- "RNA"
Sobj <- NormalizeData(Sobj)
Sobj <- FindVariableFeatures(Sobj, nfeatures = length(rownames(Sobj@assays[["RNA"]]$counts)) )
Sobj <- ScaleData(Sobj)
Sobj <- RunPCA(Sobj)
Sobj <- FindNeighbors(Sobj, dims = 1:15)
#Sobj <- FindClusters(Sobj, resolution = c(0.2, 0.5)) 
Sobj <- RunUMAP(Sobj, dims = 1:15, return.model = T)

# qsave(Sobj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_lym_Neocy_Granu_to_minor.qs")

Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_lym_Neocy_Granu_to_minor.qs")

Idents(Sobj)  <- "Inst_Cluster_n8"
DimPlot(Sobj)

head(Sobj)
table(Sobj$lym_nuo_granu_recluster, Sobj$Inst_Cluster_n8)

# 2. Get DEG, plot Vulcan plot 

# 2.1 Tcd4 DEGs 
plt_1_2 <- NULL  
plt_1_3 <- NULL  
plt_2_3 <- NULL
markers <-  c("Ptprc", "Cd3e", "Cd4") # CD45, CD3, CD4
plt_1_2 <- Plot_Volcano_DEG_sc(SeuratObj = subset(Sobj, lym_nuo_granu_recluster %in% c("Tcd4")),
                               Identity = "Group", 
                               Group_1 = comparisons$Comparison$group_1, 
                               Group_2 = comparisons$Comparison$group_2, 
                               logFC_direction = comparisons$Comparison$group_2,    
                               markers_to_label = markers,   
                               filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                               filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                               filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                               remove_mt = TRUE, 
                               genes_to_remove_from_vis = c("Hba-a1/2","Hbb"))    


# plt_1_2[[1]] # vulcan plot 
# # if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
# 
# plt_1_2[[2]] # DEG 

ggsave("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Plots/Tcd4_1_2.pdf", plot = plt_1_2[[1]] , width = 12, height = 6)

markers <-  c("Ptprc", "Cd3e", "Cd4") # CD45, CD3, CD4
plt_1_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(Sobj, lym_nuo_granu_recluster %in% c("Tcd4")),
                               Identity = "Group", 
                               Group_1 = comparisons$Comparison$group_1, 
                               Group_2 = comparisons$Comparison$group_3, 
                               logFC_direction = comparisons$Comparison$group_1,    
                               markers_to_label = markers,   
                               filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                               filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                               filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                               remove_mt = TRUE, 
                               genes_to_remove_from_vis = c("Hba-a1/2","Hbb"))

# plt_1_3[[1]] # vulcan plot 
# # if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
# 
# plt_1_3[[2]] # DEG 

ggsave("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Plots/Tcd4_1_3.pdf", plot = plt_1_3[[1]] , width = 12, height = 6)




markers <-  c("Ptprc", "Cd3e", "Cd4") # CD45, CD3, CD4
plt_2_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(Sobj, lym_nuo_granu_recluster %in% c("Tcd4")),
                               Identity = "Group", 
                               Group_1 = comparisons$Comparison$group_2, 
                               Group_2 = comparisons$Comparison$group_3, 
                               logFC_direction = comparisons$Comparison$group_2,    
                               markers_to_label = markers,   
                               filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                               filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                               filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                               remove_mt = TRUE, 
                               genes_to_remove_from_vis = c("Hba-a1/2","Hbb"))

# plt_2_3[[1]] # vulcan plot 
# # if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
# 
# plt_2_3[[2]]

ggsave("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Plots/Tcd4_2_3.pdf", plot = plt_2_3[[1]] , width = 12, height = 6)



# no padj filtering     
DEG_1_2 <- plt_1_2[[2]]    
DEG_1_3 <- plt_1_3[[2]]   
DEG_2_3 <- plt_2_3[[2]]

comparisons$Comparison

DEG_list <- list(
  Neo_RSV__Neo_RSV_IFN = DEG_1_2[DEG_1_2$logFC <= -0.5 | DEG_1_2$logFC >= 0.5, ]$feature,
  Neo_RSV__Adu_RSV = DEG_1_3[DEG_1_3$logFC <= -0.5 | DEG_1_3$logFC >= 0.5, ]$feature,
  Neo_RSV_IFN__Adu_RSV = DEG_2_3[DEG_2_3$logFC <= -0.5 | DEG_2_3$logFC >= 0.5, ]$feature
)

library(ggplot2)
library(UpSetR)
upset(fromList(DEG_list), order.by = "freq", nsets = length(DEG_list)) 

x <- upset(fromList(DEG_list), nsets = length(DEG_list))
x$New_data[1:3, 1:3]
dim(x$New_data) #

x1 <- unique(unlist(DEG_list, use.names = FALSE))
gene_1 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

gene_2 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 0 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

gene_3 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 0 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 1, ] %>% rownames() %>% as.numeric()]

RSV_IFN_gene_effects <- c(gene_1, gene_2, gene_3)

RSV_IFN_genes_df <- data.frame(genes = RSV_IFN_gene_effects,
                               CellType = "Tcd4",
                               description = "genes most likely related to the RSV_IFN effect")

comparisons$Comparison




# Filtering top genes to share with Diego
DEG_1_2_top <- rbind(DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_2_top$CellType <- "Tcd4"


DEG_1_3_top <- rbind(DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_3_top$CellType <- "Tcd4"


DEG_2_3_top <- rbind(DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_2_3_top$CellType <- "Tcd4"

library(openxlsx)
sheets_list <- list("RSV_IFN_genes_hit" = RSV_IFN_genes_df, 
                    "NeoRSV_NeoRSVIFN" = DEG_1_2_top, 
                    "NeoRSV_AduRSV" = DEG_1_3_top,
                    "NeoRSVIFN_AduRSV" = DEG_2_3_top)

write.xlsx(sheets_list, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/Tcd4_DEG_AduRSV_NeoRSV_NeoRSVIFN.xlsx")






# 2.1 Tcd8 DEGs 
plt_1_2 <- NULL  
plt_1_3 <- NULL  
plt_2_3 <- NULL

markers <- c("Ptprc", "Cd3e", "Cd8a") # CD45, CD3, CD8a
plt_1_2 <- Plot_Volcano_DEG_sc(SeuratObj = subset(Sobj, lym_nuo_granu_recluster %in% c("Tcd8")),
                               Identity = "Group", 
                               Group_1 = comparisons$Comparison$group_1, 
                               Group_2 = comparisons$Comparison$group_2, 
                               logFC_direction = comparisons$Comparison$group_2,    
                               markers_to_label = markers,   
                               filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                               filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                               filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                               remove_mt = TRUE, 
                               genes_to_remove_from_vis = c("Hba-a1/2","Hbb"))

# plt_1_2[[1]] # vulcan plot 
# # if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
# 
# plt_1_2[[2]] # DEG 

ggsave("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Plots/Tcd8_1_2.pdf", plot = plt_1_2[[1]] , width = 12, height = 6) 



markers <- c("Ptprc", "Cd3e", "Cd8a") # CD45, CD3, CD8a
plt_1_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(Sobj, lym_nuo_granu_recluster %in% c("Tcd8")),
                               Identity = "Group", 
                               Group_1 = comparisons$Comparison$group_1, 
                               Group_2 = comparisons$Comparison$group_3, 
                               logFC_direction = comparisons$Comparison$group_1,    
                               markers_to_label = markers,   
                               filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                               filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                               filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                               remove_mt = TRUE,
                               genes_to_remove_from_vis = c("Hba-a1/2","Hbb"))

# plt_1_3[[1]] # vulcan plot 
# # if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
# 
# plt_1_3[[2]] # DEG 

ggsave("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Plots/Tcd8_1_3.pdf", plot = plt_1_3[[1]] , width = 12, height = 6)




markers <- c("Ptprc", "Cd3e", "Cd8a") # CD45, CD3, CD8a
plt_2_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(Sobj, lym_nuo_granu_recluster %in% c("Tcd8")),
                               Identity = "Group", 
                               Group_1 = comparisons$Comparison$group_2, 
                               Group_2 = comparisons$Comparison$group_3, 
                               logFC_direction = comparisons$Comparison$group_2,    
                               markers_to_label = markers,   
                               filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                               filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                               filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                               remove_mt = TRUE,
                               genes_to_remove_from_vis = c("Hba-a1/2","Hbb"))

# plt_2_3[[1]] # vulcan plot 
# # if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
# 
# plt_2_3[[2]]

ggsave("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Plots/Tcd8_2_3.pdf", plot = plt_2_3[[1]] , width = 12, height = 6)



# no padj filtering       
DEG_1_2 <- plt_1_2[[2]]    
DEG_1_3 <- plt_1_3[[2]]   
DEG_2_3 <- plt_2_3[[2]]

comparisons$Comparison

DEG_list <- list(
  Neo_RSV__Neo_RSV_IFN = DEG_1_2[DEG_1_2$logFC <= -0.5 | DEG_1_2$logFC >= 0.5, ]$feature,
  Neo_RSV__Adu_RSV = DEG_1_3[DEG_1_3$logFC <= -0.5 | DEG_1_3$logFC >= 0.5, ]$feature,
  Neo_RSV_IFN__Adu_RSV = DEG_2_3[DEG_2_3$logFC <= -0.5 | DEG_2_3$logFC >= 0.5, ]$feature
)

library(ggplot2)
library(UpSetR)
upset(fromList(DEG_list), order.by = "freq", nsets = length(DEG_list)) 

x <- upset(fromList(DEG_list), nsets = length(DEG_list))
x$New_data[1:3, 1:3]
dim(x$New_data) #


x1 <- unique(unlist(DEG_list, use.names = FALSE))
gene_1 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

gene_2 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 0 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

gene_3 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 0 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 1, ] %>% rownames() %>% as.numeric()]

RSV_IFN_gene_effects <- c(gene_1, gene_2, gene_3)

RSV_IFN_genes_df <- data.frame(genes = RSV_IFN_gene_effects,
                               CellType = "Tcd8",
                               description = "genes most likely related to the RSV_IFN effect")

comparisons$Comparison




# Filtering top genes to share with Diego
DEG_1_2_top <- rbind(DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_2_top$CellType <- "Tcd8"


DEG_1_3_top <- rbind(DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_3_top$CellType <- "Tcd8"


DEG_2_3_top <- rbind(DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_2_3_top$CellType <- "Tcd8"

library(openxlsx)
sheets_list <- list("RSV_IFN_genes_hit" = RSV_IFN_genes_df, 
                    "NeoRSV_NeoRSVIFN" = DEG_1_2_top, 
                    "NeoRSV_AduRSV" = DEG_1_3_top,
                    "NeoRSVIFN_AduRSV" = DEG_2_3_top)

write.xlsx(sheets_list, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/Tcd8_DEG_AduRSV_NeoRSV_NeoRSVIFN.xlsx")





# 2.1 Bcell DEGs 
plt_1_2 <- NULL  
plt_1_3 <- NULL  
plt_2_3 <- NULL
markers <-  c("Ptprc", "Cd19", "Cd69",    "Tnfsf13b", "Tnfsf13",    "Tnfrsf13c", "Tnfrsf13b", "Tnfrsf17",   "Cd27",   "Sdc1" )
plt_1_2 <- Plot_Volcano_DEG_sc(SeuratObj = subset(Sobj, lym_nuo_granu_recluster %in% c("Bcell")),
                               Identity = "Group", 
                               Group_1 = comparisons$Comparison$group_1, 
                               Group_2 = comparisons$Comparison$group_2, 
                               logFC_direction = comparisons$Comparison$group_2,    
                               markers_to_label = markers,   
                               filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                               filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                               filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                               remove_mt = TRUE, 
                               genes_to_remove_from_vis = c("Hba-a1/2","Hbb"))    


# plt_1_2[[1]] # vulcan plot 
# # if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
# 
# plt_1_2[[2]] # DEG 

ggsave("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Plots/Bcell_1_2.pdf", plot = plt_1_2[[1]] , width = 12, height = 6)

markers <-  c("Ptprc", "Cd19", "Cd69",    "Tnfsf13b", "Tnfsf13",    "Tnfrsf13c", "Tnfrsf13b", "Tnfrsf17",   "Cd27",   "Sdc1" )
plt_1_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(Sobj, lym_nuo_granu_recluster %in% c("Bcell")),
                               Identity = "Group", 
                               Group_1 = comparisons$Comparison$group_1, 
                               Group_2 = comparisons$Comparison$group_3, 
                               logFC_direction = comparisons$Comparison$group_1,    
                               markers_to_label = markers,   
                               filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                               filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                               filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                               remove_mt = TRUE,
                               genes_to_remove_from_vis = c("Hba-a1/2","Hbb"))

# plt_1_3[[1]] # vulcan plot 
# # if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
# 
# plt_1_3[[2]] # DEG 

ggsave("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Plots/Bcell_1_3.pdf", plot = plt_1_3[[1]] , width = 12, height = 6)




markers <-  c("Ptprc", "Cd19", "Cd69",    "Tnfsf13b", "Tnfsf13",    "Tnfrsf13c", "Tnfrsf13b", "Tnfrsf17",   "Cd27",   "Sdc1" )
plt_2_3 <- Plot_Volcano_DEG_sc(SeuratObj = subset(Sobj, lym_nuo_granu_recluster %in% c("Bcell")),
                               Identity = "Group", 
                               Group_1 = comparisons$Comparison$group_2, 
                               Group_2 = comparisons$Comparison$group_3, 
                               logFC_direction = comparisons$Comparison$group_2,    
                               markers_to_label = markers,   
                               filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                               filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                               filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                               remove_mt = TRUE, 
                               genes_to_remove_from_vis = c("Hba-a1/2","Hbb"))

# plt_2_3[[1]] # vulcan plot 
# # if "Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)" just increase plot window 
# 
# plt_2_3[[2]]

ggsave("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Plots/Bcell_2_3.pdf", plot = plt_2_3[[1]] , width = 12, height = 6)



# no padj filtering     
DEG_1_2 <- plt_1_2[[2]]    
DEG_1_3 <- plt_1_3[[2]]   
DEG_2_3 <- plt_2_3[[2]]

comparisons$Comparison

DEG_list <- list(
  Neo_RSV__Neo_RSV_IFN = DEG_1_2[DEG_1_2$logFC <= -0.5 | DEG_1_2$logFC >= 0.5, ]$feature,
  Neo_RSV__Adu_RSV = DEG_1_3[DEG_1_3$logFC <= -0.5 | DEG_1_3$logFC >= 0.5, ]$feature,
  Neo_RSV_IFN__Adu_RSV = DEG_2_3[DEG_2_3$logFC <= -0.5 | DEG_2_3$logFC >= 0.5, ]$feature
)

library(ggplot2)
library(UpSetR)
upset(fromList(DEG_list), order.by = "freq", nsets = length(DEG_list)) 

x <- upset(fromList(DEG_list), nsets = length(DEG_list))
x$New_data[1:3, 1:3]
dim(x$New_data) #

x1 <- unique(unlist(DEG_list, use.names = FALSE))
gene_1 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

gene_2 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 1 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 0 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 0, ] %>% rownames() %>% as.numeric()]

gene_3 <- x1[x$New_data[x$New_data$Neo_RSV__Neo_RSV_IFN %in% 0 &
                          x$New_data$Neo_RSV_IFN__Adu_RSV %in% 1 &
                          x$New_data$Neo_RSV__Adu_RSV %in% 1, ] %>% rownames() %>% as.numeric()]

RSV_IFN_gene_effects <- c(gene_1, gene_2, gene_3)

RSV_IFN_genes_df <- data.frame(genes = RSV_IFN_gene_effects,
                               CellType = "Bcell",
                               description = "genes most likely related to the RSV_IFN effect")

comparisons$Comparison




# Filtering top genes to share with Diego
DEG_1_2_top <- rbind(DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_1_2[order(DEG_1_2$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_2_top$CellType <- "Bcell"


DEG_1_3_top <- rbind(DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_1_3[order(DEG_1_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_1_3_top$CellType <- "Bcell"


DEG_2_3_top <- rbind(DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% head(n = 20),
                     DEG_2_3[order(DEG_2_3$logFC, decreasing = TRUE), ] %>% tail(n = 20))
DEG_2_3_top$CellType <- "Bcell"

library(openxlsx)
sheets_list <- list("RSV_IFN_genes_hit" = RSV_IFN_genes_df, 
                    "NeoRSV_NeoRSVIFN" = DEG_1_2_top, 
                    "NeoRSV_AduRSV" = DEG_1_3_top,
                    "NeoRSVIFN_AduRSV" = DEG_2_3_top)

write.xlsx(sheets_list, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/Bcell_DEG_AduRSV_NeoRSV_NeoRSVIFN.xlsx")





### Extracting all intersections from upsetplot into the spreadsheet  ------
















 
## Tmp. notes 
c("Ptprc", "Cd3e", "Cd4", #TCD4
"Ptprc", "Cd3e", "Cd8a", #TCD8
"Ptprc", "Cd19", "Cd69", #B cell
"Tnfsf13b", "Tnfsf13", #B cell
"Tnfrsf13c", "Tnfrsf13b", "Tnfrsf17",#B cell
"Cd27",#B cell
"Sdc1"#B cell
)



