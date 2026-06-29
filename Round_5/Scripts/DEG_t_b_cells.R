## ====== Date: Jul 18, 2025  ===== 


### Load and define functions ---------

library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(qs)


### Load and define variables  ---------

# This subset contains only "Granulocyte", "Nuocyte", and "Lymphocyte" from our cell type call, which was good enough and should contain all the T and B cells.
Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_lym_Neocy_Granu_to_minor.qs") # this obj came from /mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Scripts/recluster_MinorCelltypes.R


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



### Differential Expression Analysis (DE) -------- 

# 1. Tcd4 DEGs across 3 groups (∴ 3 DEs)  

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


# plt_1_2[[1]] #volcano plot
# plt_1_2[[2]] #DEG list



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
# plt_1_3[[2]] # DEG 




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
# plt_2_3[[2]] # DEG list





## Getting DEG lists and intesected DEGs 

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

library(UpSetR)
upset(fromList(DEG_list), order.by = "freq", nsets = length(DEG_list))


library(ggplot2)
library(UpSetR)
library(dplyr)

# Convert DEG_list into the format used by UpSetR
x <- upset(fromList(DEG_list), nsets = length(DEG_list))

# Extract all unique genes
x1 <- unique(unlist(DEG_list, use.names = FALSE))

# Extract the new data from UpSetR (intersection matrix)
intersections_df <- x$New_data

# Get all rows (gene intersections) where any comparisons exist
combo_rows <- which(rowSums(intersections_df) > 0)

# Create list to store gene vectors for each unique binary pattern
combination_genes <- list()

# Loop through all unique binary combinations
intersection_matrix <- intersections_df[combo_rows, , drop = FALSE]
unique_combinations <- unique(intersection_matrix)

# Convert unique combinations to character vector for naming
get_combo_name <- function(row) {
  cols_on <- names(row)[which(row == 1)]
  if (length(cols_on) == 0) return("None")
  paste(cols_on, collapse = "__")
}

# Loop through all unique combinations to collect gene vectors
for (i in seq_len(nrow(unique_combinations))) {
  combo <- unique_combinations[i, ]
  # Find matching rows in x$New_data
  match_rows <- apply(intersection_matrix, 1, function(row) all(row == combo))
  gene_indices <- as.numeric(rownames(intersection_matrix)[match_rows])
  genes <- x1[gene_indices]
  
  # Create name for the column
  combo_name <- get_combo_name(combo)
  
  # Store
  combination_genes[[combo_name]] <- genes
}

# Determine the maximum gene list length
max_length <- max(sapply(combination_genes, length))

# Pad gene vectors with NA to same length
padded_genes <- lapply(combination_genes, function(g) c(g, rep(NA, max_length - length(g))))

# Create final data frame
genes_df <- as.data.frame(padded_genes, stringsAsFactors = FALSE)

# Count non-NA genes in each column
gene_counts <- sapply(genes_df, function(col) sum(!is.na(col)))

# Sort columns by count (descending)
sorted_columns <- names(sort(gene_counts, decreasing = TRUE))

# Reorder genes_df
genes_df <- genes_df[, sorted_columns]

# View sorted result
print(genes_df)

genes_df$CellType <- "Tcd4"
genes_df$description <- "Intersected genes from differential expression across groups"


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

# Start with your static sheets
sheets_list <- list(
  "Intersected_DEG" = genes_df, 
  "NeoRSV_NeoRSVIFN" = DEG_1_2_top, 
  "NeoRSV_AduRSV" = DEG_1_3_top,
  "NeoRSVIFN_AduRSV" = DEG_2_3_top
)


write.xlsx(sheets_list, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/Tcd4_DEG_AduRSV_NeoRSV_NeoRSVIFN_all_DEG_intersections.xlsx")






# 1. Tcd8 DEGs across 3 groups (∴ 3 DEs)  

plt_1_2 <- NULL  
plt_1_3 <- NULL  
plt_2_3 <- NULL
markers <-  c("Ptprc", "Cd3e", "Cd4") # CD45, CD3, CD4
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


# plt_1_2[[1]] #volcano plot
# plt_1_2[[2]] #DEG list



markers <-  c("Ptprc", "Cd3e", "Cd4") # CD45, CD3, CD4
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
# plt_1_3[[2]] # DEG 




markers <-  c("Ptprc", "Cd3e", "Cd4") # CD45, CD3, CD4
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
# plt_2_3[[2]] # DEG list





## Getting DEG lists and intesected DEGs 

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

library(UpSetR)
upset(fromList(DEG_list), order.by = "freq", nsets = length(DEG_list))


library(ggplot2)
library(UpSetR)
library(dplyr)

# Convert DEG_list into the format used by UpSetR
x <- upset(fromList(DEG_list), nsets = length(DEG_list))

# Extract all unique genes
x1 <- unique(unlist(DEG_list, use.names = FALSE))

# Extract the new data from UpSetR (intersection matrix)
intersections_df <- x$New_data

# Get all rows (gene intersections) where any comparisons exist
combo_rows <- which(rowSums(intersections_df) > 0)

# Create list to store gene vectors for each unique binary pattern
combination_genes <- list()

# Loop through all unique binary combinations
intersection_matrix <- intersections_df[combo_rows, , drop = FALSE]
unique_combinations <- unique(intersection_matrix)

# Convert unique combinations to character vector for naming
get_combo_name <- function(row) {
  cols_on <- names(row)[which(row == 1)]
  if (length(cols_on) == 0) return("None")
  paste(cols_on, collapse = "__")
}

# Loop through all unique combinations to collect gene vectors
for (i in seq_len(nrow(unique_combinations))) {
  combo <- unique_combinations[i, ]
  # Find matching rows in x$New_data
  match_rows <- apply(intersection_matrix, 1, function(row) all(row == combo))
  gene_indices <- as.numeric(rownames(intersection_matrix)[match_rows])
  genes <- x1[gene_indices]
  
  # Create name for the column
  combo_name <- get_combo_name(combo)
  
  # Store
  combination_genes[[combo_name]] <- genes
}

# Determine the maximum gene list length
max_length <- max(sapply(combination_genes, length))

# Pad gene vectors with NA to same length
padded_genes <- lapply(combination_genes, function(g) c(g, rep(NA, max_length - length(g))))

# Create final data frame
genes_df <- as.data.frame(padded_genes, stringsAsFactors = FALSE)

# Count non-NA genes in each column
gene_counts <- sapply(genes_df, function(col) sum(!is.na(col)))

# Sort columns by count (descending)
sorted_columns <- names(sort(gene_counts, decreasing = TRUE))

# Reorder genes_df
genes_df <- genes_df[, sorted_columns]

# View sorted result
print(genes_df)

genes_df$CellType <- "Tcd8"
genes_df$description <- "Intersected genes from differential expression across groups"


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

# Start with your static sheets
sheets_list <- list(
  "Intersected_DEG" = genes_df, 
  "NeoRSV_NeoRSVIFN" = DEG_1_2_top, 
  "NeoRSV_AduRSV" = DEG_1_3_top,
  "NeoRSVIFN_AduRSV" = DEG_2_3_top
)


write.xlsx(sheets_list, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/Tcd8_DEG_AduRSV_NeoRSV_NeoRSVIFN_all_DEG_intersections.xlsx")






# 2. Tcd8 DEGs across 3 groups (∴ 3 DEs)  

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


# plt_1_2[[1]] #volcano plot
# plt_1_2[[2]] #DEG list



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
# plt_1_3[[2]] # DEG 




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
# plt_2_3[[2]] # DEG list





## Getting DEG lists and intesected DEGs 

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

library(UpSetR)
upset(fromList(DEG_list), order.by = "freq", nsets = length(DEG_list))


library(ggplot2)
library(UpSetR)
library(dplyr)

# Convert DEG_list into the format used by UpSetR
x <- upset(fromList(DEG_list), nsets = length(DEG_list))

# Extract all unique genes
x1 <- unique(unlist(DEG_list, use.names = FALSE))

# Extract the new data from UpSetR (intersection matrix)
intersections_df <- x$New_data

# Get all rows (gene intersections) where any comparisons exist
combo_rows <- which(rowSums(intersections_df) > 0)

# Create list to store gene vectors for each unique binary pattern
combination_genes <- list()

# Loop through all unique binary combinations
intersection_matrix <- intersections_df[combo_rows, , drop = FALSE]
unique_combinations <- unique(intersection_matrix)

# Convert unique combinations to character vector for naming
get_combo_name <- function(row) {
  cols_on <- names(row)[which(row == 1)]
  if (length(cols_on) == 0) return("None")
  paste(cols_on, collapse = "__")
}

# Loop through all unique combinations to collect gene vectors
for (i in seq_len(nrow(unique_combinations))) {
  combo <- unique_combinations[i, ]
  # Find matching rows in x$New_data
  match_rows <- apply(intersection_matrix, 1, function(row) all(row == combo))
  gene_indices <- as.numeric(rownames(intersection_matrix)[match_rows])
  genes <- x1[gene_indices]
  
  # Create name for the column
  combo_name <- get_combo_name(combo)
  
  # Store
  combination_genes[[combo_name]] <- genes
}

# Determine the maximum gene list length
max_length <- max(sapply(combination_genes, length))

# Pad gene vectors with NA to same length
padded_genes <- lapply(combination_genes, function(g) c(g, rep(NA, max_length - length(g))))

# Create final data frame
genes_df <- as.data.frame(padded_genes, stringsAsFactors = FALSE)

# Count non-NA genes in each column
gene_counts <- sapply(genes_df, function(col) sum(!is.na(col)))

# Sort columns by count (descending)
sorted_columns <- names(sort(gene_counts, decreasing = TRUE))

# Reorder genes_df
genes_df <- genes_df[, sorted_columns]

# View sorted result
print(genes_df)

genes_df$CellType <- "Tcd8"
genes_df$description <- "Intersected genes from differential expression across groups"


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

# Start with your static sheets
sheets_list <- list(
  "Intersected_DEG" = genes_df, 
  "NeoRSV_NeoRSVIFN" = DEG_1_2_top, 
  "NeoRSV_AduRSV" = DEG_1_3_top,
  "NeoRSVIFN_AduRSV" = DEG_2_3_top
)


write.xlsx(sheets_list, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/Tcd8_DEG_AduRSV_NeoRSV_NeoRSVIFN_all_DEG_intersections.xlsx")





# 3. Bcell DEGs across 3 groups (∴ 3 DEs)  

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


# plt_1_2[[1]] #volcano plot
# plt_1_2[[2]] #DEG list



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
# plt_1_3[[2]] # DEG 




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
# plt_2_3[[2]] # DEG list





## Getting DEG lists and intesected DEGs 

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

library(UpSetR)
upset(fromList(DEG_list), order.by = "freq", nsets = length(DEG_list))


library(ggplot2)
library(UpSetR)
library(dplyr)

# Convert DEG_list into the format used by UpSetR
x <- upset(fromList(DEG_list), nsets = length(DEG_list))

# Extract all unique genes
x1 <- unique(unlist(DEG_list, use.names = FALSE))

# Extract the new data from UpSetR (intersection matrix)
intersections_df <- x$New_data

# Get all rows (gene intersections) where any comparisons exist
combo_rows <- which(rowSums(intersections_df) > 0)

# Create list to store gene vectors for each unique binary pattern
combination_genes <- list()

# Loop through all unique binary combinations
intersection_matrix <- intersections_df[combo_rows, , drop = FALSE]
unique_combinations <- unique(intersection_matrix)

# Convert unique combinations to character vector for naming
get_combo_name <- function(row) {
  cols_on <- names(row)[which(row == 1)]
  if (length(cols_on) == 0) return("None")
  paste(cols_on, collapse = "__")
}

# Loop through all unique combinations to collect gene vectors
for (i in seq_len(nrow(unique_combinations))) {
  combo <- unique_combinations[i, ]
  # Find matching rows in x$New_data
  match_rows <- apply(intersection_matrix, 1, function(row) all(row == combo))
  gene_indices <- as.numeric(rownames(intersection_matrix)[match_rows])
  genes <- x1[gene_indices]
  
  # Create name for the column
  combo_name <- get_combo_name(combo)
  
  # Store
  combination_genes[[combo_name]] <- genes
}

# Determine the maximum gene list length
max_length <- max(sapply(combination_genes, length))

# Pad gene vectors with NA to same length
padded_genes <- lapply(combination_genes, function(g) c(g, rep(NA, max_length - length(g))))

# Create final data frame
genes_df <- as.data.frame(padded_genes, stringsAsFactors = FALSE)

# Count non-NA genes in each column
gene_counts <- sapply(genes_df, function(col) sum(!is.na(col)))

# Sort columns by count (descending)
sorted_columns <- names(sort(gene_counts, decreasing = TRUE))

# Reorder genes_df
genes_df <- genes_df[, sorted_columns]

# View sorted result
print(genes_df)

genes_df$CellType <- "Bcell"
genes_df$description <- "Intersected genes from differential expression across groups"


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

# Start with your static sheets
sheets_list <- list(
  "Intersected_DEG" = genes_df, 
  "NeoRSV_NeoRSVIFN" = DEG_1_2_top, 
  "NeoRSV_AduRSV" = DEG_1_3_top,
  "NeoRSVIFN_AduRSV" = DEG_2_3_top
)


write.xlsx(sheets_list, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/Bcell_DEG_AduRSV_NeoRSV_NeoRSVIFN_all_DEG_intersections.xlsx")


