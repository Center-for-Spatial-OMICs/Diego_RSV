## ====== Date: Aug 13, 2025  ===== 


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
    pCutoff = 0.05,
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




### Create scObj for epithelial cells based on our general cell annotation ------
# Sobj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj_list_celltype_niche.qs")

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

table(All_Sobj$major_celltype)

markers <- c(
  "Epcam",
  "Scgb1a1",
  "Pdpn", "Cav1",
  "Sftpc",
  "Muc5ac",
  "Tuba1a", "Foxj1",
  "Krt5",
  "Calca", "Chga"
)

table(rownames(All_Sobj@assays$Nanostring$counts) %in% markers)
intersect(rownames(All_Sobj@assays$Nanostring$counts), markers)
# Cav1 (AT1 (podoplanin, caveolin)) 
# Epcam (General epithelial marker)

Idents(All_Sobj) <- "major_celltype"
DotPlot(All_Sobj, features = unique(c(
  # Epithelial cells markers
  "Epcam",                    # General epithelial marker
  
  # Subtypes
  "Scgb1a1",                  # Club_Cells (CC10/SCGB1A1)
  "Pdpn", "Cav1",             # AT1 (podoplanin, caveolin)
  "Sftpc",                    # AT2
  "Muc5ac",                   # Goblet_cells
  "Tuba1a", "Foxj1",          # Ciliated_cells (α-Ac-Tub = Tuba1a, Foxj1)
  "Krt5",                     # Basal_cells
  "Calca", "Chga"             # Neuroendocrine_cells
))) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Maybe it will not happen. Too few markers. And Cav1 will be probably underneth Epcam ... so there's not even 2 minor cell types based on Diego's markers 

# Let me share with him the epithelial subclustering I did so maybe they can annotate it based on the 1k panel we have for his data 
# Code for the reclustering from /mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Scripts/recluster_celltypes_DEGs.R


obj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_epithelial.qs")

library(presto)
markers <- wilcoxauc(obj, 'Resolved_clusters')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 10) %>% 
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


obj@meta.data[obj@meta.data$Resolved_clusters %in% 
                c("6_res_1",
                  "5_res_0.2",
                  "2_res_0.5", 
                  "10_res_1"), ]

# Good subclusters based on /mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Scripts/recluster_celltypes_DEGs.R
obj_sub <- subset(obj, Resolved_clusters %in% c("6_res_1",
                                     "5_res_0.2",
                                     "2_res_0.5", 
                                     "10_res_1"))

# Plotting only good epithelial subcluster markers - adding Epcam and Cav1 by the end
DefaultAssay(obj_sub) <- "RNA"
Idents(obj_sub) <- "Resolved_clusters"
DotPlot(obj_sub, 
        #features = unique(top_markers$feature)
        features =  unique(c(top_markers$feature, "Epcam","Cav1"))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14)) + ggtitle("")


p <- dittoBarPlot(
  obj_sub,
  var = "Resolved_clusters",
  group.by = "Group",
  #color.panel = ,
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)
print(p)

fills <- ggplot_build(p)$data[[1]]$fill
col_palette <- unique(fills)

col_palette <- setNames(
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
  c("10_res_1", "2_res_0.5", "5_res_0.2", "6_res_1")
)

## All groups 
p <- dittoBarPlot(
  obj_sub,
  var = "Resolved_clusters",
  group.by = "Group",
  color.panel = col_palette,
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)
print(p)

plot_and_data <- dittoBarPlot(
  obj_sub,
  var = "Resolved_clusters",
  group.by = "Group",
  data.out = TRUE
)

# Extract data frame used for plotting
plot_data <- plot_and_data$data
head(plot_data)

# Create the plot with value labels
library(ggplot2)
library(scales)  # for percent formatting

# Assuming plot_data has columns: grouping, percent, label (or the appropriate column names)
ggplot(plot_data, aes(x = grouping, y = count, fill = label)) +
  geom_bar(stat = "identity", position = "stack") +  # stacked barplot with absolute counts
  scale_fill_manual(values = col_palette) +          # apply your color palette here
  geom_text(
    aes(label = count),
    position = position_stack(vjust = 0.5),  # label centered in each stacked segment
    color = "black",
    size = 3
  ) +
  labs(x = "Group", y = "Count", fill = "Resolved_clusters") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(colour = "black", size = 12, vjust = -1), 
    axis.title.y = element_text(colour = "black", size = 12), 
    axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)
  )


## Comparison groups 
obj_sub_2 <- subset(obj_sub, Group %in% c(comparisons$Comparison$group_1,
                                          comparisons$Comparison$group_2,
                                          comparisons$Comparison$group_3))

p <- dittoBarPlot(
  obj_sub_2,
  var = "Resolved_clusters",
  group.by = "Group",
  color.panel = col_palette,
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)
print(p)

plot_and_data <- dittoBarPlot(
  obj_sub_2,
  var = "Resolved_clusters",
  group.by = "Group",
  data.out = TRUE
)

# Extract data frame used for plotting
plot_data <- plot_and_data$data
head(plot_data)

# Create the plot with value labels
library(ggplot2)
library(scales)  # for percent formatting

# Assuming plot_data has columns: grouping, percent, label (or the appropriate column names)
ggplot(plot_data, aes(x = grouping, y = count, fill = label)) +
  geom_bar(stat = "identity", position = "stack") +  # stacked barplot with absolute counts
  scale_fill_manual(values = col_palette) +          # apply your color palette here
  geom_text(
    aes(label = count),
    position = position_stack(vjust = 0.5),  # label centered in each stacked segment
    color = "black",
    size = 3
  ) +
  labs(x = "Group", y = "Count", fill = "Resolved_clusters") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(colour = "black", size = 12, vjust = -1), 
    axis.title.y = element_text(colour = "black", size = 12), 
    axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)
  )



# export gene markers table 
subcluster_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 20) %>% 
  ungroup()  %>% as.data.frame()



## Export epithelail subcluster top markers + 1k panel
# Load Excel-writing library
library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Example objects (replace with your real ones)
sheet1_data <- subcluster_markers
sheet2_data <- rownames(obj_sub_2@assays$RNA$counts)

# Add first object (as a DataFrame) to sheet 1
addWorksheet(wb, "Subcluster_Markers")
writeData(wb, "Subcluster_Markers", sheet1_data)

# Add second object (rownames vector) to sheet 2
# Convert vector to a data.frame
addWorksheet(wb, "1k_CosMx_panel_genes")
writeData(wb, "1k_CosMx_panel_genes", data.frame(Gene = sheet2_data))

# Save the workbook
saveWorkbook(
  wb,
  file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/epi_subcluster_markers_top20.xlsx",
  overwrite = TRUE
)






















