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




group_levels <- c(
  # Adult
  "Adult Control",
  "Adult RSV infected",
  "Adult reinfected",
  
  # Neonate - baseline
  "Neonate control",
  
  # Neonate - IFN upon infection
  "Neonate RSV infected (NO IFN)",
  "Neonate IFN and RSV infected",
  
  # Neonate - IFN upon reinfection
  "Neonate (NO IFN) reinfected",
  "Neonate IFN and RSV reinfected"
)


# Gene names matching the 1k mouse gene panel - from Diego
cell_types_markers_diego <- list(
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



### Load data ------

# Merged obj 
SeuObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/data/Processed/merged_Sobj_UMAP.qs")


### Lymphocytes - subtypes (Tcd4, Tcd8, Bcell) ---------
library(Seurat)
library(dplyr)

# Extract raw counts matrix from Nanostring assay
counts <- SeuObj@assays$Nanostring$counts

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

# Add annotation to Seurat meta data
SeuObj@meta.data$lyn_minor_celltype <- annot

SeuObj_sub <- subset(SeuObj, lyn_minor_celltype %in% c("Bcell", "Tcd4", "Tcd8") & major_celltype %in% "Lymphocyte")

counts_mat <- SeuObj_sub@assays$Nanostring$counts
meta_df    <- SeuObj_sub@meta.data

# create new Seurat object
SeuObj_sub <- CreateSeuratObject(
  counts = counts_mat,
  meta.data = meta_df
)
SeuObj_sub@assays$Nanostring <- SeuObj_sub@assays$RNA
SeuObj_sub@assays$RNA <- NULL
DefaultAssay(SeuObj_sub) <- "Nanostring"


table(SeuObj_sub$major_celltype)
table(SeuObj_sub$lyn_minor_celltype,
      SeuObj_sub$Group)


# Seurat workflow 
library(Seurat)
library(qs)
library(Matrix)
assay_name <- "Nanostring"
# QC
min_features <- 15
min_counts   <- 30
# Dimensionality
dimensions <- 1:15
# Clustering
cluster_res_02 <- 0.2
cluster_res_05 <- 0.5

obj <- SeuObj_sub
DefaultAssay(obj) <- assay_name

# QC
counts_mat <- GetAssayData(obj, slot = "counts")
obj$nCount   <- Matrix::colSums(counts_mat)
obj$nFeature <- Matrix::colSums(counts_mat > 0)
obj <- subset(obj, subset = nFeature >= min_features & nCount >= min_counts)

# Workflow
obj <- obj |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = nrow(obj@assays[[assay_name]]$counts)) |>
  ScaleData() |>
  RunPCA() |>
  FindNeighbors(dims = dimensions) |>
  RunUMAP(dims = dimensions) |>
  FindClusters(resolution = c(cluster_res_02, cluster_res_05))


DimPlot(obj, 
        #ols = celltype_cols,
        group.by = "lyn_minor_celltype" 
) & NoAxes()


DimPlot(obj, 
        #ols = celltype_cols,
        group.by = "Nanostring_snn_res.0.2" 
       ) & NoAxes()



library(presto)
group_col <- "lyn_minor_celltype" #lyn_minor_celltype, #Nanostring_snn_res.0.2
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, group_col)
markers %>% head()
dim(markers)


top_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  dplyr::top_n(wt = logFC, n = 30) %>%
  ungroup() 


library(ggplot2)
library(viridis)

Idents(obj) <- "lyn_minor_celltype"

DotPlot(
  obj,
  features = unique(c(top_markers$feature))
) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Markers")


library(presto)
group_col <- "Nanostring_snn_res.0.2" #lyn_minor_celltype, #Nanostring_snn_res.0.2
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, group_col)
markers %>% head()
dim(markers)


top_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  dplyr::top_n(wt = logFC, n = 10) %>%
  ungroup() 


library(ggplot2)
library(viridis)

Idents(obj) <- "Nanostring_snn_res.0.2"

DotPlot(
  obj,
  features = unique(c(top_markers$feature))
) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Markers")


library(dplyr)
library(ggplot2)

# =========================
# PARAMETERS (set these)
# =========================
group_col      <- "Group"
celltype_col   <- "lyn_minor_celltype"

# desired order (set to NULL to keep default)
group_levels    <- group_levels
celltype_levels <- c("Bcell", "Tcd4", "Tcd8")

# named colors
celltype_cols <- c("Bcell" = "#F0E442", 
                   "Tcd4" = "#E69F00", 
                   "Tcd8" = "#D55E00")

# =========================
# EXTRACT META DATA
# =========================
df <- obj@meta.data %>%
  dplyr::select(all_of(c(group_col, celltype_col))) %>%
  na.omit()

# =========================
# APPLY ORDERING
# =========================
if (!is.null(group_levels)) {
  df[[group_col]] <- factor(df[[group_col]], levels = group_levels)
}

if (!is.null(celltype_levels)) {
  df[[celltype_col]] <- factor(df[[celltype_col]], levels = celltype_levels)
}

# =========================
# COMPUTE COUNTS + PROPORTIONS
# =========================
df_summary <- df %>%
  group_by(.data[[group_col]], .data[[celltype_col]]) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(.data[[group_col]]) %>%
  mutate(prop = n / sum(n))

# =========================
# ABSOLUTE COUNTS BARPLOT
# =========================
p_counts <- ggplot(df_summary,
                   aes(x = .data[[group_col]],
                       y = n,
                       fill = .data[[celltype_col]])) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_cols, drop = FALSE) +
  theme_classic() +
  labs(x = group_col, y = "Cell count", fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_counts)

# =========================
# PROPORTION BARPLOT
# =========================
p_props <- ggplot(df_summary,
                  aes(x = .data[[group_col]],
                      y = prop,
                      fill = .data[[celltype_col]])) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_cols, drop = FALSE) +
  theme_classic() +
  labs(x = group_col, y = "Proportion", fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_props)



#qsave(obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/lymphocytes/BcellTcell_Sobj.qs")




### Dendritic cells - subtypes () ---------
library(Seurat)
library(dplyr)

# Extract raw counts matrix from Nanostring assay
counts <- SeuObj@assays$Nanostring$counts

# Initialize annotation vector with default "other"
annot <- rep("other", ncol(counts))

# Hierarchical annotation based on transcript counts

has_Itgax <- counts["Itgax", ] > 0
has_Bst2 <- counts["Bst2", ] > 0
#has_Siglech <- counts["Siglech", ] > 0 #not in the panel
has_Cd8a <- counts["Cd8a", ] > 0
has_Itgam <- counts["Itgam", ] > 0

# Assign pDC
annot[has_Itgax & has_Bst2] <- "pDC"
# Assign cDC
annot[has_Itgax & has_Cd8a] <- "cDC"



# Add annotation to Seurat meta data
SeuObj@meta.data$den_minor_celltype <- annot
table(SeuObj$major_celltype, SeuObj$den_minor_celltype)


SeuObj_sub <- subset(SeuObj, den_minor_celltype %in% c("pDC", "cDC") & major_celltype %in% c("Dendritic_cell", "Lymphocyte"))

table(SeuObj_sub$major_celltype, SeuObj_sub$den_minor_celltype)

counts_mat <- SeuObj_sub@assays$Nanostring$counts
meta_df    <- SeuObj_sub@meta.data

# create new Seurat object
SeuObj_sub <- CreateSeuratObject(
  counts = counts_mat,
  meta.data = meta_df
)
SeuObj_sub@assays$Nanostring <- SeuObj_sub@assays$RNA
SeuObj_sub@assays$RNA <- NULL
DefaultAssay(SeuObj_sub) <- "Nanostring"


table(SeuObj_sub$major_celltype)
table(SeuObj_sub$den_minor_celltype,
      SeuObj_sub$Group)


# Seurat workflow 
library(Seurat)
library(qs)
library(Matrix)
assay_name <- "Nanostring"
# QC
min_features <- 15
min_counts   <- 30
# Dimensionality
dimensions <- 1:15
# Clustering
cluster_res_02 <- 0.2
cluster_res_05 <- 0.5

obj <- SeuObj_sub
DefaultAssay(obj) <- assay_name

# QC
counts_mat <- GetAssayData(obj, slot = "counts")
obj$nCount   <- Matrix::colSums(counts_mat)
obj$nFeature <- Matrix::colSums(counts_mat > 0)
obj <- subset(obj, subset = nFeature >= min_features & nCount >= min_counts)

# Workflow
obj <- obj |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = nrow(obj@assays[[assay_name]]$counts)) |>
  ScaleData() |>
  RunPCA() |>
  FindNeighbors(dims = dimensions) |>
  RunUMAP(dims = dimensions) |>
  FindClusters(resolution = c(cluster_res_02, cluster_res_05))


DimPlot(obj, 
        #ols = celltype_cols,
        group.by = "den_minor_celltype" 
) & NoAxes()


DimPlot(obj, 
        #ols = celltype_cols,
        group.by = "Nanostring_snn_res.0.2" 
) & NoAxes()



library(presto)
group_col <- "den_minor_celltype" #den_minor_celltype, #Nanostring_snn_res.0.2
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, group_col)
markers %>% head()
dim(markers)


top_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  dplyr::top_n(wt = logFC, n = 30) %>%
  ungroup() 


library(ggplot2)
library(viridis)

Idents(obj) <- "den_minor_celltype"

DotPlot(
  obj,
  features = unique(c(top_markers$feature))
) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Markers")


library(presto)
group_col <- "Nanostring_snn_res.0.2" #lyn_minor_celltype, #Nanostring_snn_res.0.2
obj@assays$RNA <- obj@assays$Nanostring
markers <- wilcoxauc(obj, group_col)
markers %>% head()
dim(markers)


top_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  dplyr::top_n(wt = logFC, n = 10) %>%
  ungroup() 


library(ggplot2)
library(viridis)

Idents(obj) <- "Nanostring_snn_res.0.2"

DotPlot(
  obj,
  features = unique(c(top_markers$feature))
) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Markers")


library(dplyr)
library(ggplot2)

# =========================
# PARAMETERS (set these)
# =========================
group_col      <- "Group"
celltype_col   <- "den_minor_celltype"

# desired order (set to NULL to keep default)
group_levels    <- group_levels
celltype_levels <- c("pDC", "cDC")

# named colors
celltype_cols <- c("pDC" = "#FF7F00", 
                   "cDC" = "#7F3F00")

# =========================
# EXTRACT META DATA
# =========================
df <- obj@meta.data %>%
  dplyr::select(all_of(c(group_col, celltype_col))) %>%
  na.omit()

# =========================
# APPLY ORDERING
# =========================
if (!is.null(group_levels)) {
  df[[group_col]] <- factor(df[[group_col]], levels = group_levels)
}

if (!is.null(celltype_levels)) {
  df[[celltype_col]] <- factor(df[[celltype_col]], levels = celltype_levels)
}

# =========================
# COMPUTE COUNTS + PROPORTIONS
# =========================
df_summary <- df %>%
  group_by(.data[[group_col]], .data[[celltype_col]]) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(.data[[group_col]]) %>%
  mutate(prop = n / sum(n))

# =========================
# ABSOLUTE COUNTS BARPLOT
# =========================
p_counts <- ggplot(df_summary,
                   aes(x = .data[[group_col]],
                       y = n,
                       fill = .data[[celltype_col]])) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_cols, drop = FALSE) +
  theme_classic() +
  labs(x = group_col, y = "Cell count", fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_counts)

# =========================
# PROPORTION BARPLOT
# =========================
p_props <- ggplot(df_summary,
                  aes(x = .data[[group_col]],
                      y = prop,
                      fill = .data[[celltype_col]])) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_cols, drop = FALSE) +
  theme_classic() +
  labs(x = group_col, y = "Proportion", fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_props)



# qsave(obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/dendritic/pDCcDC_Sobj.qs")



### Epithelial cells - subtypes () ---------
# Need to check with them if it's worth it or not to sub-cluster epithelial. They're already so scarse across TMAs... 
# We may try a more specific marker and look into those epithekail cells. Idk ...