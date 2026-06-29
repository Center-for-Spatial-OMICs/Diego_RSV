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


### Load data ------

# Merged obj 
SeuObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/data/Processed/merged_Sobj_UMAP.qs")

# Seurat obj list (spatial layers)
obj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/Sobj_list_TMA_ID_fixed.qs")
obj_list$`NA` <- NULL
obj_list$`NA.1` <- NULL


### Visualization ------
cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- SeuObj$major_celltype %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)

celltype_levels <- c(
  # Epithelial lineage
  "Alveolar_bipotent_progenitor",
  "AT2_cell",
  "AT1_cell",
  "Epithelial_cell",
  
  # Myeloid / innate immune
  "Monocyte_progenitor",
  "Interstitial_macrophage",
  "Alveolar_macrophage",
  "Dendritic_cell",
  "Granulocyte",
  "Nuocyte",
  
  # Lymphoid / adaptive immune
  "Lymphocyte",
  
  # Structural / vascular
  "Endothelial_cell",
  "Stromal_cell"
)

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



# UMAP plot, color by cell type ----
DimPlot(SeuObj, 
        group.by = "major_celltype", 
        cols = celltype_cols) & NoAxes()


# Feature plot, cell type markers (lymphocytes, dendritic, epithelial) ----
library(Seurat)
library(ggplot2)
library(viridis)

# Bcell
FeaturePlot(
  SeuObj,
  features = "Cd19",
  order = TRUE
) & NoAxes()

# Tcd4
FeaturePlot(
  SeuObj,
  features = "Cd4",
  order = TRUE
) & NoAxes()

# Tcd8
FeaturePlot(
  SeuObj,
  features = "Cd8a",
  order = TRUE
) & NoAxes()


#pDCs
FeaturePlot(
  SeuObj,
  features = "Itgax",
  order = TRUE
) & NoAxes()

#cDCs - not too clear
FeaturePlot(
  SeuObj,
  features = "Itgam",
  order = TRUE
) & NoAxes()

#epithelial
FeaturePlot(
  SeuObj,
  features = "Epcam",
  order = TRUE
) & NoAxes()


# Dotplot, top cell type markers ----
library(presto)

# =========================
# PARAMETERS (set these)
# =========================
group_col <- "major_celltype"

# =========================
# MARKER DETECTION
# =========================
SeuObj@assays$RNA <- SeuObj@assays$Nanostring
markers <- wilcoxauc(SeuObj, group_col)
markers %>% head()
dim(markers)

# =========================
# FILTER TOP MARKERS (FC + padj)
# =========================
top_markers <- markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  dplyr::top_n(wt = logFC, n = 3) %>%
  ungroup() 


# =========================
# DOTPLOT
# =========================
library(ggplot2)
library(viridis)
SeuObj$major_celltype <- factor(
  SeuObj$major_celltype,
  levels = celltype_levels
)

Idents(SeuObj) <- "major_celltype"

DotPlot(
  SeuObj,
  features = unique(c(top_markers$feature))
) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
  scale_colour_viridis(option = "magma") +
  guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Markers")





# Spatial plot, best TMAs for each group ----
table(SeuObj$TMA_4, SeuObj$Group)

library(Seurat)
library(ggplot2)

library(patchwork)

# Adult control
obj <- obj_list[["Ap17-0331 A1 2.1"]]
obj$major_celltype <- factor(
  obj$major_celltype,
  levels = celltype_levels
)
p1 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "major_celltype"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) & NoAxes() &
  ggtitle("Adult control")


# Adult infected 
obj <- obj_list[["Ap16-3780 A1 1.1"]]
obj$major_celltype <- factor(
  obj$major_celltype,
  levels = celltype_levels
)
p2 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "major_celltype"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) & NoAxes() &
  ggtitle("Adult infected")


# Neonate control 
obj <- obj_list[["Ap16-0327 A1 1"]]
obj$major_celltype <- factor(
  obj$major_celltype,
  levels = celltype_levels
)
p3 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "major_celltype"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) & NoAxes() &
  ggtitle("Neonate control")


# Neonate infected 
obj <- obj_list[["Ap16-3768 A1 1"]]
obj$major_celltype <- factor(
  obj$major_celltype,
  levels = celltype_levels
)
p4 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "major_celltype"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) & NoAxes() &
  ggtitle("Neonate infected")


# Neonate infected + IFN
obj <- obj_list[["Ap16-2657 A2 1"]]
obj$major_celltype <- factor(
  obj$major_celltype,
  levels = celltype_levels
)
p5 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "major_celltype"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) & NoAxes()  &
  ggtitle("Neonate infected + IFN")


# =========================
# PATCHWORK COMBINE
# =========================
(p1 | p2 | p3 | p4 | p5) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")





# Cell type prop. barplots, group by group, not TMA wise ----
library(dplyr)
library(ggplot2)

# =========================
# PARAMETERS (set these)
# =========================
group_col      <- "Group"
celltype_col   <- "major_celltype"

# desired order (set to NULL to keep default)
group_levels    <- group_levels
celltype_levels <- celltype_levels

# named colors
celltype_cols <- celltype_cols

# =========================
# EXTRACT META DATA
# =========================
df <- SeuObj@meta.data %>%
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
