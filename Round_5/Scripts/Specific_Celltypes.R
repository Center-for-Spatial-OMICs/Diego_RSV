## Load packages -------
library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(qs)



## Load data -------
sobj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj_list_celltype_niche.qs")

sobj_list$Adult_1@assays$Nanostring$counts %>% rownames()
sobj_list$Adult_2@assays$Nanostring$counts %>% rownames()
sobj_list$Neonate_1@assays$Nanostring$counts %>% rownames()
sobj_list$Neonate_2@assays$Nanostring$counts %>% rownames()

gene_panel <- sobj_list$Adult_1@assays$Nanostring$counts %>% rownames() 

# Load all 4 slides merged in one single obj 
All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_merged_Sobj.qs")

comparisons <- list(
  "Comparison 1" = list(
    group_1 = "Neonate IFN and RSV infected",
    group_2 = "Adult RSV infected"
  ),
  "Comparison 2" = list(
    group_1 = "Neonate IFN and RSV reinfected",
    group_2 = "Adult reinfected"
  ),
  "Comparison 3" = list(
    group_1 = "Neonate RSV infected (NO IFN)",
    group_2 = "Neonate IFN and RSV infected"
  )
)

## Adjust Cell Type Anno. -------
# From Diego's word document 

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

markers_vec <- unlist(cell_types, use.names = FALSE)

length(intersect(gene_panel, 
                 markers_vec)) # 20 

compare_vectors <- function(gene_panel, markers) {
  # Get the variable names as strings
  gene_panel_name <- deparse(substitute(gene_panel))
  markers_name <- deparse(substitute(markers))
  
  common_genes <- intersect(gene_panel, markers)
  only_gene_panel <- setdiff(gene_panel, markers)
  only_markers <- setdiff(markers, gene_panel)
  
  comparison_table <- data.frame(
    Category = c(
      paste0("In Both (", gene_panel_name, " & ", markers_name, ")"),
      paste0("Only in ", gene_panel_name),
      paste0("Only in ", markers_name),
      paste0("Total in ", gene_panel_name),
      paste0("Total in ", markers_name),
      "Identical vectors?"
    ),
    Example_Output = c(
      length(common_genes), 
      length(only_gene_panel), 
      length(only_markers), 
      length(gene_panel), 
      length(markers), 
      identical(gene_panel, markers)
    )
  )
  return(comparison_table)
}

compare_vectors(gene_panel, markers_vec)
intersected_genes <- intersect(gene_panel, markers_vec)

filter_cell_types <- function(cell_types, intersected_genes) {
  # Recursive function to filter
  filter_recursive <- function(x) {
    if (is.list(x)) {
      # Recurse into sublists
      filtered <- lapply(x, filter_recursive)
      # Remove empty elements (lists or vectors of length 0)
      filtered <- filtered[sapply(filtered, function(y) {
        if (is.list(y)) length(y) > 0 else length(y) > 0
      })]
      return(filtered)
    } else {
      # Only keep elements in intersected_genes
      filtered_vec <- x[x %in% intersected_genes]
      return(filtered_vec)
    }
  }
  # Filter and drop empty top-level slots
  result <- filter_recursive(cell_types)
  # If the result is a list, drop empty slots at the top level as well
  if (is.list(result)) {
    result <- result[sapply(result, function(y) {
      if (is.list(y)) length(y) > 0 else length(y) > 0
    })]
  }
  return(result)
}

filtered_cell_types <- filter_cell_types(cell_types, intersected_genes)



## IGNORE
# library(presto)
# obj <- sobj_list$Adult_1
# obj@assays$RNA <- obj@assays$Nanostring
# DefaultAssay(obj) <- "RNA"
# obj <- NormalizeData(obj)
# 
# markers <- wilcoxauc(obj, 'major_celltype')
# markers %>% head()
# dim(markers) 
# # Filter DEGs based only in FC
# markers %>%
#   group_by(group) %>%
#   dplyr::filter(logFC > 0 & padj <= 0.05) %>%
#   #dplyr::top_n(wt = logFC, n = 10) %>% 
#   ungroup() %>% 
#   data.frame() -> top_markers
# 
# 
# head(top_markers)
# intersect(markers$feature, markers_vec)
# 
# top_markers[top_markers$feature %in% markers_vec[], ]$group %>% table()




## Vlnplot 
Idents(obj) <- "major_celltype"
VlnPlot(obj, features = markers_vec[1:4], pt.size = 0)



## Dotplot  

# Ordering cell type 
sobj_list$Adult_1$major_celltype %>% table() %>% names()
ordered_celltype <- c(
  "Alveolar_bipotent_progenitor",
  "AT1_cell",
  "AT2_cell",
  "Epithelial_cell",
  "Endothelial_cell",
  "Stromal_cell",
  "Monocyte_progenitor",
  "Granulocyte",
  "Alveolar_macrophage",
  "Interstitial_macrophage",
  "Dendritic_cell",
  "Lymphocyte",
  "Nuocyte"
)

## Ordering cell type, adding @Data slot
for (sample_name in c("Adult_1", "Adult_2", "Neonate_1", "Neonate_2")) {
  obj <- sobj_list[[sample_name]]
  obj@assays$RNA <- obj@assays$Nanostring
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj)
  
  # Set Idents to factor with fixed levels
  obj$major_celltype <- factor(obj$major_celltype, levels = ordered_celltype)
  Idents(obj) <- "major_celltype"
  
  # Save back to list if needed
  sobj_list[[sample_name]] <- obj
}

library(ggplot2)
library(viridis)
obj <- sobj_list$Adult_1
DefaultAssay(obj) <- "RNA"

Idents(obj) <- "major_celltype"
DotPlot(obj, 
        #features = unique(top_markers$feature)
        features =  unique(c(markers_vec[]))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16)) + ggtitle("Adult_1")


## Plot all them 
library(ggplot2)
library(viridis)
library(gridExtra) # for grid.arrange

# Store plots in a list
plot_list <- list()
for (sample_name in c("Adult_1", "Neonate_1","Adult_2", "Neonate_2")) {
  obj <- sobj_list[[sample_name]]
  DefaultAssay(obj) <- "RNA"
  Idents(obj) <- "major_celltype"
  
  p <- DotPlot(obj, features = unique(markers_vec)) +
    geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
    scale_colour_viridis(option = "magma") +
    guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          axis.text.y = element_text(size = 16)) +
    ggtitle(sample_name)
  
  plot_list[[sample_name]] <- p
}

# Arrange all plots in a 2x2 grid
grid.arrange(grobs = plot_list, ncol = 2, nrow = 2)





## Markers expression similarity across cell types across datasets 
## Each marker will have a similarity score that says wheter this gene is expressed in the same way across cell types across datasets 

library(Seurat)
library(pheatmap)
# 1. Collect average expression matrices (genes x cell types) for each dataset
expr_list <- list()
for (sample_name in c("Adult_1", "Adult_2", "Neonate_1", "Neonate_2")) {
  obj <- sobj_list[[sample_name]]
  Idents(obj) <- "major_celltype"
  avg_expr <- AverageExpression(obj, features = unique(markers_vec), group.by = "major_celltype", layer = "data")$RNA
  expr_list[[sample_name]] <- avg_expr
}

genes <- rownames(expr_list[[1]])
cell_types <- colnames(expr_list[[1]])
datasets <- names(expr_list)

# 2. For each gene, build a matrix of (cell types x datasets) and compute pairwise correlation
gene_similarity <- sapply(genes, function(g) {
  # Build expression matrix: cell types (rows) x datasets (columns)
  mat <- sapply(expr_list, function(m) m[g, cell_types])
  # Compute pairwise Pearson correlations between columns (datasets)
  cor_mat <- cor(mat, method = "pearson")
  # Return mean of upper triangle (excluding diagonal)
  mean(cor_mat[upper.tri(cor_mat)])
})

# 3. Visualize as a heatmap (genes x 1, color = mean correlation)
pheatmap(matrix(gene_similarity, ncol=1, dimnames = list(genes, "Mean_Pairwise_Correlation")),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Gene Expression Profile Similarity Across Datasets",
         color = colorRampPalette(rev(c("red", "white", "blue")))(100))





# Helper function to recursively extract all marker sets
extract_sets <- function(input, target_genes, prefix = "") {
  sets <- list()
  if (is.list(input) && !is.null(names(input))) {
    # Input is a named list: recurse into each element
    for (name in names(input)) {
      value <- input[[name]]
      full_name <- if (prefix == "") name else paste(prefix, name, sep = "$")
      sets <- c(sets, extract_sets(value, target_genes, full_name))
    }
  } else if (is.character(input)) {
    # Input is a character vector: filter and keep if any matches
    matches <- input[input %in% target_genes]
    if (length(matches) > 0) {
      name <- if (prefix == "") "genes" else prefix
      sets[[name]] <- matches
    }
  }
  sets
}


# # set 1 - markers provided
# markers_vec <- unique(unlist(cell_types, use.names = FALSE))
# 
# marker_sets <- extract_sets(cell_types, markers_vec)
# library(UpSetR)
# set_names <- sort(names(marker_sets))
# upset(
#   fromList(marker_sets),
#   sets = set_names,
#   keep.order = TRUE,
#   #order.by = "freq",
#   mb.ratio = c(0.6, 0.4),
#   nsets = length(marker_sets)
# )

# set 2 - markers available
markers_vec <- unique(unlist(filtered_cell_types, use.names = FALSE))

marker_sets <- extract_sets(filtered_cell_types, markers_vec)
library(UpSetR)
set_names <- sort(names(marker_sets))
upset(
  fromList(marker_sets),
  sets = set_names,
  keep.order = TRUE,
  #order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  nsets = length(marker_sets)
)


library(dplyr)
library(tidyr)

df <- tibble(
  set_name = names(marker_sets),
  markers = marker_sets
) %>%
  mutate(
    first_prefix = sub("\\$.*", "", set_name),      # Extract text before first $
    markers = sapply(markers, paste, collapse = ", ") # Collapse marker vectors to string
  ) %>%
  group_by(first_prefix) %>%
  summarise(
    all_markers = paste(unique(unlist(strsplit(markers, ", "))), collapse = ", ")
  )

df <- data.frame(df)
print(df)

df$first_prefix <- sub("_.*", "", df$first_prefix )
df$all_markers

df_3_majorCellTypes <- df


### Merged Data 
library(ggplot2)
library(viridis)
All_Sobj@assays$RNA <- All_Sobj@assays$Nanostring
DefaultAssay(All_Sobj) <- "RNA"
All_Sobj <- NormalizeData(All_Sobj)

Idents(All_Sobj) <- "major_celltype"
DotPlot(All_Sobj, 
        #features = unique(top_markers$feature)
        features =  unique(c(markers_vec[]))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16)) + ggtitle("All Slides")


### Seurat Processing - Sketch --------
DefaultAssay(All_Sobj) <- "RNA"
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
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = c(0.2, 0.5, 1.0))
obj <- RunUMAP(obj, dims = 1:30, return.model = T)

DefaultAssay(obj) <- "sketch"
obj <- ProjectData(
  object = obj,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:30,
  refdata = list(cluster_02_full = "sketch_snn_res.0.2",
                 cluster_05_full = "sketch_snn_res.0.5",
                 cluster_1_full = "sketch_snn_res.1")
)
# switch back to analyzing all cells
DefaultAssay(obj) <- "RNA"

Idents(obj) <- "cluster_05_full"
DimPlot(obj, reduction = "full.umap")

Idents(obj) <- "Slide"
DimPlot(obj, reduction = "full.umap")

cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- sobj_list[[1]]$major_celltype %>% table() %>% names()
cols <- setNames(cols[1:length(factor_levels)], factor_levels)

Idents(obj) <- "major_celltype"
DimPlot(obj, reduction = "full.umap", cols = cols)



## Plotting Diego's markers + Comparisons groups

## Comparison #1
Idents(obj) <- "Group"
DotPlot(subset(obj, Group %in% c(comparisons$`Comparison 1`$group_1,
                                 comparisons$`Comparison 1`$group_2)), 
        #features = unique(top_markers$feature)
        features =  unique(c(trimws(unlist(strsplit(df_3_majorCellTypes$all_markers[1], ",")))))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16)) + 
  # ggtitle("All Slides ~ seurat cluster")
  ggtitle("Comparison #1 ~ Dendritic cells (Diego)")


Idents(obj) <- "Group"
DotPlot(subset(obj, Group %in% c(comparisons$`Comparison 1`$group_1,
                                 comparisons$`Comparison 1`$group_2)), 
        #features = unique(top_markers$feature)
        features =  unique(c(trimws(unlist(strsplit(df_3_majorCellTypes$all_markers[2], ",")))))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16)) + 
  # ggtitle("All Slides ~ seurat cluster")
  ggtitle("Comparison #1 ~ Epithelial cells (Diego)")



Idents(obj) <- "Group"
DotPlot(subset(obj, Group %in% c(comparisons$`Comparison 1`$group_1,
                                 comparisons$`Comparison 1`$group_2)), 
        #features = unique(top_markers$feature)
        features =  unique(c(trimws(unlist(strsplit(df_3_majorCellTypes$all_markers[3], ",")))))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16)) + 
  # ggtitle("All Slides ~ seurat cluster")
  ggtitle("Comparison #1 ~ Lymphocyte cells (Diego)")


library(viridis)
niche_cols <- viridis(5)
factor_levels <- sobj_list[[1]]$niches_5 %>% table() %>% names()
niche_cols <- setNames(niche_cols[1:length(factor_levels)], factor_levels)
distinct_colors <- c(
  "#440154FF",  # deep purple
  "darkgreen",    # vivid blue
  "#FF5940",    # strong orange
  #"gray90",    
  #"#A65628",
  "#00BCD4",     # bright cyan
  "#FDE725FF"  # bright yellow-green
)
niche_cols <- setNames(distinct_colors[1:length(factor_levels)], factor_levels)

ImageDimPlot(subset(sobj_list[[1]], Group %in% comparisons$`Comparison 1`$group_2),
             group.by = 'niches_5',
             border.color = "white",
             fov = "FOV",
             border.size = 0.1,
             cols = niche_cols,
             # molecules = rownames(obj_comp1_1[1:4]),
             # mols.size = 0.5,
             # mols.cols = c("red", "blue", "yellow", "green"),
             # alpha = 0.5,
             size = 0.8,
             axes = F) +
  ggtitle("")

dendri_genes <- unique(c(trimws(unlist(strsplit(df_3_majorCellTypes$all_markers[1], ",")))))
epi_genes <- unique(c(trimws(unlist(strsplit(df_3_majorCellTypes$all_markers[2], ",")))))
lynph_genes <- unique(c(trimws(unlist(strsplit(df_3_majorCellTypes$all_markers[3], ",")))))


## CD4
col_pal <- viridis(10, option = "D", direction = -1)
DefaultBoundary(sobj_list[[1]]@images$FOV) <- "centroids"
ImageFeaturePlot(
  subset(sobj_list[[1]], Group %in% comparisons$`Comparison 1`$group_2),
  features = lynph_genes[3],
  fov = "FOV",
  size = 0.5,
  #molecules = lynph_genes[1],
  border.size = 0.01,
  cols = c("darkblue", col_pal[1])
) + ggtitle(paste0(comparisons$`Comparison 1`$group_2))

col_pal <- viridis(10, option = "D", direction = -1)
DefaultBoundary(sobj_list[[2]]@images$FOV) <- "centroids"
ImageFeaturePlot(
  subset(sobj_list[[2]], Group %in% comparisons$`Comparison 1`$group_2),
  features = lynph_genes[3],
  fov = "FOV",
  size = 0.5,
  #molecules = lynph_genes[1],
  border.size = 0.01,
  cols = c("darkblue", col_pal[1])
) + ggtitle(paste0(comparisons$`Comparison 1`$group_2))


col_pal <- viridis(10, option = "D", direction = -1)
DefaultBoundary(sobj_list[[3]]@images$FOV) <- "centroids"
ImageFeaturePlot(
  subset(sobj_list[[3]], Group %in% comparisons$`Comparison 1`$group_1),
  features = lynph_genes[3],
  fov = "FOV",
  size = 0.5,
  #molecules = lynph_genes[1],
  border.size = 0.01,
  cols = c("darkblue", col_pal[1])
) + ggtitle(paste0(comparisons$`Comparison 1`$group_1))


## Cav1
col_pal <- viridis(10, option = "D", direction = -1)
DefaultBoundary(sobj_list[[1]]@images$FOV) <- "centroids"
ImageFeaturePlot(
  subset(sobj_list[[1]], Group %in% comparisons$`Comparison 1`$group_2),
  features = epi_genes[2],
  fov = "FOV",
  size = 0.5,
  #molecules = lynph_genes[1],
  border.size = 0.01,
  cols = c("darkblue", col_pal[1])
) + ggtitle(paste0(comparisons$`Comparison 1`$group_2))

col_pal <- viridis(10, option = "D", direction = -1)
DefaultBoundary(sobj_list[[2]]@images$FOV) <- "centroids"
ImageFeaturePlot(
  subset(sobj_list[[2]], Group %in% comparisons$`Comparison 1`$group_2),
  features = epi_genes[2],
  fov = "FOV",
  size = 0.5,
  #molecules = lynph_genes[1],
  border.size = 0.01,
  cols = c("darkblue", col_pal[1])
) + ggtitle(paste0(comparisons$`Comparison 1`$group_2))


col_pal <- viridis(10, option = "D", direction = -1)
DefaultBoundary(sobj_list[[3]]@images$FOV) <- "centroids"
ImageFeaturePlot(
  subset(sobj_list[[3]], Group %in% comparisons$`Comparison 1`$group_1),
  features = epi_genes[2],
  fov = "FOV",
  size = 0.5,
  #molecules = lynph_genes[1],
  border.size = 0.01,
  cols = c("darkblue", col_pal[1])
) + ggtitle(paste0(comparisons$`Comparison 1`$group_1))



## Bst2
col_pal <- viridis(10, option = "D", direction = -1)
DefaultBoundary(sobj_list[[1]]@images$FOV) <- "centroids"
ImageFeaturePlot(
  subset(sobj_list[[1]], Group %in% comparisons$`Comparison 1`$group_2),
  features = dendri_genes[2],
  fov = "FOV",
  size = 0.5,
  #molecules = lynph_genes[1],
  border.size = 0.01,
  cols = c("darkblue", col_pal[1])
) + ggtitle(paste0(comparisons$`Comparison 1`$group_2))

col_pal <- viridis(10, option = "D", direction = -1)
DefaultBoundary(sobj_list[[2]]@images$FOV) <- "centroids"
ImageFeaturePlot(
  subset(sobj_list[[2]], Group %in% comparisons$`Comparison 1`$group_2),
  features = dendri_genes[2],
  fov = "FOV",
  size = 0.5,
  #molecules = lynph_genes[1],
  border.size = 0.01,
  cols = c("darkblue", col_pal[1])
) + ggtitle(paste0(comparisons$`Comparison 1`$group_2))


col_pal <- viridis(10, option = "D", direction = -1)
DefaultBoundary(sobj_list[[3]]@images$FOV) <- "centroids"
ImageFeaturePlot(
  subset(sobj_list[[3]], Group %in% comparisons$`Comparison 1`$group_1),
  features = dendri_genes[2],
  fov = "FOV",
  size = 0.5,
  #molecules = lynph_genes[1],
  border.size = 0.01,
  cols = c("darkblue", col_pal[1])
) + ggtitle(paste0(comparisons$`Comparison 1`$group_1))


## Get Morans'I score
exp <- as.data.frame(sobj_list[[1]]@assays$RNA$data[, ])
coords <- GetTissueCoordinates(sobj_list[[1]])

# res <- RunMoransI(exp, pos)

library(ape)
library(dplyr)
library(fields)
exp_test <- exp[1:4, 1:100]
coords_test <- coords[coords$cell %in% colnames(exp_test), ]
identical(coords_test$cell, colnames(exp_test)) #TRUE
rownames(coords_test) <- coords_test$cell
coords_test$cell <- NULL

head(coords_test)
exp_test[1:4, 1:4]

# Create distance matrix and weights # this part takes tooooo long
library(fields)
dist_mat <- rdist(as.matrix(coords_test))
w <- 1 / dist_mat
diag(w) <- 0  # Remove self-weight

# Now run Moran's I
library(ape)
geneA_expr <- as.numeric(exp_test[1, ])  # Expression of the first gene across cells
moran_geneA <- Moran.I(geneA_expr, w)
moran_geneA$expected #monrans score




# testing with only one TMA
exp <- as.data.frame(subset(sobj_list[[1]], Group %in% comparisons$`Comparison 1`$group_2)@assays$RNA$data[, ])
coords <- GetTissueCoordinates(subset(sobj_list[[1]], Group %in% comparisons$`Comparison 1`$group_2))

identical(coords$cell, colnames(exp)) #TRUE
rownames(coords) <- coords$cell
coords$cell <- NULL

head(coords)
exp[1:4, 1:4]

# Create distance matrix and weights # this part takes tooooo long
library(fields)
dist_mat <- rdist(as.matrix(coords))
w <- 1 / dist_mat
diag(w) <- 0  # Remove self-weight

# Now run Moran's I
library(ape)
geneA_expr <- as.numeric(exp[1, ])  # Expression of the first gene across cells
moran_geneA <- Moran.I(geneA_expr, w)
moran_geneA$observed #monrans score - close to zero means no spatial correlation; p < 0.05 indicates spatial correlation significance





### Felipe's - Cell Type Proximity Score ----------------

library(dplyr)
library(FNN)


calculate_proximity_scores <- function(
    seurat_obj,
    x_col = "",
    y_col = "",
    celltype_col = "",
    k = 5
) {
  # Extract spatial data and cell type information from the Seurat object
  spatial_data <- as.data.frame(seurat_obj@meta.data[, c(x_col, y_col, celltype_col)])
  colnames(spatial_data) <- c("CenterX", "CenterY", "CellType")
  
  # Calculate centroids for each cell type
  centroids <- spatial_data %>%
    group_by(CellType) %>%
    summarise(
      CenterX = mean(CenterX),
      CenterY = mean(CenterY)
    )
  
  # Compute pairwise distances between centroids
  dist_matrix <- as.matrix(dist(centroids[, c("CenterX", "CenterY")]))
  rownames(dist_matrix) <- centroids$CellType
  colnames(dist_matrix) <- centroids$CellType
  
  # Find k-nearest neighbors for each cell
  nearest_neighbors <- get.knnx(
    data = spatial_data[, c("CenterX", "CenterY")],
    query = spatial_data[, c("CenterX", "CenterY")],
    k = k + 1  # Include the cell itself
  )
  
  # Extract neighbor indices (excluding the cell itself)
  neighbor_indices <- nearest_neighbors$nn.index[, -1]
  
  # Initialize a matrix to store neighbor counts
  cell_types <- unique(spatial_data$CellType)
  interaction_counts <- matrix(0, nrow = length(cell_types), ncol = length(cell_types))
  rownames(interaction_counts) <- cell_types
  colnames(interaction_counts) <- cell_types
  
  # Count interactions
  for (i in 1:nrow(spatial_data)) {
    cell_type <- as.character(spatial_data$CellType[i])
    if (!(cell_type %in% rownames(interaction_counts))) {
      stop(paste("Cell type not found in interaction_counts:", cell_type))
    }
    neighbor_types <- as.character(spatial_data$CellType[neighbor_indices[i, ]])
    interaction_counts[cell_type, ] <- interaction_counts[cell_type, ] + 
      table(factor(neighbor_types, levels = cell_types))
  }
  
  # Calculate the total number of cells for each cell type
  cell_type_counts <- table(spatial_data$CellType)
  cell_type_counts <- cell_type_counts[rownames(interaction_counts)]
  
  # Normalize interaction_counts by cell_type_counts (row-wise division)
  proximity_scores <- sweep(interaction_counts, 1, cell_type_counts, "/")
  
  return(proximity_scores)
}



## proximity score for one TMA
obj_adu1_comp1_g2_core1 <- subset(sobj_list[[1]], Group %in% comparisons$`Comparison 1`$group_2 & 
                                    TMA %in% "Ap16-3779 1")

out <- calculate_proximity_scores(
    seurat_obj = obj_adu1_comp1_g2_core1,
    x_col = "CenterX_global_px",
    y_col = "CenterY_global_px",
    celltype_col = "major_celltype",
    k = 5
)

pheatmap::pheatmap(out)


## proximity score for multiple TMAs 
TMA_names <- table(sobj_list[[1]]$TMA) %>% names()
proximity_list <- list()

# Loop through each TMA
for (i in seq_along(TMA_names[1:3])) {
  # Subset Seurat object for this TMA
  obj_TMA_i <- subset(sobj_list[[1]], TMA == TMA_names[i])
  
  # Calculate proximity scores
  out <- calculate_proximity_scores(
    seurat_obj = obj_TMA_i,
    x_col = "CenterX_global_px",
    y_col = "CenterY_global_px",
    celltype_col = "major_celltype",
    k = 5
  )
  
  # Store in list with TMA name as key
  proximity_list[[TMA_names[i]]] <- out
}


## Combine and compare across TMAs 
# INPUT: proximity_list - a list of proximity score matrices, one per TMA

# 1. Find all unique cell types across all TMAs
all_cell_types <- sort(unique(unlist(lapply(proximity_list, rownames))))

# 2. Function to align a matrix to the full set of cell types
align_matrix <- function(mat, all_types) {
  aligned <- matrix(0, nrow = length(all_types), ncol = length(all_types),
                    dimnames = list(all_types, all_types))
  common <- intersect(rownames(mat), all_types)
  aligned[common, common] <- mat[common, common]
  return(aligned)
}

# 3. Align all matrices
proximity_aligned <- lapply(proximity_list, align_matrix, all_types = all_cell_types)

# 4. Compute element-wise average (consensus proximity matrix)
sum_matrix <- Reduce("+", proximity_aligned)
average_matrix <- sum_matrix / length(proximity_aligned)

# 5. Compute similarity matrix (correlation between flattened matrices)
flattened <- lapply(proximity_aligned, function(mat) as.vector(mat))
flattened_mat <- do.call(cbind, flattened)
colnames(flattened_mat) <- names(proximity_aligned) # Optional: name columns
similarity_matrix <- cor(flattened_mat)

# OUTPUTS:
# average_matrix: the consensus proximity score matrix
# similarity_matrix: TMA x TMA matrix of proximity pattern similarity

# Example: print results
print(average_matrix)
print(similarity_matrix)


pheatmap::pheatmap(average_matrix)

pheatmap::pheatmap(similarity_matrix)


sobj_list[[1]]$highlight_celltype <- "Other"
sobj_list[[1]]@meta.data[sobj_list[[1]]@meta.data$major_celltype %in% c("Endothelial_cell"), ]$highlight_celltype <- "Endothelial_cell"

ImageDimPlot(subset(sobj_list[[1]], TMA %in% c(TMA_names[1:3])),
             group.by = 'highlight_celltype',
             boundaries = "centroids",
             border.color = "white",
             fov = "FOV",
             border.size = 0.1,
             cols = c("Endothelial_cell" = "darkred", "Other" = "darkgray"),
             # molecules = rownames(obj_comp1_1[1:4]),
             # mols.size = 0.5,
             # mols.cols = c("red", "blue", "yellow", "green"),
             # alpha = 0.5,
             size = 0.8,
             axes = F) +
  ggtitle("")


sobj_list[[1]]$highlight_celltype <- "Other"
sobj_list[[1]]@meta.data[sobj_list[[1]]@meta.data$major_celltype %in% c("Stromal_cell"), ]$highlight_celltype <- "Stromal_cell"

ImageDimPlot(subset(sobj_list[[1]], TMA %in% c(TMA_names[1:3])),
             group.by = 'highlight_celltype',
             boundaries = "centroids",
             border.color = "white",
             fov = "FOV",
             border.size = 0.1,
             cols = c("Stromal_cell" = "darkred", "Other" = "darkgray"),
             # molecules = rownames(obj_comp1_1[1:4]),
             # mols.size = 0.5,
             # mols.cols = c("red", "blue", "yellow", "green"),
             # alpha = 0.5,
             size = 0.8,
             axes = F) +
  ggtitle("")



sobj_list[[1]]$highlight_celltype <- "Other"
sobj_list[[1]]@meta.data[sobj_list[[1]]@meta.data$major_celltype %in% c("Lymphocyte"), ]$highlight_celltype <- "Lymphocyte"

ImageDimPlot(subset(sobj_list[[1]], TMA %in% c(TMA_names[1:3])),
             group.by = 'highlight_celltype',
             boundaries = "centroids",
             border.color = "white",
             fov = "FOV",
             border.size = 0.1,
             cols = c("Lymphocyte" = "darkred", "Other" = "darkgray"),
             # molecules = rownames(obj_comp1_1[1:4]),
             # mols.size = 0.5,
             # mols.cols = c("red", "blue", "yellow", "green"),
             # alpha = 0.5,
             size = 0.8,
             axes = F) +
  ggtitle("")






## proximity score for multiple TMAs 
TMA_names <- table(sobj_list[[1]]$TMA) %>% names()
proximity_list <- list()

# Loop through each TMA
for (i in seq_along(TMA_names[])) {
  # Subset Seurat object for this TMA
  obj_TMA_i <- subset(sobj_list[[1]], TMA == TMA_names[i])
  
  # Calculate proximity scores
  out <- calculate_proximity_scores(
    seurat_obj = obj_TMA_i,
    x_col = "CenterX_global_px",
    y_col = "CenterY_global_px",
    celltype_col = "major_celltype",
    k = 5
  )
  
  # Store in list with TMA name as key
  proximity_list[[TMA_names[i]]] <- out
}


## Combine and compare across TMAs 
# INPUT: proximity_list - a list of proximity score matrices, one per TMA

# 1. Find all unique cell types across all TMAs
all_cell_types <- sort(unique(unlist(lapply(proximity_list, rownames))))

# 2. Function to align a matrix to the full set of cell types
align_matrix <- function(mat, all_types) {
  aligned <- matrix(0, nrow = length(all_types), ncol = length(all_types),
                    dimnames = list(all_types, all_types))
  common <- intersect(rownames(mat), all_types)
  aligned[common, common] <- mat[common, common]
  return(aligned)
}

# 3. Align all matrices
proximity_aligned <- lapply(proximity_list, align_matrix, all_types = all_cell_types)

# 4. Compute element-wise average (consensus proximity matrix)
sum_matrix <- Reduce("+", proximity_aligned)
average_matrix <- sum_matrix / length(proximity_aligned)

# 5. Compute similarity matrix (correlation between flattened matrices)
flattened <- lapply(proximity_aligned, function(mat) as.vector(mat))
flattened_mat <- do.call(cbind, flattened)
colnames(flattened_mat) <- names(proximity_aligned) # Optional: name columns
similarity_matrix <- cor(flattened_mat)

# OUTPUTS:
# average_matrix: the consensus proximity score matrix
# similarity_matrix: TMA x TMA matrix of proximity pattern similarity

# Example: print results
print(average_matrix)
print(similarity_matrix)


pheatmap::pheatmap(average_matrix)

pheatmap::pheatmap(similarity_matrix)







#### Proximity scores for Comparison #1 
comparisons$`Comparison 1`$group_2

## proximity score for multiple TMAs 
TMA_names <- sobj_list[[1]]@meta.data[sobj_list[[1]]@meta.data$Group %in% "Adult RSV infected", ]$TMA %>% unique()

proximity_list <- list()

# Loop through each TMA
for (i in seq_along(TMA_names[])) {
  # Subset Seurat object for this TMA
  obj_TMA_i <- subset(sobj_list[[1]], TMA == TMA_names[i])
  
  # Calculate proximity scores
  out <- calculate_proximity_scores(
    seurat_obj = obj_TMA_i,
    x_col = "CenterX_global_px",
    y_col = "CenterY_global_px",
    celltype_col = "major_celltype",
    k = 5
  )
  
  # Store in list with TMA name as key
  proximity_list[[TMA_names[i]]] <- out
}


## Combine and compare across TMAs 
# INPUT: proximity_list - a list of proximity score matrices, one per TMA

# 1. Find all unique cell types across all TMAs
all_cell_types <- sort(unique(unlist(lapply(proximity_list, rownames))))

# 2. Function to align a matrix to the full set of cell types
align_matrix <- function(mat, all_types) {
  aligned <- matrix(0, nrow = length(all_types), ncol = length(all_types),
                    dimnames = list(all_types, all_types))
  common <- intersect(rownames(mat), all_types)
  aligned[common, common] <- mat[common, common]
  return(aligned)
}

# 3. Align all matrices
proximity_aligned <- lapply(proximity_list, align_matrix, all_types = all_cell_types)

# 4. Compute element-wise average (consensus proximity matrix)
sum_matrix <- Reduce("+", proximity_aligned)
average_matrix <- sum_matrix / length(proximity_aligned)

# 5. Compute similarity matrix (correlation between flattened matrices)
flattened <- lapply(proximity_aligned, function(mat) as.vector(mat))
flattened_mat <- do.call(cbind, flattened)
colnames(flattened_mat) <- names(proximity_aligned) # Optional: name columns
similarity_matrix <- cor(flattened_mat)

# OUTPUTS:
# average_matrix: the consensus proximity score matrix
# similarity_matrix: TMA x TMA matrix of proximity pattern similarity

# Example: print results
Comp1_group2_average_matrix <- average_matrix
print(Comp1_group2_average_matrix)
pheatmap::pheatmap(Comp1_group2_average_matrix, cluster_rows = FALSE)





comparisons$`Comparison 1`$group_1

## proximity score for multiple TMAs 
TMA_names <- sobj_list[[3]]@meta.data[sobj_list[[3]]@meta.data$Group %in% "Neonate IFN and RSV infected", ]$TMA %>% unique()

proximity_list <- list()

# Loop through each TMA
for (i in seq_along(TMA_names[])) {
  # Subset Seurat object for this TMA
  obj_TMA_i <- subset(sobj_list[[3]], TMA == TMA_names[i])
  
  # Calculate proximity scores
  out <- calculate_proximity_scores(
    seurat_obj = obj_TMA_i,
    x_col = "CenterX_global_px",
    y_col = "CenterY_global_px",
    celltype_col = "major_celltype",
    k = 5
  )
  
  # Store in list with TMA name as key
  proximity_list[[TMA_names[i]]] <- out
}


## Combine and compare across TMAs 
# INPUT: proximity_list - a list of proximity score matrices, one per TMA

# 1. Find all unique cell types across all TMAs
all_cell_types <- sort(unique(unlist(lapply(proximity_list, rownames))))

# 2. Function to align a matrix to the full set of cell types
align_matrix <- function(mat, all_types) {
  aligned <- matrix(0, nrow = length(all_types), ncol = length(all_types),
                    dimnames = list(all_types, all_types))
  common <- intersect(rownames(mat), all_types)
  aligned[common, common] <- mat[common, common]
  return(aligned)
}

# 3. Align all matrices
proximity_aligned <- lapply(proximity_list, align_matrix, all_types = all_cell_types)

# 4. Compute element-wise average (consensus proximity matrix)
sum_matrix <- Reduce("+", proximity_aligned)
average_matrix <- sum_matrix / length(proximity_aligned)

# 5. Compute similarity matrix (correlation between flattened matrices)
flattened <- lapply(proximity_aligned, function(mat) as.vector(mat))
flattened_mat <- do.call(cbind, flattened)
colnames(flattened_mat) <- names(proximity_aligned) # Optional: name columns
similarity_matrix <- cor(flattened_mat)

# OUTPUTS:
# average_matrix: the consensus proximity score matrix
# similarity_matrix: TMA x TMA matrix of proximity pattern similarity

# Example: print results
Comp1_group1_average_matrix <- average_matrix
print(Comp1_group1_average_matrix)
pheatmap::pheatmap(Comp1_group1_average_matrix, cluster_rows = FALSE)
                    












### Dend, Epi, and Lynph DEG across comparisons  ----------------


# write_csv(data.frame(df_3_majorCellTypes), "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/markers_Diego_on_panel.csv")
library(readr)
markers_Diego <- read_csv("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/markers_Diego_on_panel.csv") %>% data.frame()

dendri_genes <- unique(c(trimws(unlist(strsplit(df_3_majorCellTypes$all_markers[1], ",")))))
epi_genes <- unique(c(trimws(unlist(strsplit(df_3_majorCellTypes$all_markers[2], ",")))))
lynph_genes <- unique(c(trimws(unlist(strsplit(df_3_majorCellTypes$all_markers[3], ",")))))

All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_merged_Sobj.qs")


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

All_Sobj@assays$RNA <- All_Sobj@assays$Nanostring
DefaultAssay(All_Sobj) <- "RNA"
All_Sobj <- NormalizeData(All_Sobj)

plt_dend <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Dendritic_cell")),
                        Identity = "Group", 
                        Group_1 = "Neonate IFN and RSV infected", 
                        Group_2 = "Adult RSV infected", 
                        markers_Diego = dendri_genes)

plt_dend[[2]][plt_dend[[2]]$feature %in% dendri_genes, ] #zero 


plt_epi <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Epithelial_cell")),
                                    Identity = "Group", 
                                    Group_1 = "Neonate IFN and RSV infected", 
                                    Group_2 = "Adult RSV infected", 
                               markers_Diego = epi_genes)

plt_epi
plt_epi[[2]][plt_epi[[2]]$feature %in% epi_genes, ] #zero 



plt_lyn <- Plot_Volcano_DEG_sc(SeuratObj = subset(All_Sobj, major_celltype %in% c("Lymphocyte")),
                                   Identity = "Group", 
                                   Group_1 = "Neonate IFN and RSV infected", 
                                   Group_2 = "Adult RSV infected",
                               markers_Diego = lynph_genes)

plt_lyn[[2]][plt_lyn[[2]]$feature %in% lynph_genes, ] #zero 







comparisons$`Comparison 1`$group_1
comparisons$`Comparison 1`$group_2




filtered_cell_types %>% print()

## IGNORE 
# comparing the gene list they gave us vs the actual gene on the panel
# library(dplyr)
# library(tidyr)
# library(UpSetR)
# 
# # Extract sets
# marker_sets_1 <- extract_sets(cell_types, unique(unlist(cell_types, use.names = FALSE)))
# marker_sets_2 <- extract_sets(filtered_cell_types, unique(unlist(filtered_cell_types, use.names = FALSE)))
# 
# # All unique genes
# all_genes <- unique(c(
#   unlist(marker_sets_1, use.names = FALSE),
#   unlist(marker_sets_2, use.names = FALSE)
# ))
# 
# # Create presence/absence data frame
# make_binary_matrix <- function(marker_sets, all_genes, prefix) {
#   mat <- sapply(marker_sets, function(set) all_genes %in% set)
#   colnames(mat) <- paste0(prefix, "_", colnames(mat))
#   as.data.frame(mat)
# }
# 
# df1 <- make_binary_matrix(marker_sets_1, all_genes, "before")
# df2 <- make_binary_matrix(marker_sets_2, all_genes, "after")
# 
# # Combine, ensuring no duplicate columns
# combined_df <- cbind(Gene = all_genes, df1, df2)
# # Remove duplicate columns if any
# combined_df <- combined_df[, !duplicated(colnames(combined_df))]
# 
# # Remove the Gene column for UpSetR
# upset_df <- combined_df[ , -1]
# rownames(upset_df) <- combined_df$Gene
# 
# upset_df[] <- lapply(upset_df, function(x) as.integer(x))
# 
# # Plot
# upset(upset_df, 
#       sets = colnames(upset_df),
#       order.by = "freq", 
#       mb.ratio = c(0.6, 0.4),
#       nsets = ncol(upset_df))




