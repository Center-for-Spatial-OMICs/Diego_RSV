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

### Define functions -------
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
    title = "",
    genes_to_remove_from_vis = NULL# TRUE to remove mitochondrial genes
) {
  # Subset to comparison 
  Idents(SeuratObj) <- Identity
  SeuratObj_sub <- subset(x = SeuratObj, idents = c(Group_1, Group_2), invert = FALSE)
  
  # Differential Expression Analysis
  library(presto)
  markers <- presto::wilcoxauc(SeuratObj_sub, group_by = Identity)
  
  if (nrow(markers) == 0) {
    stop("No differential expression markers found with applied filters.")
  }
  
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
  markers <- markers[markers$group %in% logFC_direction, ]
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
    title = title,
    subtitle = "" )
  
  return(list(p, DEGs))
}


library(shiny)
library(ggplot2)

scatter_app <- function(sp_pre, 
                        highlight_col = NULL, 
                        col_palette = NULL, 
                        ptsize = 0.5) {
  # Set the default boundary for the image object
  DefaultBoundary(sp_pre@images$FOV) <- "centroids"
  obj <- sp_pre
  
  # Extract tissue coordinates and clean up
  data <- GetTissueCoordinates(obj)
  names(data) <- c("x", "y", "id")
  data <- data[!duplicated(data$id), ]
  
  # Ensure the cell column exists in metadata
  obj@meta.data$cell <- rownames(obj@meta.data)
  
  # Merge coordinates with metadata using different column names
  data <- merge(
    data,
    obj@meta.data[, c(highlight_col, "cell")],
    by.x = "id",
    by.y = "cell"
  )
  
  # Set a default cool palette if none provided
  if (is.null(col_palette)) {
    n_colors <- length(unique(data[[highlight_col]]))
    col_palette <- RColorBrewer::brewer.pal(min(max(n_colors, 3), 9), "PuBuGn")
  }
  
  ui <- fluidPage(
    titlePanel("Scatter Plot with Selection Feature"),
    plotOutput("scatterPlot", brush = "plot_brush"),
    verbatimTextOutput("info"),
    actionButton("export", "Export Selected Data")
  )
  
  server <- function(input, output, session) {
    output$scatterPlot <- renderPlot({
      p <- ggplot(data, aes(x = x, y = y, color = factor(get(highlight_col)))) +
        geom_point(size = ptsize) +
        labs(title = "Scatter Plot with Selectable Points", color = "Cluster") +
        theme_minimal() +
        scale_color_manual(values = col_palette)
      p
    })
    
    output$info <- renderPrint({
      req(input$plot_brush)
      brushed_points <- brushedPoints(data, input$plot_brush)
      if (nrow(brushed_points) > 0) {
        brushed_points
      } else {
        "No points selected"
      }
    })
    
    selected_data <- reactiveVal()
    
    observeEvent(input$export, {
      req(input$plot_brush)
      selected_data(brushedPoints(data, input$plot_brush))
      showNotification("Selected data exported to R environment.", type = "message")
    })
    
    observe({
      if (!is.null(selected_data())) {
        assign("extracted_coord", selected_data(), envir = .GlobalEnv)
      }
    })
  }
  
  shinyApp(ui, server)
}

# Example usage:
# scatter_app(sp_pre, highlight_col = "your_column", col_palette = c("#00BFC4", "#39B185", "#404788"), ptsize = 0.5)




### Load data ------

# Merged obj 
SeuObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/lymphocytes/BcellTcell_Sobj.qs")

SeuObj$age <- ifelse(grepl("^Adult", SeuObj$Group),
                     "adult",
                     "neonate")

# Seurat obj list (spatial layers)
obj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/Sobj_list_TMA_ID_fixed.qs")
obj_list$`NA` <- NULL
obj_list$`NA.1` <- NULL


### Visualization ------
celltype_cols <- c("Bcell" = "#F0E442", 
                   "Tcd4" = "#E69F00", 
                   "Tcd8" = "#D55E00")

celltype_levels <- c("Bcell", "Tcd4", "Tcd8")


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
        group.by = "lyn_minor_celltype", 
        cols = celltype_cols) & NoAxes()



# DimPlot(SeuObj, 
#         #cols = celltype_cols,
#         group.by = "age"
#         ) & NoAxes()


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



# Cell type prop. barplots, group by group, not TMA wise ----
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



# DE analysis, volcano plot ----

obj <- subset(SeuObj, lyn_minor_celltype %in% "Bcell")
plt <- Plot_Volcano_DEG_sc(SeuratObj = obj,
                               Identity = "age", 
                               Group_1 = "adult", 
                               Group_2 = "neonate", 
                               logFC_direction = "neonate",    
                               #markers_to_label = markers,   
                               filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                               filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                               filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                               remove_mt = TRUE, 
                               genes_to_remove_from_vis = c("Hba-a1/2","Hbb"),
                           title = "Bcell")    

plt[[1]]
plt[[2]][plt[[2]]$logFC > 0 & plt[[2]]$padj < 0.05, ]$feature #neonate
plt[[2]][plt[[2]]$logFC < 0 & plt[[2]]$padj < 0.05, ]$feature #adult
DEG <- plt[[2]]
write.csv(DEG, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/lymphocytes/DEGs_Bcell.csv")
write.csv(DEG[DEG$padj < 0.05, ], "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/lymphocytes/DEGs_Bcell_filtered_padj05.csv")




df <- as.data.frame(prop.table(table(obj$Group)))
colnames(df) <- c("Group", "Proportion")
df$Percent <- df$Proportion * 100
df[order(-df$Proportion), ]

library(dplyr)
library(ggplot2)
library(forcats)

df <- plt[[2]]

padj_cutoff <- 0.05
top_n <- 20

plot_df <- df %>%
  filter(padj < padj_cutoff) %>%
  mutate(
    group_side = ifelse(logFC < 0, "Adult", "Neonates"),
    sig = -log10(padj)
  ) %>%
  group_by(group_side) %>%
  slice_max(order_by = abs(logFC), n = top_n, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(feature = fct_reorder(feature, logFC))

ggplot(plot_df, aes(x = feature, y = logFC, fill = sig)) +
  geom_col(width = 0.8) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_viridis_c(name = "-log10(padj)") +
  labs(
    x = NULL,
    y = "logFC",
    title = paste0("Differential Expression - ", obj$lyn_minor_celltype[1] %>% as.vector()),
    subtitle = "Left: Adult | Right: Neonates"
  ) +
  theme_classic()



obj <- subset(SeuObj, lyn_minor_celltype %in% "Tcd4")
plt <- Plot_Volcano_DEG_sc(SeuratObj = obj,
                           Identity = "age", 
                           Group_1 = "adult", 
                           Group_2 = "neonate", 
                           logFC_direction = "neonate",    
                           #markers_to_label = markers,   
                           filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                           filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                           filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                           remove_mt = TRUE, 
                           genes_to_remove_from_vis = c("Hba-a1/2","Hbb"),
                           title = "Tcd4")    

plt[[1]]

plt[[2]][plt[[2]]$logFC > 0 & plt[[2]]$padj < 0.05, ]$feature #neonate
plt[[2]][plt[[2]]$logFC < 0 & plt[[2]]$padj < 0.05, ]$feature #adult
DEG <- plt[[2]]
write.csv(DEG, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/lymphocytes/DEGs_Tcd4.csv")
write.csv(DEG[DEG$padj < 0.05, ], "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/lymphocytes/DEGs_Tcd4_filtered_padj05.csv")


obj <- subset(SeuObj, lyn_minor_celltype %in% "Tcd8")
plt <- Plot_Volcano_DEG_sc(SeuratObj = obj,
                           Identity = "age", 
                           Group_1 = "adult", 
                           Group_2 = "neonate", 
                           logFC_direction = "neonate",    
                           #markers_to_label = markers,   
                           filter_logFC = NULL,      # e.g., 0.5 for abs(logFC) >= 0.5
                           filter_padj = NULL,       # e.g., 0.05 for padj < 0.05
                           filter_auc = NULL,        # e.g., 0.80 for auc > 0.80
                           remove_mt = TRUE, 
                           genes_to_remove_from_vis = c("Hba-a1/2","Hbb"),
                           title = "Tcd8")    

plt[[1]]

plt[[2]][plt[[2]]$logFC > 0 & plt[[2]]$padj < 0.05, ]$feature #neonate
plt[[2]][plt[[2]]$logFC < -0.5 & plt[[2]]$padj < 0.05, ]$feature #adult

DEG <- plt[[2]]
write.csv(DEG, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/lymphocytes/DEGs_Tcd8.csv")
write.csv(DEG[DEG$padj < 0.05, ], "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/lymphocytes/DEGs_Tcd8_filtered_padj05.csv")


df <- as.data.frame(prop.table(table(obj$Group)))
colnames(df) <- c("Group", "Proportion")
df$Percent <- df$Proportion * 100
df[order(-df$Proportion), ]


library(dplyr)
library(ggplot2)
library(forcats)
library(dplyr)
library(ggplot2)
library(forcats)

df <- plt[[2]]

padj_cutoff <- 0.05
top_n <- 20

plot_df <- df %>%
  filter(padj < padj_cutoff) %>%
  mutate(
    group_side = ifelse(logFC < 0, "Adult", "Neonates"),
    sig = -log10(padj)
  ) %>%
  group_by(group_side) %>%
  slice_max(order_by = abs(logFC), n = top_n, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(feature = fct_reorder(feature, logFC))

ggplot(plot_df, aes(x = feature, y = logFC, fill = sig)) +
  geom_col(width = 0.8) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_viridis_c(name = "-log10(padj)") +
  labs(
    x = NULL,
    y = "logFC",
    title = paste0("Differential Expression - ", obj$lyn_minor_celltype[1] %>% as.vector()),
    subtitle = "Left: Adult | Right: Neonates"
  ) +
  theme_classic()



# Spatial plot, best TMAs for each group ----
table(SeuObj$TMA_4, SeuObj$Group)

celltype_cols_lyn <- c("Bcell" = "#F0E442", 
                   "Tcd4" = "#E69F00", 
                   "Tcd8" = "#D55E00")

cell_types <- c(
  "Alveolar_bipotent_progenitor", "Alveolar_macrophage",
  "AT1_cell", "AT2_cell",
  "Dendritic_cell", "Endothelial_cell",
  "Epithelial_cell", "Granulocyte",
  "Interstitial_macrophage", "Lymphocyte",
  "Monocyte_progenitor", "Nuocyte",
  "Stromal_cell"
)

cols <- c(
  '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
  '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8',
  '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
)

# assign names (only first 13 colors will be used)
celltype_cols_major <- setNames(cols[seq_along(cell_types)], cell_types)

celltype_cols_blend <- c(celltype_cols_major, celltype_cols_lyn)



library(Seurat)
library(ggplot2)

library(patchwork)

# Adult control
obj <- obj_list[["Ap17-0332 A1 1.1"]]
# Adding lyn. subtype cell anno
SeuObj@meta.data$key <- paste(SeuObj@meta.data$TMA_4,
                              SeuObj@meta.data$cell_id_prior_merge,
                              sep = "_")

obj@meta.data$key <- paste(obj@meta.data$TMA_4,
                           rownames(obj@meta.data),
                           sep = "_")
# Create named vector for mapping
map_vec <- SeuObj@meta.data$lyn_minor_celltype
names(map_vec) <- SeuObj@meta.data$key
# Map into obj
obj@meta.data$lyn_minor_celltype <- map_vec[obj@meta.data$key]
obj@meta.data[obj@meta.data$lyn_minor_celltype %in% NA, ]$lyn_minor_celltype <- "other" 

p1 <- ImageDimPlot(
  obj, 
  size = 1, 
  dark.background = F,
  group.by = "major_celltype"  #"lyn_minor_celltype"
) +
  scale_fill_manual(values = celltype_cols_major, drop = FALSE) & NoAxes() &
  ggtitle("Adult control")


scatter_app(obj, 
            highlight_col = "major_celltype", 
            col_palette = celltype_cols_major, 
            ptsize = 0.5)

obj_sub <- subset(obj, cells = extracted_coord$id)

# blend in celltype anno
keep_types <- c("Tcd4", "Bcell", "Tcd8")
obj_sub$celltype_blend <- ifelse(
  obj_sub$lyn_minor_celltype %in% keep_types,
  obj_sub$lyn_minor_celltype,
  obj_sub$major_celltype
)


DefaultBoundary(obj_sub@images$FOV) <- "segmentation"
p2 <- ImageDimPlot(
  obj_sub,  
  size = 1, 
  dark.background = F,
  group.by = "celltype_blend"  #"lyn_minor_celltype"
) +
  scale_fill_manual(values = celltype_cols_blend, drop = FALSE) & NoAxes() &
  ggtitle("Adult control")



# Neonate control
obj <- obj_list[["Ap16-0326 A1 1"]]
# Adding lyn. subtype cell anno
SeuObj@meta.data$key <- paste(SeuObj@meta.data$TMA_4,
                              SeuObj@meta.data$cell_id_prior_merge,
                              sep = "_")

obj@meta.data$key <- paste(obj@meta.data$TMA_4,
                           rownames(obj@meta.data),
                           sep = "_")
# Create named vector for mapping
map_vec <- SeuObj@meta.data$lyn_minor_celltype
names(map_vec) <- SeuObj@meta.data$key
# Map into obj
obj@meta.data$lyn_minor_celltype <- map_vec[obj@meta.data$key]
obj@meta.data[obj@meta.data$lyn_minor_celltype %in% NA, ]$lyn_minor_celltype <- "other" 

p1 <- ImageDimPlot(
  obj, 
  size = 1, 
  dark.background = F,
  group.by = "major_celltype"  #"lyn_minor_celltype"
) +
  scale_fill_manual(values = celltype_cols_major, drop = FALSE) & NoAxes() &
  ggtitle("Neonate control")


scatter_app(obj, 
            highlight_col = "lyn_minor_celltype", 
            col_palette = celltype_cols_lyn, 
            ptsize = 0.5)

obj_sub <- subset(obj, cells = extracted_coord$id)

# blend in celltype anno
keep_types <- c("Tcd4", "Bcell", "Tcd8")
obj_sub$celltype_blend <- ifelse(
  obj_sub$lyn_minor_celltype %in% keep_types,
  obj_sub$lyn_minor_celltype,
  obj_sub$major_celltype
)


DefaultBoundary(obj_sub@images$FOV) <- "segmentation"
p2 <- ImageDimPlot(
  obj_sub,  
  size = 1, 
  dark.background = F,
  group.by = "celltype_blend"  #"lyn_minor_celltype"
) +
  scale_fill_manual(values = celltype_cols_blend, drop = FALSE) & NoAxes() &
  ggtitle("Neonate control")




# Spatial plot, best TMAs for each group (DON'T USE IT)----
table(SeuObj$TMA_4, SeuObj$Group)

library(Seurat)
library(ggplot2)

library(patchwork)

# Adult control
obj <- obj_list[["Ap17-0332 A1 1.1"]]
# Adding lyn. subtype cell anno
SeuObj@meta.data$key <- paste(SeuObj@meta.data$TMA_4,
                              SeuObj@meta.data$cell_id_prior_merge,
                              sep = "_")

obj@meta.data$key <- paste(obj@meta.data$TMA_4,
                           rownames(obj@meta.data),
                           sep = "_")
# Create named vector for mapping
map_vec <- SeuObj@meta.data$lyn_minor_celltype
names(map_vec) <- SeuObj@meta.data$key
# Map into obj
obj@meta.data$lyn_minor_celltype <- map_vec[obj@meta.data$key]
obj@meta.data[obj@meta.data$lyn_minor_celltype %in% NA, ]$lyn_minor_celltype <- "other" 

p1 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "lyn_minor_celltype"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) & NoAxes() &
  ggtitle("Adult control")



# Adult infected 
obj <- obj_list[["Ap16-3781 A1 2.1"]]

# Adding lyn. subtype cell anno
SeuObj@meta.data$key <- paste(SeuObj@meta.data$TMA_4,
                              SeuObj@meta.data$cell_id_prior_merge,
                              sep = "_")

obj@meta.data$key <- paste(obj@meta.data$TMA_4,
                           rownames(obj@meta.data),
                           sep = "_")
# Create named vector for mapping
map_vec <- SeuObj@meta.data$lyn_minor_celltype
names(map_vec) <- SeuObj@meta.data$key
# Map into obj
obj@meta.data$lyn_minor_celltype <- map_vec[obj@meta.data$key]
obj@meta.data[obj@meta.data$lyn_minor_celltype %in% NA, ]$lyn_minor_celltype <- "other" 

p2 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "lyn_minor_celltype"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) & NoAxes() &
  ggtitle("Adult infected")


# Neonate control 
obj <- obj_list[["Ap16-0326 A1 1"]]

# Adding lyn. subtype cell anno
SeuObj@meta.data$key <- paste(SeuObj@meta.data$TMA_4,
                              SeuObj@meta.data$cell_id_prior_merge,
                              sep = "_")

obj@meta.data$key <- paste(obj@meta.data$TMA_4,
                           rownames(obj@meta.data),
                           sep = "_")
# Create named vector for mapping
map_vec <- SeuObj@meta.data$lyn_minor_celltype
names(map_vec) <- SeuObj@meta.data$key
# Map into obj
obj@meta.data$lyn_minor_celltype <- map_vec[obj@meta.data$key]
obj@meta.data[obj@meta.data$lyn_minor_celltype %in% NA, ]$lyn_minor_celltype <- "other" 

p3 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "lyn_minor_celltype"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) & NoAxes() &
  ggtitle("Neonate control")


# Neonate infected 
obj <- obj_list[["Ap16-3768 A1 1"]]
obj$major_celltype <- factor(
  obj$major_celltype,
  levels = celltype_levels
)

# Adding lyn. subtype cell anno
SeuObj@meta.data$key <- paste(SeuObj@meta.data$TMA_4,
                              SeuObj@meta.data$cell_id_prior_merge,
                              sep = "_")

obj@meta.data$key <- paste(obj@meta.data$TMA_4,
                           rownames(obj@meta.data),
                           sep = "_")
# Create named vector for mapping
map_vec <- SeuObj@meta.data$lyn_minor_celltype
names(map_vec) <- SeuObj@meta.data$key
# Map into obj
obj@meta.data$lyn_minor_celltype <- map_vec[obj@meta.data$key]
obj@meta.data[obj@meta.data$lyn_minor_celltype %in% NA, ]$lyn_minor_celltype <- "other" 

p4 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "lyn_minor_celltype"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) & NoAxes() &
  ggtitle("Neonate infected")


# Neonate infected + IFN
obj <- obj_list[["Ap16-3770 A1 1"]]

# Adding lyn. subtype cell anno
SeuObj@meta.data$key <- paste(SeuObj@meta.data$TMA_4,
                              SeuObj@meta.data$cell_id_prior_merge,
                              sep = "_")

obj@meta.data$key <- paste(obj@meta.data$TMA_4,
                           rownames(obj@meta.data),
                           sep = "_")
# Create named vector for mapping
map_vec <- SeuObj@meta.data$lyn_minor_celltype
names(map_vec) <- SeuObj@meta.data$key
# Map into obj
obj@meta.data$lyn_minor_celltype <- map_vec[obj@meta.data$key]
obj@meta.data[obj@meta.data$lyn_minor_celltype %in% NA, ]$lyn_minor_celltype <- "other" 

p5 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "lyn_minor_celltype"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) & NoAxes()  &
  ggtitle("Neonate infected + IFN")


# Adult reinfected
obj <- obj_list[["Ap16-3894 A1 3"]]

# Adding lyn. subtype cell anno
SeuObj@meta.data$key <- paste(SeuObj@meta.data$TMA_4,
                              SeuObj@meta.data$cell_id_prior_merge,
                              sep = "_")

obj@meta.data$key <- paste(obj@meta.data$TMA_4,
                           rownames(obj@meta.data),
                           sep = "_")
# Create named vector for mapping
map_vec <- SeuObj@meta.data$lyn_minor_celltype
names(map_vec) <- SeuObj@meta.data$key
# Map into obj
obj@meta.data$lyn_minor_celltype <- map_vec[obj@meta.data$key]
obj@meta.data[obj@meta.data$lyn_minor_celltype %in% NA, ]$lyn_minor_celltype <- "other" 

p6 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "lyn_minor_celltype"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) & NoAxes()  &
  ggtitle("Adult reinfected")


# =========================
# PATCHWORK COMBINE
# =========================
(p1 | p2 | p3 | p4 | p5 | p6) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")






# Testing  --------------------------------------------------

# Build niche assay
obj <- BuildNicheAssay(
  object = obj,
  fov = names(obj@images)[1],          # or replace with your FOV name explicitly
  group.by = "major_celltype",     # cell type column
  neighbors.k = 100,                    # number of neighboring cells
  niches.k =   4                       # number of niche clusters to define
)


pn1 <- ImageDimPlot(
  obj, 
  size = 1,
  group.by = "niches"
) +
  #scale_fill_manual(values = celltype_cols, drop = FALSE) & 
  NoAxes()  &
  ggtitle("")

print(p6 | pn1)

grepl("Cd69", rownames(obj@assays$Nanostring$counts)) %>% table()

table(obj$niches, obj$lyn_minor_celltype)
df <- data.frame(table(obj$niches, obj$major_celltype))
df[df$Var1 %in% "2", ]

dittoBarPlot(
  subset(obj, lyn_minor_celltype %in% "other", invert = T),
  var = "lyn_minor_celltype",
  group.by = "niches",
  #color.panel = minor_celltype_concate_colors,
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)


genes <- c("Ptprc", "Cd19", "Cd69")
b_idx <- obj$lyn_minor_celltype == "Bcell"

mat <- obj@assays$Nanostring$counts[genes, , drop = FALSE]
expr_all <- colSums(mat > 0) == length(genes)

obj$Bcell_Ptprc_Cd19_Cd69 <- ifelse(b_idx & expr_all, "Yes", "No")

ImageDimPlot(
  obj,
  size = 1,
  cols = c("gray21", "red"),
  group.by = "Bcell_Ptprc_Cd19_Cd69"
) +
  ggtitle(obj$Group[1])