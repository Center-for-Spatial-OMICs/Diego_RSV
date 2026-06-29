## ====================================
##   Project: Diego's CosMx RSV, StJude
##.  Objective: Merge exploratory analysis into one report with the proper groups the colaborators are most interested - see {Project Scope}
##   Author: Maycon Marção
##   Date: March 19, 2025
##   Version: Round_3 analysis 
## ====================================

## ------ Project Scope ------

### Directions from Diego: 
# The main comparison we are trying to do are:
# Neonate IFN infected vs Adult infected
# Neonate IFN reinfected vs Adult reinfected
# Neonate NO IFN infected vs Neonate IFN infected
# 
# Is there any chance to have the following groups in one PDF for us to compare?
#   
# Neonate IFN infected
# Adult infected
# Neonate IFN reinfected
# Adult reinfected
# Neonate NO IFN infected
# Neonate IFN infected
# 
# Also, while at the proportions by TMA the colors and letters vary from file to file. Is there any way to keep those consistent?


## 1. We must have our report in only 1 pdf
## 2. Adjust it to show only the groups they've highlighted 
## 3. Use only one color code through all analysis 


## ------ Load Libraries ------
library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(shiny)
library(plotly)
library(ggplot2)
library(RColorBrewer)
library(Seurat)


## ------ Define Functions ------



## ------ Definge Variables ------- 

# # DEG from GeoMx - I assume they don't need this - this DEG  didn't take account within groups it just compares Adult 1 vs Neonate 1 
# library(readr)
# Adu1_Neo1_GeoMx_DEG_YL <- read_csv("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adu1_Neo1_GeoMx_DEG_YL.csv") %>% data.frame()
# 
# names(Adu1_Neo1_GeoMx_DEG_YL)[names(Adu1_Neo1_GeoMx_DEG_YL) == "...1"] <- "gene"
# 
# dim(Adu1_Neo1_GeoMx_DEG_YL)
# head(Adu1_Neo1_GeoMx_DEG_YL)
# 
# deg_filt_neonate <- Adu1_Neo1_GeoMx_DEG_YL %>% filter(adj.P.Val <= 0.05 &
#                                                         logFC > 0)
# deg_filt_adult <- Adu1_Neo1_GeoMx_DEG_YL %>% filter(adj.P.Val <= 0.05 &
#                                                       logFC < 0)


## ------ Loadings from processed data -------

# Most recent obj I have for this analysis from Mar 21, 2025
obj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Objects/all_4_slides_processed.qs")

## ------ Prepare the Data ------

## Load seurat objects
# Note: SeuratObj merged from /mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Scripts/merge_SeuobjList.R
All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Objects/all_4_slides_one_Seuobj_and_Metadata.qs")
head(All_Sobj@meta.data)

# Remove control TMA 
obj <- subset(All_Sobj, fov %in% as.character(165:180) & Slide == "Adult_1", invert = TRUE)


## Run clustering (insitutype) (run it for 10 and 15 clusters)
unsup <- insitutype(
  x = as.matrix(t(obj@assays$Nanostring$counts)),
  neg = obj@meta.data$nCount_negprobes,
  reference_profiles = NULL,
  bg = NULL,
  #n_clusts = 10,
  n_clusts = 15,
  n_phase1 = 200,
  n_phase2 = 500,
  n_phase3 = 2000,
  n_starts = 1,
  max_iters = 5)

# unsup_cl10 <- unsup
unsup_cl15 <- unsup

obj$insitutype_cluster <- 'other'
obj$insitutype_cluster[names(unsup$clust)] <- unsup$clust
# names(obj@meta.data)[names(obj@meta.data) == "insitutype_cluster"] <- paste0("insitutype_cluster_", "10", "_clusters")
names(obj@meta.data)[names(obj@meta.data) == "insitutype_cluster"] <- paste0("insitutype_cluster_", "15", "_clusters")
# obj_bk <- obj





## ------ Data Analysis ------
obj <- obj %>% NormalizeData()


## Evaluating n of clusters/cluster resolution
# Define groups to compare 
Groups <- obj@meta.data$Group %>% table() %>% names() 
# Ordering clusters
obj@meta.data$insitutype_cluster_10_clusters <- factor(obj@meta.data$insitutype_cluster_10_clusters, levels = names(table(obj@meta.data$insitutype_cluster_10_clusters)))
# Ordering groups
obj@meta.data$Group <- factor(obj@meta.data$Group, levels = Groups)

# Ordering clusters
obj@meta.data$insitutype_cluster_15_clusters <- factor(obj@meta.data$insitutype_cluster_15_clusters, levels = names(table(obj@meta.data$insitutype_cluster_15_clusters)))


# unsup <- unsup_cl10
cols = brewer.pal(10, 'Paired')
unsup <- unsup_cl15
cols <- cols[seq_along(unique(unsup$clust))]
names(cols) <- unique(unsup$clust)
fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup, col = cols[unsup$clust])
# fp_cl10 <- fp
fp_cl15 <- fp


print(fp_cl10) 
print(fp_cl15) # seems better (more 90% pure clusters)

table(obj$insitutype_cluster_10_clusters)
table(obj$insitutype_cluster_15_clusters)

# Checking the most constant cluster - clustree 
library(data.table)
library(dplyr)
library(clustree)
meta <- obj@meta.data
head(meta[, c("insitutype_cluster_10_clusters", "insitutype_cluster_15_clusters")])
# Subset
cluster_df_sub <- meta[, c("insitutype_cluster_10_clusters", "insitutype_cluster_15_clusters")]
# Rename it
names(cluster_df_sub) <- c("insitutype_row_1", "insitutype_row_2") # same order as in names(cluster_df_sub). It should follow the crescent resolution of your clusters 

# Run clustree
library(clustree)
clustree_out <- clustree(cluster_df_sub[, ], prefix = "insitutype_row_", node_colour = "sc3_stability")

print(clustree_out)
# Conclusion: most constant clusters (from insitutype_cluster_15_clusters) are g, o, b, a, k, m, n.


## Set colors 
cols <-
  c(
    '#8DD3C7',
    '#BEBADA',
    '#FB8072',
    '#80B1D3',
    '#FDB462',
    '#B3DE69',
    '#FCCDE5',
    '#D9D9D9',
    '#BC80BD',
    '#CCEBC5',
    '#FFED6F',
    '#E41A1C',
    '#377EB8',
    '#4DAF4A',
    '#984EA3',
    '#FF7F00',
    '#FFFF33',
    '#A65628',
    '#F781BF',
    '#999999'
  )

factor_levels <- obj$insitutype_cluster_15_clusters %>% table() %>% names()
cols <- setNames(cols[1:length(factor_levels)], factor_levels)


## Other processing steps ----
## Remove unlabeled cells - won't effect the analysis, there're just 2348 cells that got no FOV attached to them when I was preparing the objects. 
obj <- subset(obj, TMA_2 %in% "", invert = TRUE)

## Get UMAP calculated 
obj <- obj %>%
  #NormalizeData() %>% #already did it
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  #FindClusters() %>% #already have insitutype 
  RunUMAP(dims = 1:15)

# qsave(obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Objects/all_4_slides_processed.qs")

Idents(obj) <- "Group"
DimPlot(obj)

Idents(obj) <- "Slide"
DimPlot(obj)

Idents(obj) <- "insitutype_cluster_15_clusters"
DimPlot(obj)


## Proportion plot ----
library(dittoSeq)
dittoBarPlot(obj, "insitutype_cluster_15_clusters", 
             group.by = "TMA_2",
             split.by = "Slide",
             color.panel = cols, 
             main = "Clusters by TMAs (ROIs)")

library(dittoSeq)
dittoBarPlot(obj, "insitutype_cluster_15_clusters", 
             group.by = "Slide",
             color.panel = cols,
             main = "Clusters by Slides")

library(dittoSeq)
dittoBarPlot(obj, "insitutype_cluster_15_clusters", 
             group.by = "Group",
             color.panel = cols,
             main = "Clusters by Group")

# Neonate IFN infected vs Adult infected
library(dittoSeq)
bplt_comp1 <- dittoBarPlot(subset(obj, Group %in% c("Neonate IFN and RSV infected",
                                      "Adult RSV infected")), "insitutype_cluster_15_clusters", 
             group.by = "Group",
             color.panel = cols,
             main = "Clusters by Group")

print(bplt_comp1)

# Neonate IFN reinfected vs Adult reinfected
library(dittoSeq)
bplt_comp2 <- dittoBarPlot(subset(obj, Group %in% c("Neonate IFN and RSV reinfected",
                                      "Adult reinfected")), "insitutype_cluster_15_clusters", 
             group.by = "Group",
             color.panel = cols,
             main = "Clusters by Group")

print(bplt_comp2)

# Neonate NO IFN infected vs Neonate IFN infected
library(dittoSeq)
bplt_comp3 <- dittoBarPlot(subset(obj, Group %in% c("Neonate RSV infected (NO IFN)",
                                                    "Neonate IFN and RSV infected")), "insitutype_cluster_15_clusters", 
                           group.by = "Group",
                           color.panel = cols,
                           main = "Clusters by Group")

print(bplt_comp3)


## Spatial plot (ImageDimPlot) ----

# Neonate IFN infected vs Adult infected
obj_comp1_1 <- subset(obj, Group %in% c("Neonate IFN and RSV infected"))


## Cropped plots  ----
seurat_obj <- obj_comp1_1
dim(seurat_obj)

Idents(seurat_obj) <- "fov"
#seurat_obj <- subset(x = seurat_obj, downsample = 100)
dim(seurat_obj)

# Getting coordinates - only keep one cell coordinate
data <- GetTissueCoordinates(seurat_obj)
names(data) <- c("x", "y", "id")
data <- data[!duplicated(data$id), ]

seurat_obj@meta.data$cell <- rownames(seurat_obj@meta.data)

dim(data) #66740
data <- merge(data, seurat_obj@meta.data[, c("insitutype_cluster_15_clusters", "cell")], by.x = "id", by.y = "cell")  # Different column names
dim(data) #66740

head(data)

## Define and run shiny app
library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Scatter Plot with Selection Feature"),
  plotOutput("scatterPlot", brush = "plot_brush"),  # Enable brushing for selection
  verbatimTextOutput("info"),
  actionButton("export", "Export Selected Data")  # Button to export selected data
)



server <- function(input, output) {
  output$scatterPlot <- renderPlot({
    ggplot(data, aes(x = x, y = y, color = factor(insitutype_cluster_15_clusters))) +
      geom_point(size = 0.5) +  # Adjust point size as needed
      scale_color_brewer(palette = "Paired") +  # Use the 'Paired' color palette
      labs(title = "Scatter Plot with Selectable Points", color = "Cluster") +
      theme_minimal()
  })
  
  output$info <- renderPrint({
    req(input$plot_brush)  # Ensure there is a brush selection
    
    # Get the coordinates of the brushed area
    brushed_points <- brushedPoints(data, input$plot_brush)
    
    if (nrow(brushed_points) > 0) {
      brushed_points  # Display the selected points
    } else {
      "No points selected"
    }
  })
  
  # Reactive value to store selected data
  selected_data <- reactiveVal()
  
  observeEvent(input$export, {
    req(input$plot_brush)  # Ensure there is a brush selection
    
    # Get the brushed points and store them in a reactive variable
    selected_data(brushedPoints(data, input$plot_brush))
    
    # Optionally print a message to confirm export
    showNotification("Selected data exported to R environment.", type = "message")
  })
  
  # Make the selected data accessible in the R environment with the name 'extracted_coord'
  observe({
    if (!is.null(selected_data())) {
      assign("extracted_coord", selected_data(), envir = .GlobalEnv)
    }
  })
}

shinyApp(ui, server)

head(extracted_coord)

DefaultBoundary(obj_comp1_1[["FOV.2"]]) <- "segmentation"
splt_comp1_1 <- ImageDimPlot(subset(obj_comp1_1, cells = extracted_coord$id),
                             group.by = 'insitutype_cluster_15_clusters',
                             border.color = "white",
                             fov = "FOV.2",
                             border.size = 0.1,
                             cols = cols,
                             # molecules = rownames(obj_comp1_1[1:4]),
                             # mols.size = 0.5,
                             # mols.cols = c("red", "blue", "yellow", "green"),
                             # alpha = 0.5,
                             size = 1.5,
                             axes = TRUE) +
  ggtitle("Neonate IFN and RSV infected")

splt_comp1_1



## To plot molecules over - it's not following my subset 
# obj_comp1_1@images$zoom <- CreateFOV(
#   coords = extracted_coord[, c("id", "x", "y")],
#   assay = "RNA",  # or whichever assay you're using
#   key = "zoom"
# )
# 
# DefaultFOV(obj_comp1_1) <- "zoom"

# splt_comp1_1 <- ImageDimPlot(subset(obj_comp1_1, cells = extracted_coord$id),
#                     group.by = 'insitutype_cluster_15_clusters', 
#                     border.color = "white", 
#                     fov = "zoom",
#                     border.size = 0.1,
#                     cols = cols, 
#                     molecules = rownames(obj_comp1_1[1:4]),
#                     mols.size = 0.5,
#                     mols.cols = c("red", "blue", "yellow", "green"),
#                     alpha = 0.5,
#                     size = 1.5,
#                     axes = TRUE) +
#   ggtitle("Neonate IFN and RSV infected")
# 
# splt_comp1_1

## ----



DefaultBoundary(obj_comp1_1[["FOV.2"]]) <- "centroids"
splt_comp1_1 <- ImageDimPlot(obj_comp1_1,
                             group.by = 'insitutype_cluster_15_clusters',
                             border.color = "white",
                             fov = "FOV.2",
                             border.size = 0.1,
                             cols = cols,
                             # molecules = rownames(obj_comp1_1[1:4]),
                             # mols.size = 0.5,
                             # mols.cols = c("red", "blue", "yellow", "green"),
                             # alpha = 0.5,
                             size = 1.5,
                             axes = F) +
  ggtitle("Neonate IFN and RSV infected") +
  NoLegend() +
  #theme_void() +  # Use theme_void() for a completely blank background
  theme(
    # plot.background = element_rect(fill = "white", color = NA),
    # panel.background = element_rect(fill = "white", color = NA),
    # plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    # legend.background = element_rect(fill = "white", color = NA),
    # legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  theme(axis.text = element_blank(),  # Remove axis text
        axis.title = element_blank(),  # Remove axis titles
        axis.ticks = element_blank())  # Remove axis ticks

splt_comp1_1


# Plot multiple TMAs 
obj_comp1_2 <- subset(obj, Group %in% c("Adult RSV infected"))

DefaultBoundary(obj_comp1_1[["FOV"]]) <- "centroids"
splt_comp1_2 <- ImageDimPlot(obj_comp1_2,
                             group.by = 'insitutype_cluster_15_clusters',
                             border.color = "white",
                             fov = "FOV",
                             border.size = 0.1,
                             cols = cols,
                             # molecules = rownames(obj_comp1_1[1:4]),
                             # mols.size = 0.5,
                             # mols.cols = c("red", "blue", "yellow", "green"),
                             # alpha = 0.5,
                             size = 0.8,
                             axes = F) +
  ggtitle("Adult RSV infected") +
  NoLegend() +
  #theme_void() +  # Use theme_void() for a completely blank background
  theme(
    # plot.background = element_rect(fill = "white", color = NA),
    # panel.background = element_rect(fill = "white", color = NA),
    # plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    # legend.background = element_rect(fill = "white", color = NA),
    # legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  theme(axis.text = element_blank(),  # Remove axis text
        axis.title = element_blank(),  # Remove axis titles
        axis.ticks = element_blank())  # Remove axis ticks

splt_comp1_2


# # Plot a TMA at a time 
# DefaultBoundary(obj_comp1_2[["FOV"]]) <- "centroids"
# tmas <- obj_comp1_2@meta.data$TMA %>% table() %>% names()
# 
# splt_comp1_2 <- ImageDimPlot(subset(obj_comp1_2, TMA %in% tmas[1]),
#                              group.by = 'insitutype_cluster_15_clusters',
#                              border.color = "white",
#                              fov = "FOV",
#                              border.size = 0.1,
#                              cols = cols,
#                              size = 1.5,
#                              axes = FALSE) +  # Set axes to FALSE
#   ggtitle("Adult RSV infected") +
#   NoLegend() +
#   #theme_void() +  # Use theme_void() for a completely blank background
#   theme(
#     # plot.background = element_rect(fill = "white", color = NA),
#     # panel.background = element_rect(fill = "white", color = NA),
#     # plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     # legend.background = element_rect(fill = "white", color = NA),
#     # legend.key = element_rect(fill = "white", color = NA),
#     legend.text = element_text(size = 10),
#     legend.title = element_text(size = 12, face = "bold")
#   ) +
#   theme(axis.text = element_blank(),  # Remove axis text
#         axis.title = element_blank(),  # Remove axis titles
#         axis.ticks = element_blank())  # Remove axis ticks
# 
# splt_comp1_2


## Plot - alternate clusters 
DefaultBoundary(obj_comp1_1[["FOV"]]) <- "centroids"
obj_comp1_2@meta.data$select_cluster  <- NA
obj_comp1_2@meta.data$select_cluster <- obj_comp1_2@meta.data$insitutype_cluster_15_clusters
obj_comp1_2@meta.data[!obj_comp1_2@meta.data$insitutype_cluster_15_clusters %in% "l", ]$select_cluster <- "other"

splt_comp1_2 <- ImageDimPlot(obj_comp1_2,
                             group.by = 'select_cluster',
                             border.color = "white",
                             fov = "FOV",
                             border.size = 0.1,
                             cols = "red",
                             # molecules = rownames(obj_comp1_1[1:4]),
                             # mols.size = 0.5,
                             # mols.cols = c("red", "blue", "yellow", "green"),
                             # alpha = 0.5,
                             size = 0.8,
                             axes = F) +
  ggtitle("Adult RSV infected") +
  NoLegend() +
  #theme_void() +  # Use theme_void() for a completely blank background
  theme(
    # plot.background = element_rect(fill = "white", color = NA),
    # panel.background = element_rect(fill = "white", color = NA),
    # plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    # legend.background = element_rect(fill = "white", color = NA),
    # legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  theme(axis.text = element_blank(),  # Remove axis text
        axis.title = element_blank(),  # Remove axis titles
        axis.ticks = element_blank())  # Remove axis ticks

splt_comp1_2








### ------ Plot Report ------ 
library(Seurat)
library(dplyr)
library(ggplot2)
library(dittoSeq)
library(viridis)
library(presto)


plot_report <- function(obj, group_1, group_2, ImageFeaturePlot_ptSize_g1, ImageFeaturePlot_ptSize_g2, ImageFeaturePlot_borSize, save_plot_outdir) {
  if (!is(obj, "Seurat")) {
    stop("Input object must be a Seurat object")
  }
  
  if (!all(c(group_1, group_2) %in% unique(obj$Group))) {
    stop("Specified groups not found in the object's metadata")
  }
  
  if (!dir.exists(save_plot_outdir)) {
    dir.create(save_plot_outdir, recursive = TRUE)
  }
  
  Start <- Sys.time()
  
  Groups <- c(group_1, group_2)
  
  # Ordering clusters
  obj@meta.data$insitutype_cluster_15_clusters <- factor(obj@meta.data$insitutype_cluster_15_clusters, levels = names(table(obj@meta.data$insitutype_cluster_15_clusters)))
  # Ordering groups
  obj@meta.data$Group <- factor(obj@meta.data$Group, levels = Groups)
  
  
  pdf(file = file.path(paste0(save_plot_outdir, group_1,"_vs_", group_2,"_report.pdf")), width = 12, height = 10)
  
  ## Set color map 
  cols <-
    c(
      '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
      '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
      '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
    )
  
  factor_levels <- obj$insitutype_cluster_15_clusters %>% table() %>% names()
  cols <- setNames(cols[1:length(factor_levels)], factor_levels)
  
  ## Subset objects 
  obj_comp1_1 <- subset(obj, Group %in% c(group_1))
  obj_comp1_2 <- subset(obj, Group %in% c(group_2))
  
  ## Plot: Clusters proportion by Group
  plt <- dittoBarPlot(subset(obj, Group %in% c(group_1, group_2)), 
                      var = "insitutype_cluster_15_clusters", 
                      group.by = "Group",
                      color.panel = cols,
                      main = paste0("Cluster Proportion By Group") )
  print(plt)
  
  ## Plot: All clusters over spatial plot - group_1
  plt <- ImageDimPlot(subset(obj_comp1_1, Group == group_1), 
                      group.by = 'insitutype_cluster_15_clusters', 
                      cols = cols, 
                      size = ImageFeaturePlot_ptSize_g1,
                      border.size = ImageFeaturePlot_borSize,
                      axes = F) +
    ggtitle(paste0("Group: ", group_1))
  print(plt)
  
  ## Plot: All clusters over spatial plot - group_2
  plt <- ImageDimPlot(subset(obj_comp1_2, Group == group_2), 
                      group.by = 'insitutype_cluster_15_clusters', 
                      cols = cols,
                      size = ImageFeaturePlot_ptSize_g2,
                      border.size = ImageFeaturePlot_borSize,
                      axes = F) +
    ggtitle(paste0("Group: ", group_2))
  print(plt)
  
  ## Plot: Clusters over spatial plot - one cluster at a time - group_1
  for (cluster in levels(obj$insitutype_cluster_15_clusters)) {
    plt <- ImageDimPlot(subset(obj_comp1_1, Group == group_1 & insitutype_cluster_15_clusters == cluster), 
                        group.by = 'insitutype_cluster_15_clusters', 
                        cols = cols[cluster],
                        size = ImageFeaturePlot_ptSize_g1,
                        border.size = ImageFeaturePlot_borSize,
                        axes = F) +
      ggtitle(paste0("Group: ", group_1, " - Cluster: ", cluster))
    print(plt)
  }
  
  ## Plot: Clusters over spatial plot - one cluster at a time - group_2
  for (cluster in levels(obj$insitutype_cluster_15_clusters)) {
    plt <- ImageDimPlot(subset(obj_comp1_2, Group == group_2 & insitutype_cluster_15_clusters == cluster), 
                        group.by = 'insitutype_cluster_15_clusters', 
                        cols = cols[cluster],
                        size = ImageFeaturePlot_ptSize_g2,
                        border.size = ImageFeaturePlot_borSize,
                        axes = F) +
      ggtitle(paste0("Group: ", group_2, " - Cluster: ", cluster))
    print(plt)
  }
  
  ## Run Find Markers 
  obj_sub <- subset(obj, Group %in% c(group_1, group_2))
  obj_sub@assays$RNA <- obj_sub@assays$Nanostring
  markers <- wilcoxauc(obj_sub, 'insitutype_cluster_15_clusters')
  
  top_markers <- markers %>%
    group_by(group) %>%
    dplyr::filter(logFC > 0 & padj <= 0.05) %>%
    dplyr::top_n(wt = logFC, n = 5) %>%
    ungroup() %>%
    pull(feature) %>%
    unique()
  
  ## Plot: Dotplot of Cluster Markers across clusters
  Idents(obj_sub) <- "insitutype_cluster_15_clusters"
  plt <- DotPlot(obj_sub, features = unique(c(top_markers))) +
    geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
    scale_colour_viridis(option = "magma") +
    guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle("Cluster Markers across clusters")
  print(plt)
  
  ## Plot: Cluster Markers across Groups
  Idents(obj_sub) <- "Group"
  plt <- DotPlot(obj_sub, features = unique(c(top_markers))) +
    geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
    scale_colour_viridis(option = "magma") +
    guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) + 
    ggtitle("Cluster Markers across groups")
  print(plt)
  
  ## High Specific Cluster makers on spatial plot 
  top_markers <- markers[order(markers$auc, decreasing = TRUE), ]$feature[1:20]
  top_markers <- top_markers[!top_markers %in% "Hba-a1/2"]
  
  ## Plot: Top markers Feature Spatial plot - group 1
  for (i in 1:length(top_markers)) {
    if (top_markers[i] %in% rownames(obj_comp1_1[["Nanostring"]])) {
      plt <- ImageFeaturePlot(subset(obj_comp1_1, Group == group_1), 
                              border.color = 'white', 
                              size = ImageFeaturePlot_ptSize_g1,
                              border.size = ImageFeaturePlot_borSize,
                              features = c(top_markers[i])) +
        ggtitle(paste0("Group: ", group_1, " - Top Cluster Markers ", top_markers[i]))
      print(plt)
    } else {
      warning(paste("Feature", top_markers[i], "not found. Skipping plot."))
    }
  }
  
  ## Plot: Top markers Feature Spatial plot - group 2
  for (i in 1:length(top_markers)) {
    if (top_markers[i] %in% rownames(obj_comp1_2[["Nanostring"]])) {
      plt <- ImageFeaturePlot(subset(obj_comp1_2, Group == group_2), 
                              border.color = 'white', 
                              size = ImageFeaturePlot_ptSize_g2,
                              border.size = ImageFeaturePlot_borSize,
                              features = c(top_markers[i])) +
        ggtitle(paste0("Group: ", group_2, " - Top Cluster Markers ", top_markers[i]))
      print(plt)
    } else {
      warning(paste("Feature", top_markers[i], "not found. Skipping plot."))
    }
  }
  
  dev.off()
  End <- Sys.time()
  plot_time <- print(End - Start)
  return(plot_time)
}


## Usage
#  Remove clusters e and f - they account just for a couple of cells for THIS comparison

obj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Objects/all_4_slides_processed.qs")

obj <- subset(obj, insitutype_cluster_15_clusters %in% c("e", "f"), invert = TRUE)

obj$insitutype_cluster_15_clusters <- droplevels(obj$insitutype_cluster_15_clusters[obj$insitutype_cluster_15_clusters != "e" & obj$insitutype_cluster_15_clusters != "f"])

table(obj$insitutype_cluster_15_clusters)
table(obj$Group)

# Comparison #1: Neonate IFN and RSV infected vs Adult RSV infected
plot_report(obj = obj, 
            group_1 = "Neonate IFN and RSV infected", 
            group_2 = "Adult RSV infected", 
            ImageFeaturePlot_ptSize_g1 = 2,
            ImageFeaturePlot_ptSize_g2 = 0.7,
            ImageFeaturePlot_borSize = 0.01,
            save_plot_outdir = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Plots/")
  

# Comparison #2: Neonate IFN and RSV reinfected vs Adult reinfected
plot_report(obj = obj, 
            group_1 = "Neonate IFN and RSV reinfected", 
            group_2 = "Adult reinfected", 
            ImageFeaturePlot_ptSize_g1 = 0.3,
            ImageFeaturePlot_ptSize_g2 = 0.7,
            ImageFeaturePlot_borSize = 0.01,
            save_plot_outdir = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Plots/")


# Comparison #3: Neonate RSV infected (NO IFN) vs Neonate IFN and RSV infected
plot_report(obj = obj, 
            group_1 = "Neonate RSV infected (NO IFN)", 
            group_2 = "Neonate IFN and RSV infected", 
            ImageFeaturePlot_ptSize_g1 = 0.7,
            ImageFeaturePlot_ptSize_g2 = 2,
            ImageFeaturePlot_borSize = 0.01,
            save_plot_outdir = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_3/Plots/")








## Testing view config (to get optimal ptSize) ------------------ 
## Comparison #1 g1 
group_1 = "Neonate IFN and RSV infected"
group_2 = "Adult RSV infected"

obj_comp1_1 <- subset(obj, Group %in% c(group_1))
obj_comp1_2 <- subset(obj, Group %in% c(group_2))

plt <- ImageDimPlot(subset(obj_comp1_1, Group == group_1), 
                    group.by = 'insitutype_cluster_15_clusters', 
                    cols = cols, 
                    size = 2,
                    border.size = 0.01,
                    axes = F) +
  ggtitle(paste0("Group: ", group_1))
print(plt)


## Comparison #1 g2
plt <- ImageDimPlot(subset(obj_comp1_2, Group == group_2), 
                    group.by = 'insitutype_cluster_15_clusters', 
                    cols = cols, 
                    size = 0.7,
                    border.size = 0.01,
                    axes = F) +
  ggtitle(paste0("Group: ", group_1))
print(plt)


## Comparison #2 g1
group_1 = "Neonate IFN and RSV reinfected"
group_2 = "Adult reinfected"

obj_comp2_1 <- subset(obj, Group %in% c(group_1))
obj_comp2_2 <- subset(obj, Group %in% c(group_2))

plt <- ImageDimPlot(subset(obj_comp2_1, Group == group_1), 
                    group.by = 'insitutype_cluster_15_clusters', 
                    cols = cols, 
                    size = 0.3,
                    border.size = 0.01,
                    axes = F) +
  ggtitle(paste0("Group: ", group_1))
print(plt)


## Comparison #2 g2 
plt <- ImageDimPlot(subset(obj_comp2_2, Group == group_2), 
                    group.by = 'insitutype_cluster_15_clusters', 
                    cols = cols, 
                    size = 0.7,
                    border.size = 0.01,
                    axes = F) +
  ggtitle(paste0("Group: ", group_2))
print(plt)



## Comparison #3 g1 
group_1 = "Neonate RSV infected (NO IFN)"
group_2 = "Neonate IFN and RSV infected"

obj_comp3_1 <- subset(obj, Group %in% c(group_1))
obj_comp3_2 <- subset(obj, Group %in% c(group_2))

plt <- ImageDimPlot(subset(obj_comp3_1, Group == group_1), 
                    group.by = 'insitutype_cluster_15_clusters', 
                    cols = cols, 
                    size = 0.7,
                    border.size = 0.01,
                    axes = F) +
  ggtitle(paste0("Group: ", group_1))
print(plt)

## Comparison #3 g2 
plt <- ImageDimPlot(subset(obj_comp3_2, Group == group_2), 
                    group.by = 'insitutype_cluster_15_clusters', 
                    cols = cols, 
                    size = 2,
                    border.size = 0.01,
                    axes = F) +
  ggtitle(paste0("Group: ", group_2))
print(plt)



















