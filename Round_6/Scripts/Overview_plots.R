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
library(patchwork)


## Load objects -----------
# ScObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/ScObj_umap.qs")

obj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/Sobj_list_TMA_ID_fixed.qs")
obj_list$`NA` <- NULL
obj_list$`NA.1` <- NULL

## Plot QC metris spatially 
library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)

# Apply DefaultBoundary setting
obj_list <- lapply(obj_list, function(x) {
  DefaultBoundary(x[["FOV"]]) <- "centroids"
  return(x)
})

# Get unique groups from the combined metadata of all objects
groups <- unique(sapply(obj_list, function(x) x@meta.data$Group[1]))

plot_list <- list()

for (grp in groups) {
  # Filter objects that belong to this group
  objs_in_group <- obj_list[sapply(obj_list, function(x) x@meta.data$Group[1] == grp)]
  
  # For each object in the group, create a plot for the specified features
  plots_per_obj <- lapply(objs_in_group, function(obj) {
    ImageFeaturePlot(obj, features = c('nCount_RNA', 'nFeature_RNA'), max.cutoff = c("q99")) &
      scale_fill_viridis() &
      coord_flip() &
      ggtitle(obj@meta.data$Group[1], subtitle = obj@meta.data$TMA_4[1]) &
      theme(plot.title = element_text(size = 14, hjust = 0))
  })
  
  # Store the list of plots per group
  plot_list[[grp]] <- plots_per_obj
}

# To view plots for a particular group, e.g.:
# pdf gets a weird pixel-color resolution; png doens't worki neither
# png(filename = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Plots/spatialPlot_nfeature_ncount_gp1.png",
#     width = 3000, height = 1000, res = 300)
# 
# print(wrap_plots(plot_list[[1]], ncol = 5))
# 
# dev.off()

print(wrap_plots(plot_list[[1]], ncol = 3))

print(wrap_plots(plot_list[[2]], ncol = 3))

print(wrap_plots(plot_list[[3]], ncol = 3))



print(wrap_plots(plot_list[[4]], ncol = 4))

print(wrap_plots(plot_list[[5]], ncol = 3))

print(wrap_plots(plot_list[[6]], ncol = 4))

print(wrap_plots(plot_list[[7]], ncol = 4))

print(wrap_plots(plot_list[[8]], ncol = 4))



## Plot QC metrics vlnplot ~ TMA ~ Group 
all_meta <- lapply(obj_list, function(obj) {
  meta <- obj@meta.data
  return(meta)
})
all_meta <- do.call(plyr::rbind.fill, all_meta)

# nFeature_RNA ~ gropups 
pdf(file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Plots/nfeature_facet_groups.pdf", width = 18, height = 10)

plt <- ggplot(all_meta, aes(x = nFeature_RNA, y = TMA_4, fill = TMA_3)) +
  geom_violin() +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 3,
    color = "black",
    fill = "white"
  ) +
  xlim(0, 500) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),            # smaller y tick labels
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),  # smaller y title
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, color = "black", face = "bold"),
    plot.margin = margin(10, 10, 10, 40),
    panel.spacing = unit(0.5, "lines")
  ) +
  NoLegend() +
  facet_grid(~ Group)

print(plt)

dev.off()



# nCount_RNA ~ gropups 
pdf(file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Plots/ncount_facet_groups.pdf", width = 18, height = 10)

plt <- ggplot(all_meta, aes(x = nCount_RNA, y = TMA_4, fill = TMA_3)) +
  geom_violin() +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 3,
    color = "black",
    fill = "white"
  ) +
  xlim(0, 500) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),            # smaller y tick labels
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),  # smaller y title
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, color = "black", face = "bold"),
    plot.margin = margin(10, 10, 10, 40),
    panel.spacing = unit(0.5, "lines")
  ) +
  NoLegend() +
  facet_grid(~ Group)

print(plt)

dev.off()


# ncount vs nfeature scatter plot
pdf(file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Plots/ncount_vs_nfeature.pdf", width = 36, height = 6)

library(ggplot2)
plt <- ggplot(all_meta, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_bin2d(bins = 100) +  # Adjust bins for resolution
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  facet_grid(~ Group) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 40, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 40),            # smaller y tick labels
    axis.title.x = element_text(size = 40, face = "bold"),
    axis.title.y = element_text(size = 40, face = "bold"),  # smaller y title
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 22, color = "black", face = "bold"),
    plot.margin = margin(10, 10, 10, 40),
    panel.spacing = unit(0.5, "lines"),
    strip.text.x = element_text(size = 18, face = "bold")
  ) 
 

print(plt)

dev.off()

## Plot cell type prop  barplot ~ TMA ~ Group
all_meta <- lapply(obj_list, function(obj) {
  meta <- obj@meta.data
  return(meta)
})
all_meta <- do.call(plyr::rbind.fill, all_meta)

cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- all_meta$major_celltype %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)

desired_order <- c(
  "Adult Control",
  "Adult RSV infected",
  "Adult reinfected",
  "Neonate control",
  "Neonate RSV infected (NO IFN)",
  "Neonate (NO IFN) reinfected",
  "Neonate IFN and RSV infected",
  "Neonate IFN and RSV reinfected"
)

library(ggplot2)
library(patchwork)
plot_list <- lapply(desired_order, function(grp) {
  all_meta_sub <- all_meta[all_meta$Group == grp, ]
  ggplot(all_meta_sub, aes(x = TMA_4, fill = major_celltype)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = celltype_cols) +
    labs(
      x = "",
      y = "Proportion of Cell Types",
      fill = "Cell Type",
      title = paste("Cell Type Proportions by TMA Core\nGroup:", grp)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16, color = "black"),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 14, color = "black"),
      axis.text.x = element_text(face = "plain", size = 10, color = "black", angle = 90),
      axis.text.y = element_text(face = "plain", size = 12, color = "black"),
      legend.title = element_text(face = "bold", size = 14, color = "black"),
      legend.text = element_text(face = "plain", size = 12, color = "black")
    )
})

pdf(file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Plots/celltype_prop_barplots.pdf", width = 8, height = 6)
for(i in 1:length(plot_list)){
print(plot_list[[i]])
}
dev.off()



for (i in seq_along(obj_list)) {
  obj <- obj_list[[i]]
  
  # Find marker genes present in the data
  cd4_markers <- intersect(cell_types$Lymphocytes$T_CD4_Cells$Markers, rownames(obj@assays$Nanostring$counts))
  cd8_markers <- intersect(cell_types$Lymphocytes$T_CD8_Cells$Markers, rownames(obj@assays$Nanostring$counts))
  bcell_markers <- intersect(cell_types$Lymphocytes$Activated_B_Cells$Markers, rownames(obj@assays$Nanostring$counts))
  
  counts <- obj@assays$Nanostring$counts
  
  # Initialize annotation vector with default "other"
  annot <- rep("other", ncol(counts))
  
  # Hierarchical annotation based on transcript counts
  has_Ptprc <- counts["Ptprc", ] > 0
  has_Cd3e <- counts["Cd3e", ] > 0
  has_Cd4 <- counts["Cd4", ] > 0
  has_Cd8a <- counts["Cd8a", ] > 0
  has_Cd19 <- counts["Cd19", ] > 0
  
  # Assign Tcd4
  annot[has_Ptprc & has_Cd3e & has_Cd4] <- "Tcd4"
  # Assign Tcd8
  annot[has_Ptprc & has_Cd3e & has_Cd8a] <- "Tcd8"
  # Assign B cell
  annot[has_Ptprc & !has_Cd3e & has_Cd19] <- "Bcell"
  
  # Add annotation to Seurat meta data
  obj@meta.data$lyn_minor_celltype <- annot
  
  # Show table and plot
  cat("Object", i, "\n")
  print(table(obj@meta.data$lyn_minor_celltype))
  barplot(table(obj@meta.data$lyn_minor_celltype), main=paste("Cell type annotation for obj", i), ylab="Count")
  
  # Save back to obj_list (if you want to keep the annotation)
  obj_list[[i]] <- obj
}


all_meta <- lapply(obj_list, function(obj) {
  meta <- obj@meta.data
  return(meta)
})
all_meta <- do.call(plyr::rbind.fill, all_meta)

all_meta <- all_meta %>% filter(lyn_minor_celltype != "other")


desired_order <- c(
  "Adult Control",
  "Adult RSV infected",
  "Adult reinfected",
  "Neonate control",
  "Neonate RSV infected (NO IFN)",
  "Neonate (NO IFN) reinfected",
  "Neonate IFN and RSV infected",
  "Neonate IFN and RSV reinfected"
)

library(ggplot2)
library(patchwork)
plot_list <- lapply(desired_order, function(grp) {
  all_meta_sub <- all_meta[all_meta$Group == grp, ]
  ggplot(all_meta_sub, aes(x = TMA_4, fill = lyn_minor_celltype)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = c("#F0E442", "#E69F00", "#D55E00")) +
    labs(
      x = "",
      y = "Proportion of Cell Types",
      fill = "Cell Type",
      title = paste("Cell Type Proportions by TMA Core\nGroup:", grp)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16, color = "black"),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 14, color = "black"),
      axis.text.x = element_text(face = "plain", size = 10, color = "black", angle = 90),
      axis.text.y = element_text(face = "plain", size = 12, color = "black"),
      legend.title = element_text(face = "bold", size = 14, color = "black"),
      legend.text = element_text(face = "plain", size = 12, color = "black")
    )
})

print(plot_list[[1]])


pdf(file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Plots/TcellBcell_prop_barplots.pdf", width = 8, height = 6)
for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}
dev.off()




## Plot cell type prop ∼ RSV infection effect ------------
if (!"lyn_minor_celltype" %in% colnames(obj_list[[1]]@meta.data)) {
  stop("Column 'lyn_minor_celltype' not found in obj_list[[1]]. Execution halted.")
}

## major cell types -------
# Collect cell type tables and group info
celltype_tables <- lapply(obj_list, function(obj) table(obj@meta.data$major_celltype))
group_vec <- sapply(obj_list, function(obj) unique(obj@meta.data$Group)[1])

# Convert tables to long-format data frame with group info
df <- bind_rows(
  lapply(seq_along(celltype_tables), function(i) {
    tibble(
      Object = paste0("Obj", i),
      Group = group_vec[i],
      CellType = names(celltype_tables[[i]]),
      Count = as.numeric(celltype_tables[[i]])
    )
  })
)

# Remove "other" cell types
df_long <- df %>% filter(CellType != "other")

# Wide format (if needed)
df_wide <- df_long %>%
  tidyr::pivot_wider(names_from = CellType, values_from = Count, values_fill = 0)


desired_order <- c(
  "Adult Control",
  "Adult RSV infected",
  "Adult reinfected",
  "Neonate control",
  "Neonate RSV infected (NO IFN)",
  "Neonate (NO IFN) reinfected",
  "Neonate IFN and RSV infected",
  "Neonate IFN and RSV reinfected"
)


celltypes_of_interest <- c(table(obj_list[[1]]@meta.data$major_celltype) %>% names())


# Make sure Group is a factor with correct order
df_long$Group <- factor(df_long$Group, levels = desired_order)
df_long$CellType <- factor(df_long$CellType, levels = celltypes_of_interest)

# 1. For each core (object), calculate % of each cell type of interest
df_pct <- df_long %>%
  filter(Group %in% desired_order, CellType %in% celltypes_of_interest) %>%
  group_by(Object, Group) %>%
  mutate(TotalCells = sum(Count)) %>%
  ungroup() %>%
  mutate(Proportion = Count / TotalCells) %>%
  select(Object, Group, CellType, Proportion)

# 2. Grouped boxplot: Each group, boxplots for each cell type (side-by-side)
ggplot(df_pct, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = CellType), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              alpha = 0.5, size = 1) +
  scale_fill_manual(values = celltype_cols) +
  labs(title = "Cell Type Proportions per Group (B and T cells)",
       y = "Cell Type Proportion",
       x = "Group") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "black"),
    axis.title.x = element_text(face = "plain", size = 16),
    axis.title.y = element_text(face = "bold", size = 14, color = "black"),
    axis.text.x = element_text(face = "bold", size = 18, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "plain", size = 14, color = "black"),
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text = element_text(face = "plain", size = 12, color = "black")
  )


## B and T cells -------
# Collect cell type tables and group info
celltype_tables <- lapply(obj_list, function(obj) table(obj@meta.data$lyn_minor_celltype))
group_vec <- sapply(obj_list, function(obj) unique(obj@meta.data$Group)[1])

# Convert tables to long-format data frame with group info
df <- bind_rows(
  lapply(seq_along(celltype_tables), function(i) {
    tibble(
      Object = paste0("Obj", i),
      Group = group_vec[i],
      CellType = names(celltype_tables[[i]]),
      Count = as.numeric(celltype_tables[[i]])
    )
  })
)

# Remove "other" cell types
df_long <- df %>% filter(CellType != "other")

# Wide format (if needed)
df_wide <- df_long %>%
  tidyr::pivot_wider(names_from = CellType, values_from = Count, values_fill = 0)


desired_order <- c(
  "Adult Control",
  "Adult RSV infected",
  "Adult reinfected",
  "Neonate control",
  "Neonate RSV infected (NO IFN)",
  "Neonate (NO IFN) reinfected",
  "Neonate IFN and RSV infected",
  "Neonate IFN and RSV reinfected"
)


celltypes_of_interest <- c("Bcell", "Tcd8", "Tcd4" )


# Make sure Group is a factor with correct order
df_long$Group <- factor(df_long$Group, levels = desired_order)
df_long$CellType <- factor(df_long$CellType, levels = celltypes_of_interest)

# 1. For each core (object), calculate % of each cell type of interest
df_pct <- df_long %>%
  filter(Group %in% desired_order, CellType %in% celltypes_of_interest) %>%
  group_by(Object, Group) %>%
  mutate(TotalCells = sum(Count)) %>%
  ungroup() %>%
  mutate(Proportion = Count / TotalCells) %>%
  select(Object, Group, CellType, Proportion)

# 2. Grouped boxplot: Each group, boxplots for each cell type (side-by-side)
ggplot(df_pct, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = CellType), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("#F0E442", "#E69F00", "#D55E00")) +
  labs(title = "Cell Type Proportions per Group (B and T cells)",
       y = "Cell Type Proportion",
       x = "Group") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "black"),
    axis.title.x = element_text(face = "plain", size = 16),
    axis.title.y = element_text(face = "bold", size = 14, color = "black"),
    axis.text.x = element_text(face = "bold", size = 18, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(face = "plain", size = 14, color = "black"),
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text = element_text(face = "plain", size = 12, color = "black")
  )

# 3. Kruskal-Wallis test for each cell type (across groups)
kruskal_results <- lapply(celltypes_of_interest, function(ct) {
  df_ct <- df_pct %>% filter(CellType == ct)
  test <- kruskal.test(Proportion ~ Group, data = df_ct)
  data.frame(CellType = ct, Statistic = test$statistic, p.value = test$p.value)
})
kruskal_df <- bind_rows(kruskal_results)
print("Kruskal-Wallis test results:")
print(kruskal_df)

# 4. Dunn's post-hoc test for each cell type (if Kruskal-Wallis is significant)
dunn_results <- lapply(celltypes_of_interest, function(ct) {
  df_ct <- df_pct %>% filter(CellType == ct)
  kw_test <- kruskal.test(Proportion ~ Group, data = df_ct)
  if (kw_test$p.value < 0.05) {
    test <- dunnTest(Proportion ~ Group, data = df_ct, method = "bonferroni")
    out <- as.data.frame(test$res)
    out$CellType <- ct
    out
  } else {
    data.frame(Comparison = NA, Z = NA, P.unadj = NA, P.adj = NA, CellType = ct)
  }
})
dunn_df <- bind_rows(dunn_results)
print("Dunn's post-hoc test results (if Kruskal-Wallis significant):")
print(dunn_df)




