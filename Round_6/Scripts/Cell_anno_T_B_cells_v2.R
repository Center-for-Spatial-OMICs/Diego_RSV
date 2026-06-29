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

obj_list$`NA` <- NULL
obj_list$`NA.1` <- NULL


## Diego's cell markers  ---------
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




## Hierarchical annotation based on transcript counts ---------

# Step1: annotate minor Lymphocytes 
# Loop over each object in the list
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

# Step2: get a visualization 
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
df_long$Group <- factor(df_long$Group, levels = desired_order)

# Plot: Stacked barplot, colored by cell type, faceted by group
ggplot(
  #df_long[df_long$Group %in% c("Adult Control", "Adult RSV infected"), ],
  df_long[df_long$Group %in% c("Neonate control", "Neonate RSV infected (NO IFN)", "Neonate IFN and RSV infected"), ],
  aes(x = Object, y = Count, fill = CellType)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group, scales = "free_x") +
  labs(title = "Stacked Barplot: Cell Type Annotation by Group", y = "Cell Count") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())



# Filter to desired groups
df_prop <- df_long %>%
  filter(Group %in% desired_order) %>%
  group_by(Object) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

ggplot(
  #df_prop[df_prop$Group %in% c("Adult Control", "Adult RSV infected"), ],
  df_prop[df_prop$Group %in% c("Neonate control", "Neonate RSV infected (NO IFN)", "Neonate IFN and RSV infected"), ],
  aes(x = Object, y = Count, fill = CellType)
) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill" makes it proportional
  facet_wrap(~ Group, scales = "free_x") +
  labs(title = "Proportion Plot: Cell Type Annotation by Group",
       y = "Proportion of Cell Types") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())




# Load required libraries
library(dplyr)
library(ggplot2)
library(FSA) # for Dunn's test

celltype_of_interest <- "Tcd4" #"Bcell" ##"Tcd8"

# 1. For each core (object), calculate % of cell type of interest
df_pct <- df_long %>%
  filter(Group %in% desired_order) %>%
  group_by(Object, Group) %>%
  mutate(TotalCells = sum(Count)) %>%
  ungroup() %>%
  filter(CellType == celltype_of_interest) %>%
  mutate(Proportion = Count / TotalCells) %>%
  select(Object, Group, Proportion)

# 2. Boxplot: Compare proportions across groups
ggplot(df_pct, aes(x = Group, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  labs(title = paste("Proportion of", celltype_of_interest, "Across Groups"),
       y = paste("Proportion of", celltype_of_interest),
       x = "Group") +
  theme_minimal()

# 3. Kruskal-Wallis test
kw_test <- kruskal.test(Proportion ~ Group, data = df_pct)
print(kw_test)

# 4. Dunn's post-hoc test (if Kruskal-Wallis is significant)
if (kw_test$p.value < 0.05) {
  dunn_test <- dunnTest(Proportion ~ Group, data = df_pct, method = "bonferroni")
  print(dunn_test)
}


# now plot it by TMAs 
obj <- obj_list[[1]]
DefaultAssay(obj) <- "Nanostring"
Idents(obj) <- "lyn_minor_celltype"
DotPlot(obj,
        features =  unique(c(cd4_markers, cd8_markers, bcell_markers))
        #features =  unique(c(top_markers_auc))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("")


DefaultAssay(obj) <- "Nanostring"
Idents(obj) <- "major_celltype"
DotPlot(obj,
        features =  unique(c(cd4_markers, cd8_markers, bcell_markers))
        #features =  unique(c(top_markers_auc))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("")




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
  labs(title = "Cell Type Proportions per Group",
       y = "Proportion",
       x = "Group") +
  theme_minimal()

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









## tests --------



# Load required libraries
library(dplyr)
library(ggplot2)
library(FSA) # for Dunn's test


# 1. Prepare proportion data per object and cell type
df_prop_summary <- df_long %>%
  filter(Group %in% desired_order) %>%
  group_by(Object, Group, CellType) %>%
  summarise(Proportion = Count / sum(Count), .groups = "drop")

# 2. Kruskal-Wallis test for each cell type
cell_types_to_test <- unique(df_prop_summary$CellType)
kruskal_results <- lapply(cell_types_to_test, function(ct) {
  data_ct <- df_prop_summary %>% filter(CellType == ct)
  test <- kruskal.test(Proportion ~ Group, data = data_ct)
  data.frame(
    CellType = ct,
    Statistic = test$statistic,
    p.value = test$p.value
  )
})
kruskal_df <- bind_rows(kruskal_results)
print("Kruskal-Wallis test results:")
print(kruskal_df)

# 3. Dunn's post-hoc test for each cell type (if needed)
dunn_results <- lapply(cell_types_to_test, function(ct) {
  data_ct <- df_prop_summary %>% filter(CellType == ct)
  if(length(unique(data_ct$Group)) > 1) {
    test <- dunnTest(Proportion ~ Group, data = data_ct, method = "bonferroni")
    out <- as.data.frame(test$res)
    out$CellType <- ct
    out
  } else {
    data.frame(Comparison = NA, Z = NA, P.unadj = NA, P.adj = NA, CellType = ct)
  }
})
dunn_df <- bind_rows(dunn_results)
print("Dunn's post-hoc test results:")
print(dunn_df)

# 4. Boxplot visualization
ggplot(df_prop_summary, aes(x = Group, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ CellType, scales = "free_y") +
  labs(title = "Cell Type Proportion Across Groups",
       y = "Proportion",
       x = "Group") +
  theme_minimal()












# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

# 1. Collect annotation tables from all objects
celltype_tables <- lapply(obj_list, function(obj) table(obj@meta.data$lyn_minor_celltype))

# 2. Convert tables to a long-format data frame
df <- bind_rows(
  lapply(seq_along(celltype_tables), function(i) {
    tibble(
      Object = paste0("Obj", i),
      CellType = names(celltype_tables[[i]]),
      Count = as.numeric(celltype_tables[[i]])
    )
  })
)

# 3. Fill missing cell types with zero (wide format for summary/heatmap)
df_wide <- df %>%
  tidyr::pivot_wider(names_from = CellType, values_from = Count, values_fill = 0)

# 4. Stacked Barplot: Cell type counts per object
df_long <- df_wide %>%
  tidyr::pivot_longer(-Object, names_to = "CellType", values_to = "Count")

df_long <- df_long[!df_long$CellType %in% "other", ]

ggplot(df_long, aes(x = Object, y = Count, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(title = "Stacked Barplot: Cell Type Annotation Across All Objects", y = "Cell Count") +
  theme_minimal()






# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

# 1. Count unique cell types per object
n_types <- sapply(obj_list, function(obj) length(unique(obj@meta.data$lyn_minor_celltype)))

# 2. Create a data frame for plotting
df_types <- data.frame(
  Object = paste0("Obj", seq_along(obj_list)),
  NumCellTypes = n_types
)

# 3. Barplot: Number of cell types per object
ggplot(df_types, aes(x = Object, y = NumCellTypes)) +
  geom_bar(stat = "identity", fill = "#377eb8") +
  labs(title = "Number of Cell Types per Object", y = "Unique Cell Types") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# 4. Histogram: Distribution of cell type counts across objects
ggplot(df_types, aes(x = NumCellTypes)) +
  geom_histogram(binwidth = 1, fill = "#e41a1c", color = "black") +
  labs(title = "Distribution of Cell Type Counts Across Objects",
       x = "Number of Unique Cell Types",
       y = "Number of Objects") +
  theme_minimal()

# 5. Summary table: How many objects have each possible number of cell types
print(table(df_types$NumCellTypes))

# 6. Top 10 most and least diverse objects
print("Top 10 most diverse objects:")
print(head(df_types[order(-df_types$NumCellTypes), ], 10))

print("Top 10 least diverse objects:")
print(head(df_types[order(df_types$NumCellTypes), ], 10))

# 7. Optional: Presence/absence heatmap of cell types across objects
# First, get all possible cell types
all_celltypes <- unique(unlist(lapply(obj_list, function(obj) unique(obj@meta.data$lyn_minor_celltype))))

# Build a presence/absence matrix
presence_matrix <- sapply(obj_list, function(obj) {
  ct <- unique(obj@meta.data$lyn_minor_celltype)
  as.integer(all_celltypes %in% ct)
})
colnames(presence_matrix) <- paste0("Obj", seq_along(obj_list))
rownames(presence_matrix) <- all_celltypes

# Plot heatmap
pheatmap(presence_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "#4daf4a"),
         main = "Presence/Absence of Cell Types Across Objects")















obj <- obj_list[[1]]

cd4_markers <- intersect(cell_types$Lymphocytes$T_CD4_Cells$Markers, obj@assays$Nanostring$counts %>% rownames())
cd8_markers <- intersect(cell_types$Lymphocytes$T_CD8_Cells$Markers, obj@assays$Nanostring$counts %>% rownames())
bcell_markers <- intersect(cell_types$Lymphocytes$Activated_B_Cells$Markers, obj@assays$Nanostring$counts %>% rownames())


counts <- obj@assays$Nanostring$counts

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
obj@meta.data$lyn_minor_celltype <- annot
table(obj@meta.data$lyn_minor_celltype)

# plot a barplot for table(obj@meta.data$lyn_minor_celltype)














# now plot it by TMAs 
obj <- NormalizeData(obj, assay = "Nanostring")
DefaultAssay(obj) <- "Nanostring"
Idents(obj) <- "lyn_minor_celltype"
DotPlot(obj,
        features =  unique(c(cd4_markers, cd8_markers, bcell_markers))
        #features =  unique(c(top_markers_auc))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("")

dittoBarPlot(
  subset(obj, lyn_minor_celltype %in% "other", invert = TRUE),
  var = "lyn_minor_celltype",
  group.by = "TMA_4",
  #color.panel = minor_celltype_concate_colors,
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)









