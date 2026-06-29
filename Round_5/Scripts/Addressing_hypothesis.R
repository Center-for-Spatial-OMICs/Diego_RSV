## ====== Date: Aug 14, 2025  ===== 


### Load and define functions ---------

library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(qs)



All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_processed.qs")

comparisons <- list(
  "Comparison" = list(
    group_1 = "Neonate RSV infected (NO IFN)",
    group_2 = "Neonate IFN and RSV infected",
    group_3 = "Adult RSV infected"
  )
)

obj <- All_Sobj


p <- dittoBarPlot(
  obj,
  var = "major_celltype",
  group.by = "Group",
  #color.panel = ,
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)
print(p)


obj_sub_adu <- subset(obj, Group %in% c("Adult Control",
                                    "Adult RSV infected"))

obj_sub_adu@meta.data$Group <- factor(obj_sub_adu@meta.data$Group, 
                                  levels = c("Adult Control",
                                             "Adult RSV infected"))

obj_sub_neo <- subset(obj, Group %in% c("Neonate control",
                                        "Neonate RSV infected (NO IFN)",
                                        "Neonate IFN and RSV infected"))
obj_sub_neo@meta.data$Group <- factor(obj_sub_neo@meta.data$Group, 
                                  levels = c("Neonate control",
                                             "Neonate RSV infected (NO IFN)",
                                             "Neonate IFN and RSV infected"))

p1 <- dittoBarPlot(
  obj_sub_adu,
  var = "major_celltype",
  group.by = "Group",
  #color.panel = ,
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)
print(p1)


p2 <- dittoBarPlot(
  obj_sub_neo,
  var = "major_celltype",
  group.by = "Group",
  #color.panel = ,
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)
print(p2)

fills <- ggplot_build(p2)$data[[1]]$fill
col_palette <- unique(fills)
col_palette <- setNames(
  c("#E69F00", "#56B4E9", "#009E73"),
  c(table(obj_sub_neo$major_celltype) %>% names())
)


# Looking for "lower frequency of B cells in Nenoates compared to adults" 
Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_lym_Neocy_Granu_to_minor.qs") 


# Adding lymphocytes minor cell type label 
obj_sub <- subset(Sobj, lym_nuo_granu_recluster %in% "Bcell")

obj_sub_adu@meta.data[obj_sub_adu@meta.data$cell %in% obj_sub@meta.data$cell, ]$major_celltype <- "Bcell"

obj_sub_neo@meta.data[obj_sub_neo@meta.data$cell %in% obj_sub@meta.data$cell, ]$major_celltype <- "Bcell"


obj_sub <- subset(Sobj, lym_nuo_granu_recluster %in% "Tcd4")

obj_sub_adu@meta.data[obj_sub_adu@meta.data$cell %in% obj_sub@meta.data$cell, ]$major_celltype <- "Tcd4"

obj_sub_neo@meta.data[obj_sub_neo@meta.data$cell %in% obj_sub@meta.data$cell, ]$major_celltype <- "Tcd4"


obj_sub <- subset(Sobj, lym_nuo_granu_recluster %in% "Tcd8")

obj_sub_adu@meta.data[obj_sub_adu@meta.data$cell %in% obj_sub@meta.data$cell, ]$major_celltype <- "Tcd8"

obj_sub_neo@meta.data[obj_sub_neo@meta.data$cell %in% obj_sub@meta.data$cell, ]$major_celltype <- "Tcd8"


# col_palette <- setNames(
#   c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#1C91D2","#0072B2", "#D55E00", "#CC79A7",
#     "#666666", "#AD7700", "#007756", "#D5C711", "#005685", "#1C91D4", "#56B4E8"),
#   c(table(obj_sub_adu$major_celltype) %>% names())
# )


obj_sub_adu_lym <- subset(obj_sub_adu, major_celltype %in% c("Tcd4", "Tcd8", "Bcell"))
p1 <- dittoBarPlot(
  obj_sub_adu_lym,
  var = "major_celltype",
  group.by = "Group",
  color.panel = c("#F0E442", "#E69F00", "#D55E00"),
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)
print(p1)

obj_sub_neo_lym <- subset(obj_sub_neo, major_celltype %in% c("Tcd4", "Tcd8", "Bcell"))
p2 <- dittoBarPlot(
  obj_sub_neo_lym,
  var = "major_celltype",
  group.by = "Group",
  color.panel = c("#F0E442", "#E69F00", "#D55E00"),
  retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)
print(p2)









fills <- ggplot_build(p)$data[[1]]$fill
col_palette <- unique(fills)

col_palette <- setNames(
  c("#E69F00", "#56B4E9", "#009E73"),
  c("9_res_0.5",
    "8_res_0.5",
    "7_res_0.5")
)

## All groups 
p <- dittoBarPlot(
  obj,
  var = "major_celltype",
  group.by = "Group",
  color.panel = col_palette,
  #retain.factor.levels = TRUE, # Ensure factor levels are respected
  main = ""
)
print(p)

plot_and_data <- dittoBarPlot(
  obj,
  var = "major_celltype",
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
