## ====== Date: Sep 11, 2025  ===== 


### Load and define functions ---------

library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(qs)

options(future.globals.maxSize = 99999 * 1024^2)
Sys.setenv("VROOM_CONNECTION_SIZE" = 999999 * 2)

### Load and define variables ---------
Sobj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj_list_celltype_niche.qs") #objects are separated by slides, not TMAs. Although I have a TMA columns to retrieve this info on their @meta.data 


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


cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- Sobj_list[[1]]$major_celltype %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)

library(viridis)
niche_cols <- viridis(5)
factor_levels <- Sobj_list[[1]]$niches_5 %>% table() %>% names()
niche_cols <- setNames(niche_cols[1:length(factor_levels)], factor_levels)





obj <- subset(Sobj_list[[1]], TMA %in% "Ap17-0333 1")
ImageFeaturePlot(obj, 
                 features = c("Ptprc", "Cd19", "Cd69",
                              "Tnfsf13b", "Tnfsf13",
                              "Tnfrsf13c", "Tnfrsf13b", "Tnfrsf17",
                              "Cd27",
                              "Sdc1")) &
  scale_fill_viridis() &
  #coord_flip() &
  ggtitle(obj@meta.data$Group[1]) &
  theme(plot.title = element_text(size = 14, hjust = 0))


obj@meta.data$highlight_anno <- "other"
obj@meta.data[obj@meta.data$major_celltype %in% "Lymphocyte", ]$highlight_anno <- "Lymphocyte"

ImageDimPlot(obj,
             group.by = 'highlight_anno',
             border.color = "white",
             fov = "FOV",
             border.size = 0.1,
             #cols = celltype_cols,
             # molecules = rownames(obj_comp1_1[1:4]),
             # mols.size = 0.5,
             # mols.cols = c("red", "blue", "yellow", "green"),
             # alpha = 0.5,
             size = 2,
             axes = F) +
  ggtitle("")


minor_lym <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_lym_Neocy_Granu_to_minor.qs")



intersect(obj@meta.data$cell_id, minor_lym@meta.data$cell_id) #yes
intersect(obj@meta.data %>% rownames(), minor_lym@meta.data$cell_id) #no
rownames(obj@meta.data) <- obj@meta.data$cell_id
intersect(obj@meta.data %>% rownames(), minor_lym@meta.data$cell_id) #yes 


Bcell_ID <- minor_lym@meta.data[minor_lym@meta.data$lym_nuo_granu_recluster %in% "Bcell", ]$cell_id

# obj_bcell <- subset(obj, cells = Bcell_ID) #errors ...


obj@meta.data$highlight_anno <- "other"
obj@meta.data[obj@meta.data$cell_id %in% Bcell_ID, ]$highlight_anno <- "Bcell"

ImageDimPlot(obj,
             group.by = 'highlight_anno',
             border.color = "white",
             fov = "FOV",
             border.size = 0.1,
             cols = c("gray", "gold"),
             # molecules = rownames(obj_comp1_1[1:4]),
             # mols.size = 0.5,
             # mols.cols = c("red", "blue", "yellow", "green"),
             # alpha = 0.5,
             size = 1.5,
             axes = F) +
  ggtitle("")




table(obj$major_celltype, 
      obj$highlight_anno)


## THIS IS THE WRITE WAY TO COME FROM THE SC-OBJECTS BACK TO THE "ORGINAL-OBJECTS"


TMAs_1 <- unique(Sobj_list[[i]]$TMA)
obj <- subset(Sobj_list[[i]], TMA %in% TMAs_1[j])

Bcell_ID <- minor_lym@meta.data[minor_lym@meta.data$lym_nuo_granu_recluster %in% "Bcell" & minor_lym@meta.data$TMA %in% TMAs_1[j], ]$cell_id

obj@meta.data$highlight_anno <- "other"
obj@meta.data[obj@meta.data$cell_id %in% Bcell_ID, ]$highlight_anno <- "Bcell"

ImageDimPlot(obj,
             group.by = 'highlight_anno',
             border.color = "white",
             fov = "FOV",
             border.size = 0.1,
             cols = c("gray", "gold"),
             # molecules = rownames(obj_comp1_1[1:4]),
             # mols.size = 0.5,
             # mols.cols = c("red", "blue", "yellow", "green"),
             # alpha = 0.5,
             size = 1.5,
             axes = F) +
  ggtitle(paste0(obj$TMA[1], " --- ", table(obj$Group) %>% names()))

table(obj$major_celltype, 
      obj$highlight_anno) # correct 




# preallocate list to hold all plots
plots_list <- list()

for (i in seq_along(Sobj_list)) {
  TMAs_1 <- unique(Sobj_list[[i]]$TMA)
  for (j in seq_along(TMAs_1)) {
    obj <- subset(Sobj_list[[i]], TMA %in% TMAs_1[j])
    Bcell_ID <- minor_lym@meta.data[
      minor_lym@meta.data$lym_nuo_granu_recluster == "Bcell" & 
        minor_lym@meta.data$TMA == TMAs_1[j], 
    ]$cell_id
    
    obj@meta.data$highlight_anno <- "other"
    obj@meta.data[obj@meta.data$cell_id %in% Bcell_ID, ]$highlight_anno <- "Bcell"
    
    p <- ImageDimPlot(
      obj,
      group.by = 'highlight_anno',
      border.color = "white",
      fov = "FOV",
      border.size = 0.1,
      cols = c("gray", "gold"),
      size = 1.5,
      axes = FALSE
    ) + ggtitle(paste0(obj$TMA[1], " --- ", names(table(obj$Group))))
    
    # Save plot to list with a unique name key
    plots_list[[paste0("i", i, "_j", j)]] <- p
  }
}

# The plots_list now contains all plots

plots_list$i1_j2



# Open PDF device — all plots go here, one per page
pdf("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Plots/looking_for_bcell_sp_structure.pdf", width = 8, height = 6) 

# Loop through all plots in the list and print each on a new page
for (plot_name in names(plots_list)) {
  print(plots_list[[plot_name]])
}

# Close the PDF device
dev.off()



## NEXT:  ------------------
# WE CAN do Felipe's code for cell proximity using only Bcells  



## BEFORE THAT, let's fix this TMAs .... -----------------
# from the looking_for_bcell_sp_structure.pdf I can locate which TMAs have more than one 
Sobj_list$Adult_1@meta.data[Sobj_list$Adult_1@meta.data$Group %in% "Adult Control", ]$fov %>% hist()

Sobj_list$Adult_2@meta.data[Sobj_list$Adult_2@meta.data$TMA %in% "Ap17−0332 A1 2" & Sobj_list$Adult_2@meta.data$Group %in% "Adult Control", ]$fov %>% hist()

Sobj_list$Adult_2@meta.data[Sobj_list$Adult_2@meta.data$TMA %in% "Ap17-0332 A1 2", ]$fov %>% hist()

obj <- subset(Sobj_list[["Adult_2"]], TMA %in% "Ap17-0332 A1 2")
ImageFeaturePlot(obj, 
                 features = c("nCount_RNA", "nFeature_RNA")) &
  scale_fill_viridis() &
  #coord_flip() &
  ggtitle(obj@meta.data$Group[1]) &
  theme(plot.title = element_text(size = 14, hjust = 0))



# Adult_2 - Hannah's ppt slides - 16 TMAs 
table(Sobj_list[["Adult_2"]]$TMA) %>% names() # 14 TMAs 
table(Sobj_list[["Adult_2"]]$TMA_2) %>% names() #I made TMA_2 to be able to match Diego's group table

Sobj_list[["Adult_2"]]@meta.data[Sobj_list[["Adult_2"]]@meta.data$TMA %in% "", ] $fov

# Checking the groups we did DE analysis with
comparisons <- list(
  "Comparison" = list(
    group_1 = "Neonate RSV infected (NO IFN)",
    group_2 = "Neonate IFN and RSV infected",
    group_3 = "Adult RSV infected"
  )
)

## 1. Correct at the TMA level ------
# looking at 
# Diego GeoMx and CosMx Sample Namining and ROIs.csv
# AAAAND 
# Full naming both projects.pptx
# just to make sure there's no type from ppt to spreadsheet
# files at /Users/mmarcao/Library/CloudStorage/OneDrive-St.JudeChildren'sResearchHospital/Jasmine_group/Diego_RSV/CosMx/common_files


Sobj_list_TMA_review <- Sobj_list

## Adult_1
# Ap17-0333
Sobj_list_TMA_review[["Adult_1"]]@meta.data$TMA_3 <- "NA"
Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 1:9, ]$TMA_3 <- "Ap17-0333 A1 1"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 10:18, ]$TMA_3 <- "Ap17-0333 A1 2"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 19:28, ]$TMA_3 <- "Ap17-0333 A1 3"


# Ap17-0332
Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 29:37, ]$TMA_3 <- "Ap17-0332 A1 1"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 38:46, ]$TMA_3 <- "Ap17-0332 A1 2"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 48:57, ]$TMA_3 <- "Ap17-0332 A1 3"


# Ap17-0331
Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 58:67, ]$TMA_3 <- "Ap17-0331 A1 1"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 68:76, ]$TMA_3 <- "Ap17-0331 A1 2"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 77:85, ]$TMA_3 <- "Ap17-0331 A1 3"


# Ap16-3781
Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 86:91, ]$TMA_3 <- "Ap16-3781 A1 1"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 111:120, ]$TMA_3 <- "Ap16-3781 A1 2"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 92:100, ]$TMA_3 <- "Ap16-3781 A1 3"


# Ap16-3780
Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 101:110, ]$TMA_3 <- "Ap16-3780 A1 1"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 121:129, ]$TMA_3 <- "Ap16-3780 A1 2"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 130:136, ]$TMA_3 <- "Ap16-3780 A1 3"


# Ap16-3779
Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 156:164, ]$TMA_3 <- "Ap16-3779 A1 1"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 146:155, ]$TMA_3 <- "Ap16-3779 A1 2"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 137:145, ]$TMA_3 <- "Ap16-3779 A1 3"


# Ap16-3896
Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 233:246, ]$TMA_3 <- "Ap16-3896 A1 1"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 201:212, ]$TMA_3 <- "Ap16-3896 A1 2"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 213:219, ]$TMA_3 <- "Ap16-3896 A1 3"


# Ap16-3894
Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 181:187, ]$TMA_3 <- "Ap16-3894 A1 1"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 188:196, ]$TMA_3 <- "Ap16-3894 A1 2"

Sobj_list_TMA_review[["Adult_1"]]@meta.data[
  Sobj_list_TMA_review[["Adult_1"]]@meta.data$fov %in% 218:232, ]$TMA_3 <- "Ap16-3894 A1 3"

table(Sobj_list_TMA_review[["Adult_1"]]@meta.data$TMA_3)

## Adult_2
Sobj_list_TMA_review[["Adult_2"]]@meta.data$TMA_3 <- "NA"

# Ap16-3781 A1
Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 1:10, ]$TMA_3 <- "Ap16-3781 A1 1"

Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 29:37, ]$TMA_3 <- "Ap16-3781 A1 2"


# Ap16-3780 A1
Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 38:49, ]$TMA_3 <- "Ap16-3780 A1 1"

Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 11:19, ]$TMA_3 <- "Ap16-3780 A1 2"


# Ap16-3894 A1
Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 59:74, ]$TMA_3 <- "Ap16-3894 A1 1"

Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 100:113, ]$TMA_3 <- "Ap16-3894 A1 2"


# Ap16-3779 A1
Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 20:28, ]$TMA_3 <- "Ap16-3779 A1 1"

Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 50:58, ]$TMA_3 <- "Ap16-3779 A1 2"


# Ap17-0331 A1
Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 75:83, ]$TMA_3 <- "Ap17-0331 A1 1"

Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 147:161, ]$TMA_3 <- "Ap17-0331 A1 2"


# Ap16-3896 A1
Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 114:128, ]$TMA_3 <- "Ap16-3896 A1 1"

Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 84:99, ]$TMA_3 <- "Ap16-3896 A1 2"


# Ap17-0333 A1
Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 129:137, ]$TMA_3 <- "Ap17-0333 A1 1"

Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 162:173, ]$TMA_3 <- "Ap17-0333 A1 2"


# Ap17-0332 A1
Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% 175:188, ]$TMA_3 <- "Ap17-0332 A1 1"

Sobj_list_TMA_review[["Adult_2"]]@meta.data[
  Sobj_list_TMA_review[["Adult_2"]]@meta.data$fov %in% c(138:146, 174), ]$TMA_3 <- "Ap17-0332 A1 2"

table(Sobj_list_TMA_review[["Adult_2"]]@meta.data$TMA_3)


## Neonate_1
Sobj_list_TMA_review[["Neonate_1"]]@meta.data$TMA_3 <- "NA"

# Ap16-2658
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 170:178, ]$TMA_3 <- "Ap16-2658 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 218:226, ]$TMA_3 <- "Ap16-2658 A1 2"


# Ap16-2656
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 179:187, ]$TMA_3 <- "Ap16-2656 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 149:160, ]$TMA_3 <- "Ap16-2656 A1 2"


# Ap16-2657
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 161:169, ]$TMA_3 <- "Ap16-2657 A1 1"


# Ap16-2655
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% c(197:199, 204:205, 211:213), ]$TMA_3 <- "Ap16-2655 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% c(200:203, 207:210, 214:217), ]$TMA_3 <- "Ap16-2655 A1 2"


# Ap16-2653
# N/A block listed; skipped, add range if needed
# Sobj_list_TMA_review[["Neonate_2"]]@meta.data[...]$TMA_3 <- "Ap16-2653 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 188:196, ]$TMA_3 <- "Ap16-2653 A1 2"


# Ap16-2654
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 135:143, ]$TMA_3 <- "Ap16-2654 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 144:148, ]$TMA_3 <- "Ap16-2654 A1 2"


# Ap16-0326
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 1:9, ]$TMA_3 <- "Ap16-0326 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 10:18, ]$TMA_3 <- "Ap16-0326 A1 2"


# Ap16-0327
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 19:28, ]$TMA_3 <- "Ap16-0327 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 29:37, ]$TMA_3 <- "Ap16-0327 A1 2"


# Ap16-3768
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% c(38:46, 51), ]$TMA_3 <- "Ap16-3768 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 47:54, ]$TMA_3 <- "Ap16-3768 A1 2"


# Ap16-0325
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 55:63, ]$TMA_3 <- "Ap16-0325 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 64:70, ]$TMA_3 <- "Ap16-0325 A1 2"


# Ap16-3770
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 71:79, ]$TMA_3 <- "Ap16-3770 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 80:88, ]$TMA_3 <- "Ap16-3770 A1 2"


# Ap16-3769
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 89:97, ]$TMA_3 <- "Ap16-3769 A1 1"

# N/A block listed; skipped, add range if available
# Sobj_list_TMA_review[["Neonate_1"]]@meta.data[...]$TMA_3 <- "Ap16-3769 A1 2"


# Ap16-3772
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 98:106, ]$TMA_3 <- "Ap16-3772 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 107:116, ]$TMA_3 <- "Ap16-3772 A1 2"


# Ap16-3771
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 117:122, ]$TMA_3 <- "Ap16-3771 A1 1"

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_1"]]@meta.data$fov %in% 123:134, ]$TMA_3 <- "Ap16-3771 A1 2"

table(Sobj_list_TMA_review[["Neonate_1"]]$TMA_3)


## Neonate_2
Sobj_list_TMA_review[["Neonate_2"]]@meta.data$TMA_3 <- "NA"

# Ap16-2653 A2
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 20:28, ]$TMA_3 <- "Ap16-2653 A2 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 11:19, ]$TMA_3 <- "Ap16-2653 A2 2"


# Ap16-2654 A2
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 1:10, ]$TMA_3 <- "Ap16-2654 A2 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 56:64, ]$TMA_3 <- "Ap16-2654 A2 2"


# Ap16-2655 A2
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 47:55, ]$TMA_3 <- "Ap16-2655 A2 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 38:46, ]$TMA_3 <- "Ap16-2655 A2 2"


# Ap16-2656 A2
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 29:37, ]$TMA_3 <- "Ap16-2656 A2 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 92:100, ]$TMA_3 <- "Ap16-2656 A2 2"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 83:91, ]$TMA_3 <- "Ap16-2656 A2 3"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 220:228, ]$TMA_3 <- "Ap16-2656 A2 4"


# Ap16-2657 A2
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 74:82, ]$TMA_3 <- "Ap16-2657 A2 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 65:73, ]$TMA_3 <- "Ap16-2657 A2 2"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 211:219, ]$TMA_3 <- "Ap16-2657 A2 3"


# Ap16-26587 A2
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 119:127, ]$TMA_3 <- "Ap16-26587 A2 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 110:118, ]$TMA_3 <- "Ap16-26587 A2 2"


# Ap17-0325 A1
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 240:253, ]$TMA_3 <- "Ap17-0325 A1 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 101:109, ]$TMA_3 <- "Ap17-0325 A1 2"


# Ap17-0326 A1
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% c(137:145, 201), ]$TMA_3 <- "Ap17-0326 A1 1"


# Ap17-0327 A1
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% c(128:136, 239), ]$TMA_3 <- "Ap17-0327 A1 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 266:277, ]$TMA_3 <- "Ap17-0327 A1 2"


# Ap16-3768 A1
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 254:265, ]$TMA_3 <- "Ap16-3768 A1 1"

# N/A block listed here
# Sobj_list_TMA_review[["Neonate_2"]]@meta.data[...]$TMA_3 <- "Ap16-3768 A1 2"


# Ap16-3769 A1
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 278:291, ]$TMA_3 <- "Ap16-3769 A1 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 155:163, ]$TMA_3 <- "Ap16-3769 A1 2"


# Ap16-3770 A1
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 146:154, ]$TMA_3 <- "Ap16-3770 A1 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 191:200, ]$TMA_3 <- "Ap16-3770 A1 2"


# Ap16-3771 A1
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 182:190, ]$TMA_3 <- "Ap16-3771 A1 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 173:181, ]$TMA_3 <- "Ap16-3771 A1 2"


# Ap16-3772 A1
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 164:172, ]$TMA_3 <- "Ap16-3772 A1 1"

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 229:238, ]$TMA_3 <- "Ap16-3772 A1 2"


# Ap16-36587 A2  (possibly a typo of 26587?)
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[
  Sobj_list_TMA_review[["Neonate_2"]]@meta.data$fov %in% 202:210, ]$TMA_3 <- "Ap16-36587 A2 1"


table(Sobj_list_TMA_review[["Neonate_2"]]@meta.data$TMA_3)

## 2. Match it again ---------


query <- table(Sobj_list_TMA_review[["Neonate_2"]]@meta.data$TMA) %>% names()
reference <- table(Sobj_list_TMA_review[["Neonate_2"]]@meta.data$TMA_3) %>% names()

length(intersect(query,
                 reference)) # right 

# Groups in query but not in reference
setdiff(query, reference)

# Groups in reference but not in query
setdiff(reference, query)

# Checking to see if there's any mismatch in the groups I've done DE analysis 
comparisons <- list(
  "Comparison" = list(
    group_1 = "Neonate RSV infected (NO IFN)",
    group_2 = "Neonate IFN and RSV infected",
    group_3 = "Adult RSV infected"
  )
)

# "Neonate RSV infected (NO IFN)" - ok
Sobj_list_TMA_review[["Adult_1"]]@meta.data[Sobj_list_TMA_review[["Adult_1"]]@meta.data$Group %in% comparisons$Comparison$group_1, ]$TMA %>% unique()

Sobj_list_TMA_review[["Adult_1"]]@meta.data[Sobj_list_TMA_review[["Adult_1"]]@meta.data$Group %in% comparisons$Comparison$group_1, ]$TMA_3 %>% unique()


# "Neonate IFN and RSV infected" - ok
Sobj_list_TMA_review[["Adult_1"]]@meta.data[Sobj_list_TMA_review[["Adult_1"]]@meta.data$Group %in% comparisons$Comparison$group_2, ]$TMA %>% unique()

Sobj_list_TMA_review[["Adult_1"]]@meta.data[Sobj_list_TMA_review[["Adult_1"]]@meta.data$Group %in% comparisons$Comparison$group_2, ]$TMA_3 %>% unique()


# "Adult RSV infected" - missing "Ap16-3780 A1 3"
Sobj_list_TMA_review[["Adult_1"]]@meta.data[Sobj_list_TMA_review[["Adult_1"]]@meta.data$Group %in% comparisons$Comparison$group_3, ]$TMA %>% unique()

Sobj_list_TMA_review[["Adult_1"]]@meta.data[Sobj_list_TMA_review[["Adult_1"]]@meta.data$Group %in% comparisons$Comparison$group_3, ]$TMA_3 %>% unique()




# "Neonate RSV infected (NO IFN)" - ok
Sobj_list_TMA_review[["Adult_2"]]@meta.data[Sobj_list_TMA_review[["Adult_2"]]@meta.data$Group %in% comparisons$Comparison$group_1, ]$TMA %>% unique()

Sobj_list_TMA_review[["Adult_2"]]@meta.data[Sobj_list_TMA_review[["Adult_2"]]@meta.data$Group %in% comparisons$Comparison$group_1, ]$TMA_3 %>% unique()


# "Neonate IFN and RSV infected" - ok
Sobj_list_TMA_review[["Adult_2"]]@meta.data[Sobj_list_TMA_review[["Adult_2"]]@meta.data$Group %in% comparisons$Comparison$group_2, ]$TMA %>% unique()

Sobj_list_TMA_review[["Adult_2"]]@meta.data[Sobj_list_TMA_review[["Adult_2"]]@meta.data$Group %in% comparisons$Comparison$group_2, ]$TMA_3 %>% unique()


# "Adult RSV infected" - ok
Sobj_list_TMA_review[["Adult_2"]]@meta.data[Sobj_list_TMA_review[["Adult_2"]]@meta.data$Group %in% comparisons$Comparison$group_3, ]$TMA %>% unique()

Sobj_list_TMA_review[["Adult_2"]]@meta.data[Sobj_list_TMA_review[["Adult_2"]]@meta.data$Group %in% comparisons$Comparison$group_3, ]$TMA_3 %>% unique()




# "Neonate RSV infected (NO IFN)" - 3/4 ok, missing 1 TMA 
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[Sobj_list_TMA_review[["Neonate_1"]]@meta.data$Group %in% comparisons$Comparison$group_1, ]$TMA %>% unique()

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[Sobj_list_TMA_review[["Neonate_1"]]@meta.data$Group %in% comparisons$Comparison$group_1, ]$TMA_3 %>% unique()


# "Neonate IFN and RSV infected" - ok
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[Sobj_list_TMA_review[["Neonate_1"]]@meta.data$Group %in% comparisons$Comparison$group_2, ]$TMA %>% unique()

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[Sobj_list_TMA_review[["Neonate_1"]]@meta.data$Group %in% comparisons$Comparison$group_2, ]$TMA_3 %>% unique()


# "Adult RSV infected" - ok
Sobj_list_TMA_review[["Neonate_1"]]@meta.data[Sobj_list_TMA_review[["Neonate_1"]]@meta.data$Group %in% comparisons$Comparison$group_3, ]$TMA %>% unique()

Sobj_list_TMA_review[["Neonate_1"]]@meta.data[Sobj_list_TMA_review[["Neonate_1"]]@meta.data$Group %in% comparisons$Comparison$group_3, ]$TMA_3 %>% unique()



# "Neonate RSV infected (NO IFN)" - ok 
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[Sobj_list_TMA_review[["Neonate_2"]]@meta.data$Group %in% comparisons$Comparison$group_1, ]$TMA %>% unique()

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[Sobj_list_TMA_review[["Neonate_2"]]@meta.data$Group %in% comparisons$Comparison$group_1, ]$TMA_3 %>% unique()


# "Neonate IFN and RSV infected" - 1/1 mismatch ...
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[Sobj_list_TMA_review[["Neonate_2"]]@meta.data$Group %in% comparisons$Comparison$group_2, ]$TMA %>% unique()

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[Sobj_list_TMA_review[["Neonate_2"]]@meta.data$Group %in% comparisons$Comparison$group_2, ]$TMA_3 %>% unique()


# "Adult RSV infected" - ok 
Sobj_list_TMA_review[["Neonate_2"]]@meta.data[Sobj_list_TMA_review[["Neonate_2"]]@meta.data$Group %in% comparisons$Comparison$group_3, ]$TMA %>% unique()

Sobj_list_TMA_review[["Neonate_2"]]@meta.data[Sobj_list_TMA_review[["Neonate_2"]]@meta.data$Group %in% comparisons$Comparison$group_3, ]$TMA_3 %>% unique()



# Plotting Bcells again - 
# Better to go back and call cell type  for each TMA individually 










