library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(readr)

library(GeomxTools)
library(NanoStringNCTools)


## Other paths
## /mnt/plummergrp/Temp/Diego-RSV # Yutian's path
## /mnt/scratch1/maycon/Diego_RSV_CosMx/Round_2 # this work env

## Load GeoMx expression mtx -------------
geomx_exp_mtx <- read_csv("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_2/Objects/ex_df_diego_rsv.csv") %>% data.frame()
head(geomx_exp_mtx)
dim(geomx_exp_mtx) #19963    93
head(geomx_exp_mtx)
rownames(geomx_exp_mtx) <- geomx_exp_mtx$...1
geomx_exp_mtx$...1 <- NULL
head(geomx_exp_mtx)

## Load GeoMx target/Nanostring obj ---------
target_demoData <- readRDS("/mnt/plummergrp/Temp/Diego-RSV/target_demoData.rds")

View(target_demoData)

target_demoData@assayData # exp mtx and it's normalized versions
target_demoData@phenoData@data # metadata, coordinates etc
head(target_demoData@phenoData@data)
metadata <- target_demoData@phenoData@data
metadata$`Slide Name` %>% table()
meta_adult_2 <- metadata[metadata$`Slide Name` %in% "Adult.2_8.30.24", ]
dim(meta_adult_2)
table(meta_adult_2$Sample_ID)

head(meta_adult_2)

## ROI ID
obj <- protocolData(target_demoData)
obj@data %>% head()
meta_ROI <- obj@data
head(meta_ROI)
meta_ROI_adult_2 <- meta_ROI[rownames(meta_adult_2), ]
table(meta_ROI_adult_2$Roi) # YEEEES

current_cat <- meta_ROI_adult_2$Roi
new_categories <- seq(1:16) %>% as.character()
  
current_categories <- names(table(obj_1_5@meta.data$spatial_cluster))
new_categories <- setNames(paste0("TMA_", seq_along(current_categories)), current_categories)

meta_ROI_adult_2$ROI_ID <- plyr::mapvalues(meta_ROI_adult_2$Roi,
                                                     from = current_cat, 
                                                     to = new_categories)

meta_ROI_adult_2$Sample_ID_2 <- meta_ROI_adult_2 %>% rownames()
meta_ROI_adult_2 <- meta_ROI_adult_2[, c("ROI_ID", "Sample_ID_2")]
head(meta_adult_2)
meta_adult_2$Sample_ID_2 <- meta_adult_2 %>% rownames()

meta_adult_2 <- merge(meta_adult_2, meta_ROI_adult_2, by="Sample_ID_2")
meta_adult_2$Sample_ID_2 <- gsub("-", ".",meta_adult_2$Sample_ID_2)

## Filter GeoMx objects
## Exp mtx - this is what Yutian had used for Differental Expression Analysis  
## **We could use other transformations as q_norm (already in target_demoData@assayData)
Exp_mtx_adult_2_GeoMx <- geomx_exp_mtx[, colnames(geomx_exp_mtx) %in% meta_adult_2$Sample_ID_2]

GeoMx_adult_2 <- list(exp_mtx = Exp_mtx_adult_2_GeoMx,
                      metadata = meta_adult_2)



## Load CosMx objects -------------
Sobj_adult_2_CosMx <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult_2_obj.rds")



## Analysis -----------
# ROI is from GeoMx 
# FOV is from CosMx


## Adult_2 CosMx FOV 
# fov_list <- list(
#   'Ap16-3781 A1 1' = 1:10,
#   'Ap16-3780 A1 1' = 38:49,
#   'Ap16-3780 A1 2' = 11:19,
#   'Ap16-3894 A1 1' = 59:74,
#   'Ap16-3779 A1 1' = 20:28,
#   'Ap16-3779 A1 2' = 50:58,
#   'Ap16-3781 A1 2' = 29:37,
#   'Ap17-0331 A1 1' = 75:83,
#   'Ap16-3896 A1 1' = 114:128,
#   'Ap16-3896 A1 2' = 84:99,
#   'Ap16-3894 A1 2' = 100:113,
#   'Ap17-0333 A1 1' = 129:137,
#   'Ap17-0332 A1 1' = 175:188,
#   'Ap17-0332 A1 2' = 138:174,
#   'Ap17-0331 A1 2' = 147:161,
#   'Ap17-0333 A1 2' = 162:173
# )

fov_roi_df <- data.frame(
  FOV_start = c(1, 38, 11, 59, 20, 50, 29, 75, 114, 84, 100, 129, 175, 138, 147, 162),
  FOV_end = c(10, 49, 19, 74, 28, 58, 37, 83, 128, 99, 113, 137, 188, 174, 161, 173),
  ROI = c("3", "2", "1", "7", "6", "5", "4", "11", "10", "9", "8", "15", "14", "13", "12", "16"),
  Identifier = c('Ap16-3781 A1 1', 'Ap16-3780 A1 1', 'Ap16-3780 A1 2', 'Ap16-3894 A1 1',
                 'Ap16-3779 A1 1', 'Ap16-3779 A1 2', 'Ap16-3781 A1 2', 'Ap17-0331 A1 1',
                 'Ap16-3896 A1 1', 'Ap16-3896 A1 2', 'Ap16-3894 A1 2', 'Ap17-0333 A1 1',
                 'Ap17-0332 A1 1', 'Ap17-0332 A1 2', 'Ap17-0331 A1 2', 'Ap17-0333 A1 2')
)

## First ROI  

# Objects 
Sobj_adult_2_CosMx <- Sobj_adult_2_CosMx %>%
  NormalizeData() 
GeoMx_adult_2
fov_roi_df

geomx_meta <- GeoMx_adult_2$metadata[GeoMx_adult_2$metadata$ROI_ID %in% 1, ]
geomx_exp <- GeoMx_adult_2$exp_mtx[, colnames(GeoMx_adult_2$exp_mtx) %in% geomx_meta$Sample_ID_2, drop = FALSE]

geomx_exp$DSP.1001660013798.H.E10.dcc %>% hist(breaks = 100)
range(geomx_exp$DSP.1001660013798.H.E10.dcc)




library(tidyverse)
library(Seurat)

# Assuming 'geomx_data' is your expression matrix from GeoMx
# and 'cosmx_seurat' is your Seurat object for CosMx data


GeoMx_adult_2 #non-normalized counts from Yutian .csv

cosmx_seurat <- Sobj_adult_2_CosMx #seurat normalized data 

fov_roi_df


## From GeoMx to CosMx 
check_expression <- function(rois, GeoMx_adult_2, cosmx_seurat, fov_roi_df) {
  result_list <- list()
  
  for (roi in rois) {
    # Filtering GeoMx by ROI-FOV
    geomx_meta <- GeoMx_adult_2$metadata[GeoMx_adult_2$metadata$ROI_ID %in% roi, ]
    geomx_exp <- GeoMx_adult_2$exp_mtx[, colnames(GeoMx_adult_2$exp_mtx) %in% geomx_meta$Sample_ID_2, drop = FALSE]
    
    geomx_data <- geomx_exp # filtered to only one ROI
    
    top_genes_geomx <- geomx_data %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      gather(key = "sample", value = "expression", -gene) %>%
      group_by(gene) %>%
      summarize(mean_expression = mean(expression)) %>%
      arrange(desc(mean_expression)) %>%
      head(100) %>%
      pull(gene)
    
    # Filtering CosMx by ROI-FOV
    fov_range <- fov_roi_df %>%
      filter(ROI == roi) %>%
      select(FOV_start, FOV_end)
    
    cosmx_cells <- subset(cosmx_seurat, fov %in% (fov_range$FOV_start:fov_range$FOV_end)) %>% colnames()
    
    # Combine results for both platforms
    for (gene in top_genes_geomx) {
     
      if (gene %in% rownames(cosmx_seurat)) {
        cosmx_gene_expr <- GetAssayData(cosmx_seurat, slot = "data")[gene, cosmx_cells]
        cosmx_mean_expr <- mean(cosmx_gene_expr, na.rm = TRUE)
        cosmx_percent_expressed <- mean(cosmx_gene_expr > 0, na.rm = TRUE) * 100
      } else {
        cosmx_mean_expr <- NA
        cosmx_percent_expressed_percell <- NA
      }
      
      result_list[[length(result_list) + 1]] <- data.frame(
        ROI = roi,
        FOV_start = fov_range$FOV_start,
        FOV_end = fov_range$FOV_end,
        Gene = gene,
        CosMx_Mean_Expr = cosmx_mean_expr,
        CosMx_Percent_Expressed = cosmx_percent_expressed
      )
    }
  }
  
  return(do.call(rbind, result_list))
}


results <- check_expression(rois = fov_roi_df$ROI, 
                              GeoMx_adult_2 = GeoMx_adult_2,
                              cosmx_seurat = cosmx_seurat,
                              fov_roi_df = fov_roi_df) %>% filter(!is.na(CosMx_Mean_Expr))




# Analyze and visualize results
ggplot(results, aes(x = CosMx_Mean_Expr, y = CosMx_Percent_Expressed, color = ROI)) +
  geom_point() +
  facet_wrap(~Gene, scales = "free") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Top GeoMx Genes Expression in CosMx Cells",
       x = "Mean Expression in CosMx",
       y = "Percent of CosMx Cells Expressing")









## Kind of analysis I could perform 

## The problem: "where I have one .dcc file expression from GeoMx for multiple cells expressions from CosMx, how can I verify that indeed the most expressed genes we see on that single .dcc file is also detected in those cells on CosMx??"

# i) doing cell type deconvolution to have cell types for this ROI then we could compare to CosMx 
# Requirements: using all ROIs from adult_2; having a reference sc data with the cell types we want; I still need to call cell types for CosMx too

# ii) GeoMx DE within ROIs and CosMx FindAllMarkers within the corresponded FOVs (eg FOV 1:10 would be called "ROI_1"). Then we could compare the DEGs from both technologies 


# iii) Look up the highest expressed genes








# Questions 

## Why on ROI have only one spatial coordinate ?
GeoMx_adult_2$metadata[GeoMx_adult_2$metadata$ROI_ID %in% 1, ]$`ROI Coordinate X` 





















