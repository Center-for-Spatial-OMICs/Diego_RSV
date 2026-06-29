library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(readr)

library(GeomxTools)
library(NanoStringNCTools)


# --- Load GeoMx exp and metadata ---

## Other paths
## /mnt/plummergrp/Temp/Diego-RSV # Yutian's path
## /mnt/scratch1/maycon/Diego_RSV_CosMx/Round_2 # this work env

## Load GeoMx expression mtx
geomx_exp_mtx <- read_csv("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_2/Objects/ex_df_diego_rsv.csv") %>% data.frame()
head(geomx_exp_mtx)
dim(geomx_exp_mtx) #19963    93
head(geomx_exp_mtx)
rownames(geomx_exp_mtx) <- geomx_exp_mtx$...1
geomx_exp_mtx$...1 <- NULL
head(geomx_exp_mtx)

## Load GeoMx target/Nanostring obj
target_demoData <- readRDS("/mnt/plummergrp/Temp/Diego-RSV/target_demoData.rds")

View(target_demoData)

target_demoData@assayData # exp mtx and it's normalized versions
target_demoData@phenoData@data # metadata, coordinates etc
head(target_demoData@phenoData@data)
metadata <- target_demoData@phenoData@data
metadata$`Slide Name` %>% table()
meta_adult_1 <- metadata[metadata$`Slide Name` %in% "Adult.1_8.30.24", ]
dim(meta_adult_1)
table(meta_adult_1$Sample_ID)

head(meta_adult_1)

## ROI ID
obj <- protocolData(target_demoData)
obj@data %>% head()
meta_ROI <- obj@data
head(meta_ROI)
meta_ROI_adult_1 <- meta_ROI[rownames(meta_adult_1), ]
table(meta_ROI_adult_1$Roi) # YEEEES

current_cat <- meta_ROI_adult_1$Roi
new_categories <- seq(1:33) %>% as.character()

meta_ROI_adult_1$ROI_ID <- plyr::mapvalues(meta_ROI_adult_1$Roi,
                                           from = current_cat, 
                                           to = new_categories)

meta_ROI_adult_1$Sample_ID_2 <- meta_ROI_adult_1 %>% rownames()
meta_ROI_adult_1 <- meta_ROI_adult_1[, c("ROI_ID", "Sample_ID_2")]
head(meta_adult_1)
meta_adult_1$Sample_ID_2 <- meta_adult_1 %>% rownames()

meta_adult_1 <- merge(meta_adult_1, meta_ROI_adult_1, by="Sample_ID_2")
meta_adult_1$Sample_ID_2 <- gsub("-", ".",meta_adult_1$Sample_ID_2)

## Filter GeoMx objects
## Exp mtx - this is what Yutian had used for Differental Expression Analysis  
## **We could use other transformations as q_norm (already in target_demoData@assayData)
Exp_mtx_adult_1_GeoMx <- geomx_exp_mtx[, colnames(geomx_exp_mtx) %in% meta_adult_1$Sample_ID_2]

GeoMx_adult_1 <- list(exp_mtx = Exp_mtx_adult_1_GeoMx,
                      metadata = meta_adult_1)

GeoMx_adult_1$exp_mtx[1:4, 1:4]

# --- Load CosMx exp and metadata ---
Adult12_Neonate12_seuobj_list <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult12_Neonate12_seuobj_list.rds")

Sobj_adult_1_CosMx <- Adult12_Neonate12_seuobj_list$Adult_1


# --- Filter only 1 FOV and it's ROIs --- 
## FOV-ROI relation from Hanna's spreadsheet at /mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Diego GeoMx and CosMx Sample Namining and ROIs.xlsx

Ap16_3894_FOV_218_232 <- subset(Sobj_adult_1_CosMx, fov %in% c(as.character(218:232)))

Adult_1_ROI_20_24_metadata <- GeoMx_adult_1$metadata[GeoMx_adult_1$metadata$ROI_ID %in%  c(as.character(20:24)), ]
Adult_1_ROI_20_24_exp_mtx <- select(GeoMx_adult_1$exp_mtx, Adult_1_ROI_20_24_metadata$Sample_ID_2)

Adult_1_ROI_ROI_20_24 <- list(Adult_1_ROI_20_24_metadata = Adult_1_ROI_20_24_metadata,
                      Adult_1_ROI_20_24_exp_mtx = Adult_1_ROI_20_24_exp_mtx)


Adult_2_ROI_7$Adult_2_ROI_7_exp_mtx$DSP.1001660013798.H.F04.dcc %>% hist(breaks = 100000)



Ap16_3894_FOV_218_232




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













