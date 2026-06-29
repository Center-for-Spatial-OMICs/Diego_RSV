## Description: Creating Seurat Object for all slides (adult 1, neonate 1, adult 2, neoante 2)

## Load packages -------
library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)


## Define variables -------

Adult_1_dir <- "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult1831439_08_11_2024_9_49_46_243/flatFiles/831439_IND_CR_FFPE_TMA_Adult1_STJ_N_R1"

Neonate_1_dir <- "/mnt/scratch2/CosMX/Diego/DiegoRSVNeonate1831439_08_11_2024_11_18_07_116/flatFiles/831439_IND_CR_FFPE_TMA_Neonate1_STJ_N_R1"

Adult_2_dir <- "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult2831439_08_11_2024_12_58_14_961/flatFiles/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1"

Neonate_2_dir <- "/mnt/scratch2/CosMX/Diego/DIegoRSVNeonate2831439_08_11_2024_13_26_50_655/flatFiles/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1"



### Create Adult_1 Seurat Object ------------
Adult_1_obj <- LoadNanostring(Adult_1_dir, fov = 'FOV')

Adult_1_tx <- fread(list.files(path = Adult_1_dir, pattern = "tx.*\\.csv\\.gz$", full.names = TRUE)) %>% as.data.frame

Adult_1_meta <- fread(list.files(path = Adult_1_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)) %>% as.data.frame

Adult_1_obj$fov <- Adult_1_meta$fov
Adult_1_obj$Area <- Adult_1_meta$Area
Adult_1_obj$Area.um2 <- Adult_1_meta$Area.um2

# match cell id in meta obj
Adult_1_meta$cell = paste0(as.character(Adult_1_meta$cell_ID), "_", Adult_1_meta$fov)

# remove control and system probes
which <- grep('Neg*|Syst*', rownames(Adult_1_obj))
Adult_1_obj <- Adult_1_obj[-c(which), ]

# Remove low quality cells
Adult_1_obj <- subset(Adult_1_obj, nCount_Nanostring > 10 & nFeature_Nanostring > 10 & Area.um2 < 500)

#Update meta
Adult_1_meta <- Adult_1_meta[Adult_1_meta$cell %in% colnames(Adult_1_obj), ]
dim(Adult_1_meta) #104893     58



# Classify different embryos by FOV numbers
Adult_1_obj$TMA <- ""
Adult_1_obj[[]] <- Adult_1_obj[[]] %>%
  mutate(TMA = case_when(
    fov %in% 1:9 ~ 'Ap17-0333 1',
    fov %in% 10:18  ~ 'Ap17-0333 2',
    fov %in% 19:28 ~ 'Ap17-0333 3',
    fov %in%  29:37 ~ 'Ap17-0332 1',
    fov %in% 38:47 ~ 'Ap17-0332 2',
    fov %in%  48:57 ~ 'Ap17-0332 3',
    fov %in%  58:67 ~ 'Ap17-0331 1',
    fov %in%  68:76 ~ 'Ap17-0331 2',
    fov %in%   77:85 ~ 'Ap17-0331 3',
    fov %in%   86:91 ~ 'Ap16-3781 1',
    fov %in%   111:123 ~ 'Ap16-3781 2',
    fov %in%   92:100 ~ 'Ap16-3781 3',
    fov %in%   101:109 ~ 'Ap16-3780 1',
    fov %in%    121:129~ 'Ap16-3780 2',
    fov %in%    130:136~ 'Ap16-3780 2',
    fov %in%    156:164~ 'Ap16-3779 1',
    fov %in%    146:155~ 'Ap16-3779 2',
    fov %in%    137:145~ 'Ap16-3779 3',
    fov %in%    234:246~ 'Ap16-3896 1',
    fov %in%    201:212~ 'Ap16-3896 2',
    fov %in%    c(213:217) ~ 'Ap16-3896 3',
    fov %in%    181:187~ 'Ap16-3894 1',
    fov %in%    188:196~ 'Ap16-3894 2',
    fov %in%    218:232~ 'Ap16-3894 3',
    
    TRUE ~ TMA  # Retain existing value or specify a default if needed
  ))

table(Adult_1_obj$TMA, useNA = 'always')

#saveRDS(Adult_1_obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult_1_obj.rds")


### Create Neonate_1 Seurat Object ------------
Neonate_1_obj <- LoadNanostring(Neonate_1_dir, fov = 'FOV')

Neonate_1_tx <- fread(list.files(path = Neonate_1_dir, pattern = "tx.*\\.csv\\.gz$", full.names = TRUE)) %>% as.data.frame

Neonate_1_meta <- fread(list.files(path = Neonate_1_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)) %>% as.data.frame

Neonate_1_obj$fov <- Neonate_1_meta$fov
Neonate_1_obj$Area <- Neonate_1_meta$Area
Neonate_1_obj$Area.um2 <- Neonate_1_meta$Area.um2

# match cell id in meta obj
Neonate_1_meta$cell = paste0(as.character(Neonate_1_meta$cell_ID), "_", Neonate_1_meta$fov)

# remove control and system probes
which <- grep('Neg*|Syst*', rownames(Neonate_1_obj))
Neonate_1_obj <- Neonate_1_obj[-c(which), ]

# Remove low quality cells
Neonate_1_obj <- subset(Neonate_1_obj, nCount_Nanostring > 10 & nFeature_Nanostring > 10 & Area.um2 < 500)

#Update meta
Neonate_1_meta <- Neonate_1_meta[Neonate_1_meta$cell %in% colnames(Neonate_1_obj), ]
dim(Neonate_1_meta) #104893     58



# Classify different embryos by FOV numbers
Neonate_1_obj$TMA <- ""
Neonate_1_obj[[]] <- Neonate_1_obj[[]] %>%
  mutate(TMA = case_when(
    fov %in%  170:178~ 'Ap16-2658 1',
    fov %in%  218:226~ 'Ap16-2658 2',
    fov %in%  179:187~ 'Ap16-2656 1',
    fov %in%  161:169~ 'Ap16-2657',
    fov %in%  c(197:199, 204:206, 211:213) ~ 'Ap16-2655 1',
    fov %in%  c(200:203, 207:210, 214:217)~ 'Ap16-2655 2',
    fov %in%  149:160~ 'Ap16-2656 2',
    fov %in%  188:196~ 'Ap16-2653',
    fov %in%  135:143~ 'Ap16-2654',
    fov %in%  144:148~ 'Ap16-2654',
    fov %in%  1:9~ 'Ap16-0326 1',
    fov %in%  10:18~ 'Ap16-0326 2',
    fov %in%  19:27~ 'Ap16-0327',
    fov %in%  29:37~ 'Ap16-0327',
    fov %in%  38:46~ 'Ap16-3768 1',
    fov %in%  47:56~ 'Ap16-3768 2',
    fov %in%  55:63~ 'Ap16-0325 1',
    fov %in%  64:70~ 'Ap16-0325 2',
    fov %in%  71:79~ 'Ap16-3770 1',
    fov %in%  80:88~ 'Ap16-3770 2',
    fov %in%  89:97~ 'Ap16-3769',
    fov %in%  98:106~ 'Ap16-3772 1',
    fov %in%  107:116~ 'Ap16-3772 2',
    fov %in%  117:122~ 'Ap16-3771 1',
    fov %in%  123:134~ 'Ap16-3771 2',

    TRUE ~ TMA  # Retain existing value or specify a default if needed
  ))

table(Neonate_1_obj$TMA, useNA = 'always')

#saveRDS(Neonate_1_obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Neonate_1_obj.rds")

### Create Adult_2 Seurat Object ------------
Adult_2_obj <- LoadNanostring(Adult_2_dir, fov = 'FOV')

Adult_2_tx <- fread(list.files(path = Adult_2_dir, pattern = "tx.*\\.csv\\.gz$", full.names = TRUE)) %>% as.data.frame

Adult_2_meta <- fread(list.files(path = Adult_2_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)) %>% as.data.frame

Adult_2_obj$fov <- Adult_2_meta$fov
Adult_2_obj$Area <- Adult_2_meta$Area
Adult_2_obj$Area.um2 <- Adult_2_meta$Area.um2

# match cell id in meta obj
Adult_2_meta$cell = paste0(as.character(Adult_2_meta$cell_ID), "_", Adult_2_meta$fov)

# remove control and system probes
which <- grep('Neg*|Syst*', rownames(Adult_2_obj))
Adult_2_obj <- Adult_2_obj[-c(which), ]

# Remove low quality cells
Adult_2_obj <- subset(Adult_2_obj, nCount_Nanostring > 10 & nFeature_Nanostring > 10 & Area.um2 < 500)

#Update meta
Adult_2_meta <- Adult_2_meta[Adult_2_meta$cell %in% colnames(Adult_2_obj), ]
dim(Adult_2_meta) #104893     58



# Classify different embryos by FOV numbers
Adult_2_obj$TMA <- ""
Adult_2_obj[[]] <- Adult_2_obj[[]] %>%
  mutate(TMA = case_when(
    fov %in% 1:10 ~ 'Ap16-3781 A1 1',
    fov %in% 38:49 ~ 'Ap16-3780 A1 1',
    fov %in% 11:19 ~ 'Ap16-3780 A1 2',
    fov %in% 59:74 ~ 'Ap16-3894 A1 1',
    fov %in% 20:28 ~ 'Ap16-3779 A1 1',
    fov %in% 50:58 ~ 'Ap16-3779 A1 2',
    fov %in% 29:37 ~ 'Ap16-3781 A1 2',
    fov %in% 75:83 ~ 'Ap17-0331 A1 1',
    fov %in% 114:128 ~ 'Ap16-3896 A1 1',
    fov %in% 84:99 ~ 'Ap16-3896 A1 2',
    fov %in% 100:113 ~ 'Ap16-3894 A1 2',
    fov %in% 129:137 ~ 'Ap17-0333 A1 1',
    fov %in% 175:188 ~ 'Ap17-0332 A1 1',
    fov %in% 138:174 ~ 'Ap17-0332 A1 2',
    fov %in% 147:161 ~ 'Ap17-0331 A1 2',
    fov %in% 162:173 ~ 'Ap17-0333 A1 2',
    
    TRUE ~ TMA  # Retain existing value or specify a default if needed
  ))

table(Adult_2_obj$TMA, useNA = 'always')

#saveRDS(Adult_2_obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult_2_obj.rds")


### Create Neonate_2 Seurat Object ------------
Neonate_2_obj <- LoadNanostring(Neonate_2_dir, fov = 'FOV')

Neonate_2_tx <- fread(list.files(path = Neonate_2_dir, pattern = "tx.*\\.csv\\.gz$", full.names = TRUE)) %>% as.data.frame

Neonate_2_meta <- fread(list.files(path = Neonate_2_dir, pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)) %>% as.data.frame


Neonate_2_obj$fov <- Neonate_2_meta$fov
Neonate_2_obj$Area <- Neonate_2_meta$Area
Neonate_2_obj$Area.um2 <- Neonate_2_meta$Area.um2

# match cell id in meta obj
Neonate_2_meta$cell = paste0(as.character(Neonate_2_meta$cell_ID), "_", Neonate_2_meta$fov)

# remove control and system probes
which <- grep('Neg*|Syst*', rownames(Neonate_2_obj))
Neonate_2_obj <- Neonate_2_obj[-c(which), ]

# Remove low quality cells
Neonate_2_obj <- subset(Neonate_2_obj, nCount_Nanostring > 10 & nFeature_Nanostring > 10 & Area.um2 < 500)

#Update meta
Neonate_2_meta <- Neonate_2_meta[Neonate_2_meta$cell %in% colnames(Neonate_2_obj), ]
dim(Neonate_2_meta) #

# Neonate_2_obj_bk <- Neonate_2_obj


# # Classify different embryos by FOV numbers
Neonate_2_obj$TMA <- ""
Neonate_2_obj[[]] <- Neonate_2_obj[[]] %>%
  mutate(TMA = case_when(
    fov %in% 1:10 ~ 'Ap16-2653 A2 1',
    fov %in% 11:19 ~ 'Ap16-2653 A2 2',
    fov %in% 20:28 ~ 'Ap16-2653 A2 3',
    fov %in% 29:37 ~ 'Ap16-2656 A2 1',
    fov %in% 38:46 ~ 'Ap16-2655 A2 1',
    fov %in% 47:55 ~ 'Ap16-2655 A2 2',
    fov %in% 56:64 ~ 'Ap16-2654 A2 1',
    fov %in% 65:73 ~ 'Ap16-2657 A2 1',
    fov %in% 74:82 ~ 'Ap16-2657 A2 2',
    fov %in% 83:91 ~ 'Ap16-2656 A2 2',
    fov %in% 92:100 ~ 'Ap16-2656 A2 3',
    fov %in% 101:109 ~ 'Ap17-0325 A1 1',
    fov %in% 240:253 ~ 'Ap17-0325 A1 2',
    fov %in% 110:118 ~ 'Ap16-26587 A2 1',
    fov %in% 119:127 ~ 'Ap16-26587 A2 2',
    fov %in% 254:265 ~ 'Ap16-3768 A1 1',
    fov %in% 266:277 ~ 'Ap17-0327 A1 1',
    fov %in% 128:136 ~ 'Ap17-0327 A1 2',
    fov %in% 137:201 ~ 'Ap17-0326 A1 1',
    fov %in% 146:154 ~ 'Ap16-3770 A1 1',
    fov %in% 155:166 ~ 'Ap16-3769 A1 1',
    fov %in% 278:291 ~ 'Ap16-3769 A1 2',
    fov %in% 164:172 ~ 'Ap16-3772 A1 1',
    fov %in% 173:181 ~ 'Ap16-3771 A1 1',
    fov %in% 182:190 ~ 'Ap16-3771 A1 2',
    fov %in% 191:200 ~ 'Ap16-3768 A1 2',
    fov %in% 202:210 ~ 'Ap16-26587 A2 3',
    fov %in% 211:219 ~ 'Ap16-2657 A2 3',
    fov %in% 220:228 ~ 'Ap16-2656 A2 4',
    fov %in% 229:238 ~ 'Ap16-3770 A1 3',
    
    TRUE ~ TMA  # Retain existing value or specify a default if needed
  ))

table(Neonate_2_obj$TMA, useNA = 'always')

#saveRDS(Neonate_2_obj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Neonate_2_obj.rds")





## Adding group labels - I had this info. after weeks ...
library(readxl)
library(stringr)
library(Seurat)
library(plyr)

group_meta <- read_excel("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Hijano_samples_GeoMX_CosMx (1)[100].xlsx") %>% data.frame()

group_meta$Accession.Number <- gsub("P", "p", group_meta$Accession.Number)

## Neonate_1_obj
obj <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Neonate_1_obj.rds")

first_10_chars <- substr(obj@meta.data$TMA, 1, 10)
first_10_chars <- gsub(" ", "", first_10_chars)
obj@meta.data$TMA_2 <- first_10_chars
obj$cell_ID <- rownames(obj@meta.data)

meta_fseu <- obj@meta.data %>% select(TMA_2, TMA, cell_ID)
meta_flabsheet <- group_meta %>% select(Accession.Number, Group.label)
names(meta_flabsheet) <- c("TMA_2", "Group")

TMA_ID_fseu <- meta_fseu$TMA_2 %>% as.factor() %>% levels()
TMA_ID_fseu <- TMA_ID_fseu[-1]
TMA_ID_flabsheet <- meta_flabsheet$TMA_2 %>% as.factor() %>% levels()

intersect(TMA_ID_fseu, TMA_ID_flabsheet)
TMA_ID_toChange <- setdiff(TMA_ID_fseu, TMA_ID_flabsheet)
print(TMA_ID_toChange)


mapping <- c("Ap16-0325" = "Ap17-0325",
             "Ap16-0326" = "Ap17-0326",
             "Ap16-0327" = "Ap17-0327")

meta_fseu$TMA_2 <- mapvalues(meta_fseu$TMA_2,
                                        from = names(mapping),
                                        to = mapping,
                                        warn_missing = FALSE)

meta_fseu <- merge(meta_fseu, meta_flabsheet, by="TMA_2")
rownames(meta_fseu) <- meta_fseu$cell_ID

obj = AddMetaData(
object = obj,
metadata = meta_fseu)

Neonate_1_obj <- obj


## Neonate_2_obj
obj <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Neonate_2_obj.rds")

first_10_chars <- substr(obj@meta.data$TMA, 1, 10)
first_10_chars <- gsub(" ", "", first_10_chars)
obj@meta.data$TMA_2 <- first_10_chars
obj$cell_ID <- rownames(obj@meta.data)

meta_fseu <- obj@meta.data %>% select(TMA_2, TMA, cell_ID)
meta_flabsheet <- group_meta %>% select(Accession.Number, Group.label)
names(meta_flabsheet) <- c("TMA_2", "Group")

TMA_ID_fseu <- meta_fseu$TMA_2 %>% as.factor() %>% levels()
TMA_ID_fseu <- TMA_ID_fseu[-1]
TMA_ID_flabsheet <- meta_flabsheet$TMA_2 %>% as.factor() %>% levels()

intersect(TMA_ID_fseu, TMA_ID_flabsheet)
TMA_ID_toChange <- setdiff(TMA_ID_fseu, TMA_ID_flabsheet)
print(TMA_ID_toChange)

mapping <- c("Ap16-26587" = "Ap16-2658") #based on the intersected ones, it's the Ap17-2658 that is missing 

meta_fseu$TMA_2 <- mapvalues(meta_fseu$TMA_2,
                             from = names(mapping),
                             to = mapping,
                             warn_missing = FALSE)

meta_fseu <- merge(meta_fseu, meta_flabsheet, by="TMA_2")
rownames(meta_fseu) <- meta_fseu$cell_ID

obj = AddMetaData(
  object = obj,
  metadata = meta_fseu)


Neonate_2_obj <- obj



## Adult_1_obj
obj <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult_1_obj.rds")

first_10_chars <- substr(obj@meta.data$TMA, 1, 10)
first_10_chars <- gsub(" ", "", first_10_chars)
obj@meta.data$TMA_2 <- first_10_chars
obj$cell_ID <- rownames(obj@meta.data)

meta_fseu <- obj@meta.data %>% select(TMA_2, TMA, cell_ID)
meta_flabsheet <- group_meta %>% select(Accession.Number, Group.label)
names(meta_flabsheet) <- c("TMA_2", "Group")

TMA_ID_fseu <- meta_fseu$TMA_2 %>% as.factor() %>% levels()
TMA_ID_fseu <- TMA_ID_fseu[-1]
TMA_ID_flabsheet <- meta_flabsheet$TMA_2 %>% as.factor() %>% levels()

intersect(TMA_ID_fseu, TMA_ID_flabsheet)
TMA_ID_toChange <- setdiff(TMA_ID_fseu, TMA_ID_flabsheet)
print(TMA_ID_toChange)

# no changes to make

meta_fseu <- merge(meta_fseu, meta_flabsheet, by="TMA_2")
rownames(meta_fseu) <- meta_fseu$cell_ID

obj = AddMetaData(
  object = obj,
  metadata = meta_fseu)


Adult_1_obj <- obj



## Adult_2_obj
obj <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult_2_obj.rds")

first_10_chars <- substr(obj@meta.data$TMA, 1, 10)
first_10_chars <- gsub(" ", "", first_10_chars)
obj@meta.data$TMA_2 <- first_10_chars
obj$cell_ID <- rownames(obj@meta.data)

meta_fseu <- obj@meta.data %>% select(TMA_2, TMA, cell_ID)
meta_flabsheet <- group_meta %>% select(Accession.Number, Group.label)
names(meta_flabsheet) <- c("TMA_2", "Group")

TMA_ID_fseu <- meta_fseu$TMA_2 %>% as.factor() %>% levels()
TMA_ID_fseu <- TMA_ID_fseu[-1]
TMA_ID_flabsheet <- meta_flabsheet$TMA_2 %>% as.factor() %>% levels()

intersect(TMA_ID_fseu, TMA_ID_flabsheet)
TMA_ID_toChange <- setdiff(TMA_ID_fseu, TMA_ID_flabsheet)
print(TMA_ID_toChange)

# no changes to make

meta_fseu <- merge(meta_fseu, meta_flabsheet, by="TMA_2")
rownames(meta_fseu) <- meta_fseu$cell_ID

obj = AddMetaData(
  object = obj,
  metadata = meta_fseu)


Adult_2_obj <- obj


## Save it
SeuObj_list <- list(Adult_1_obj = Adult_1_obj,
  Adult_2_obj = Adult_2_obj,
  Neonate_1_obj = Neonate_1_obj, 
  Neonate_2_obj = Neonate_2_obj)

names(SeuObj_list) <- c("Adult_1", "Adult_2", "Neonate_1", "Neonate_2")


## Go back to the main ppt slides to check wheter the objects still hold after ID handling (it does! )
coord <- GetTissueCoordinates(SeuObj_list[[1]])
plot <- ggplot(coord, aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs() +
  theme_minimal()
print(plot)

coord <- GetTissueCoordinates(SeuObj_list[[2]])
plot <- ggplot(coord, aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs() +
  theme_minimal()
print(plot)

coord <- GetTissueCoordinates(SeuObj_list[[3]])
plot <- ggplot(coord, aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs() +
  theme_minimal()
print(plot)

coord <- GetTissueCoordinates(SeuObj_list[[4]])
plot <- ggplot(coord, aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs() +
  theme_minimal()
print(plot)


# #saveRDS(SeuObj_list, file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult12_Neonate12_seuobj_list.rds")


