## Load packages -------
library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(qs)

## Set global var. -------
options(future.globals.maxSize = 99999 * 1024^2)


## Prepare objects  -------
# Load the Seurat object list
sobj_list <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Adult12_Neonate12_seuobj_list.rds")

coord <- GetTissueCoordinates(sobj_list[[3]])
plot <- ggplot(coord, aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs() +
  theme_minimal()
print(plot)

table(sobj_list[[1]]$Group, useNA = "always")
table(sobj_list[[2]]$Group, useNA = "always")
table(sobj_list[[3]]$Group, useNA = "always")
table(sobj_list[[4]]$Group, useNA = "always")

# This NAs on Adult 1 are the control 
table(sobj_list[[1]]$TMA_2, sobj_list[[1]]$Group, useNA = "always")

# Removing it from our analysis (and every cell that has not gotten label by the TMA IDs from /mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Scripts/Creating_SeuObj.R) 
# Function to subset out NA Group cells

subset_no_NA_Group <- function(obj) {
  subset(obj, TMA %in% "", invert = TRUE) # had to filter by $TMA otherwise it would not go through - might be a data type issue
}

# Apply to all objects in the list
sobj_list <- lapply(sobj_list, subset_no_NA_Group)

# Checking for NAs 
table(sobj_list[[1]]$Group, useNA = "always")
table(sobj_list[[2]]$Group, useNA = "always")
table(sobj_list[[3]]$Group, useNA = "always")
table(sobj_list[[4]]$Group, useNA = "always")



# Define metadata directories
metadata_dirs <- list(
  Adult_1 = "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult1831439_08_11_2024_9_49_46_243/flatFiles/831439_IND_CR_FFPE_TMA_Adult1_STJ_N_R1",
  Neonate_1 = "/mnt/scratch2/CosMX/Diego/DiegoRSVNeonate1831439_08_11_2024_11_18_07_116/flatFiles/831439_IND_CR_FFPE_TMA_Neonate1_STJ_N_R1",
  Adult_2 = "/mnt/scratch2/CosMX/Diego/DiegoRSVAdult2831439_08_11_2024_12_58_14_961/flatFiles/831439_IND_CR_FFPE_TMA_Adult2_STJ_N_R1",
  Neonate_2 = "/mnt/scratch2/CosMX/Diego/DIegoRSVNeonate2831439_08_11_2024_13_26_50_655/flatFiles/831439_IND_CR_FFPE_TMA_Neonate2_STJ_N_R1"
)

# Load metadata files and merge them into the respective Seurat objects
for (slide_name in names(metadata_dirs)) {
  # Load metadata for the current slide
  metadata_file <- list.files(path = metadata_dirs[[slide_name]], pattern = "metadata.*\\.csv\\.gz$", full.names = TRUE)
  metadata <- fread(metadata_file) %>% as.data.frame()
  
  # Add Slide and cell columns to the metadata
  metadata$Slide <- slide_name
  metadata$cell <- paste0(as.character(metadata$cell_ID), "_", metadata$fov, "_", metadata$Slide)
  
  # Set rownames of the metadata to match cell identifiers
  rownames(metadata) <- metadata$cell
  
  # Filter to seurat obj
  seurat_obj <- sobj_list[[slide_name]]
  seurat_obj@meta.data$Slide <- slide_name
  seurat_obj@meta.data$cell <- paste0(as.character(seurat_obj@meta.data$cell_ID), #"_", seurat_obj@meta.data$fov,
                                      "_", seurat_obj@meta.data$Slide)
  #rownames(seurat_obj@meta.data) <- seurat_obj@meta.data$cell
  
  metadata <- metadata[metadata$cell %in% seurat_obj@meta.data$cell, ]
  
  # Merge the metadata into the corresponding Seurat object
  sobj_list[[slide_name]] <- AddMetaData(seurat_obj, metadata)
}



## Check spatial coord. across slides -------
# Suppose your list is named sobj_list and each object has a name
obj_names <- names(sobj_list)

# Extract coordinates and add ID column for each object
coord_list <- lapply(seq_along(sobj_list), function(i) {
  coords <- as.data.frame(GetTissueCoordinates(sobj_list[[i]]))
  coords$ID <- obj_names[i]
  coords
})

# Combine all into one data frame
all_coords <- do.call(rbind, coord_list)

library(ggplot2)

ggplot(all_coords, aes(x = x, y = y, color = ID)) +
  geom_point(alpha = 0.7) +
  labs(title = "Tissue Coordinates by Object", x = "X", y = "Y", color = "Object") +
  theme_minimal()



## Conclusion: We can’t do spatial analysis (niche specifically) across different slides 




## Call Cell Types -------

# Load reference cell type obj
load("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/Lung_MCA.RData")
metadata
cellGroups
profile_matrix[1:10, 1:10] # row: mouse genes; col: cell type; content: cell type "fraction" by gene
dim(profile_matrix)


# Run cell type
# Define the function to process each Sobj
call_celltype <- function(Sobj) {
  # Run insitutypeML
  sup <- insitutypeML(
    x = t(as.matrix(Sobj@assays$Nanostring$counts)),
    neg = Sobj@meta.data$nCount_negprobes,
    # cohort = cohort, # Uncomment if you have a cohort variable
    reference_profiles = as.matrix(profile_matrix)
  )
  
  # Assign clusters
  Sobj@meta.data$insitutype_cluster <- NA
  Sobj@meta.data$insitutype_cluster <- 'other'
  Sobj$insitutype_cluster[names(sup$clust)] <- sup$clust
  names(Sobj@meta.data)[names(Sobj@meta.data) == "insitutype_cluster"] <- 
    paste0("insitutype_cluster_", "Superv", "_celltypes")
  
  # Assign major cell types
  Sobj@meta.data$major_celltype <- "Other"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Alveolar.bipotent.progenitor", ]$major_celltype <- "Alveolar_bipotent_progenitor"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% c("Alveolar.macrophage.Ear2.high", "Alveolar.macrophage.Pclaf.high"), ]$major_celltype <- "Alveolar_macrophage"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "AT1.cell", ]$major_celltype <- "AT1_cell"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "AT2.cell", ]$major_celltype <- "AT2_cell"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% c("B.cell", "IgA.producing.B.cell", "NK.cell", "T.cell.Cd8b1.high"), ]$major_celltype <- "Lymphocyte"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% c("Basophil", "Eosinophil", "Neutrophil"), ]$major_celltype <- "Granulocyte"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% c("Ciliated.cell", "Clara.cell"), ]$major_celltype <- "Epithelial_cell"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% c("Dendritic.cell.Gngt2.high", "Dendritic.cell.H2.M2.high", "Dendritic.cell.Mgl2.high", "Dendritic.cell.Naaa.high", "Dendritic.cell.Tubb5.high", "Plasmacytoid.dendritic.cell"), ]$major_celltype <- "Dendritic_cell"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% c("Endothelial.cell.Kdr.high", "Endothelial.cell.Tmem100.high", "Endothelial.cell.Vwf.high"), ]$major_celltype <- "Endothelial_cell"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Interstitial.macrophage", ]$major_celltype <- "Interstitial_macrophage"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Monocyte.progenitor", ]$major_celltype <- "Monocyte_progenitor"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes %in% c("Stromal.cell.Acta2.high", "Stromal.cell.Dcn.high", "Stromal.cell.Inmt.high"), ]$major_celltype <- "Stromal_cell"
  Sobj@meta.data[Sobj@meta.data$insitutype_cluster_Superv_celltypes == "Nuocyte", ]$major_celltype <- "Nuocyte"
  
  # Optionally print tables (remove if not needed)
  print(table(Sobj@meta.data$insitutype_cluster_Superv_celltypes, useNA = "always"))
  print(table(Sobj@meta.data$major_celltype, useNA = "always"))
  
  return(Sobj)
}

# Apply the function to each object in the list
sobj_list <- lapply(sobj_list, call_celltype)



## Call Niches -------
# Function to run niche analysis on a single Seurat object
run_niche_analysis <- function(SeuObj) {
  DefaultAssay(SeuObj) <- "Nanostring"
  SeuObj <- BuildNicheAssay(
    object = SeuObj,
    fov = "FOV",
    group.by = "major_celltype",
    niches.k = 5,
    neighbors.k = 20,
    cluster.name = "niches_5"
  )
  return(SeuObj)
}

# Apply to all objects in the list
sobj_list <- lapply(sobj_list, run_niche_analysis)



## Visualization ------- 

cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- sobj_list[[1]]$major_celltype %>% table() %>% names()
cols <- setNames(cols[1:length(factor_levels)], factor_levels)

library(viridis)
niche_cols <- viridis(5)
factor_levels <- sobj_list[[1]]$niches_5 %>% table() %>% names()
niche_cols <- setNames(niche_cols[1:length(factor_levels)], factor_levels)


## Plotting by Slide  
# Cell Type Proportion by group (prop. barplot)
library(dittoSeq)
dittoBarPlot(sobj_list[[3]], 
             "major_celltype", 
             group.by = "Group", 
             color.panel = cols,
             main = "")

# Cell Type Proportion by group (prop. barplot)
library(dittoSeq)
dittoBarPlot(sobj_list[[3]], 
             "niches_5", 
             group.by = "Group", 
             color.panel = niche_cols,
             main = "")

library(dittoSeq)
dittoBarPlot(sobj_list[[3]], 
             "niches_5", 
             group.by = "TMA_2", 
             split.by = "Group",
             color.panel = niche_cols,
             main = "")





# Niche ~ Cell Type by group (prop. barplot)
library(dittoSeq)
dittoBarPlot(sobj_list[[3]], 
             "major_celltype", 
             group.by = "niches_5", 
             split.by = "Group",
             color.panel = cols,
             main = "Neonate_1")

# loop it 
library(dittoSeq)
library(patchwork) # or use cowplot::plot_grid

# Iterate and create plots
plots <- lapply(seq_along(sobj_list), function(i) {
  dittoBarPlot(
    sobj_list[[i]],
    "major_celltype",
    group.by = "niches_5",
    split.by = "Group",
    color.panel = cols,
    main = sobj_list[i] %>% names()
  )
})

# Combine in a 2x2 grid
(plots[[1]] | plots[[2]]) /
  (plots[[3]] | plots[[4]])




## Top CellType Markers Across Cores/Slides -------
## Do cell types look right across slides ? 



## Matching Niches Across Slides  (TMA_ID strange mismatch) -------  -------

ImageDimPlot(subset(sobj_list[[3]], Group %in% "Neonate IFN and RSV reinfected"),
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


ImageDimPlot(subset(sobj_list[[3]], Group %in% "Neonate RSV infected (NO IFN)"),
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


table(sobj_list[[3]]$TMA_2,
      sobj_list[[3]]$Group)


table(sobj_list[[3]]$TMA_2)



# Triple checking on TMA_ID and actual TMA cores in the data

ImageDimPlot(sobj_list[[3]],
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

coord <- GetTissueCoordinates(sobj_list[[1]])
plot <- ggplot(coord, aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs() +
  theme_minimal()
print(plot)

coord <- GetTissueCoordinates(sobj_list[[2]])
plot <- ggplot(coord, aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs() +
  theme_minimal()
print(plot)

coord <- GetTissueCoordinates(sobj_list[[3]])
plot <- ggplot(coord, aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs() +
  theme_minimal()
print(plot)

coord <- GetTissueCoordinates(sobj_list[[4]])
plot <- ggplot(coord, aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs() +
  theme_minimal()
print(plot)

Neonate_1_obj <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_1/Objects/Neonate_1_obj.rds")
coord <- GetTissueCoordinates(Neonate_1_obj)
plot <- ggplot(coord, aes(x = x, y = y)) +
  geom_point(size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs() +
  theme_minimal()
print(plot)

library(dbscan)
result <- dbscan(coord[, c("x", "y")], eps = 10, minPts = 20)
table(result$cluster)

# Conclusion: the number of TMA cores are fewer compared to Hannah's ppt slides becasue I've filtered low quality cells already. 
# However, we still can't have more TMA cores then TMA_IDs


## Labeling single TMA cores
FOV_catcheR <- function(
    seurat_obj,
    n_clusters = 5,
    n_nstart = 25, # to run the algorithm multile times - it gets better
    save_path = NULL,
    plot_clusters = TRUE,
    seed = 42,
    object_name = "",
    export_seurat_obj = FALSE,   # NEW ARGUMENT
    # Additional parameters for future features
    feature1 = NULL,
    feature2 = NULL,
    feature3 = NULL
) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Get tissue coordinates
  data <- GetTissueCoordinates(seurat_obj)
  names(data) <- c("x", "y", "id")
  data <- data[!duplicated(data$id), ]
  
  # Perform K-means clustering
  kmeans_result <- kmeans(data[, c("x", "y")], centers = n_clusters, nstart = n_nstart)
  
  # Add cluster information to the data
  data$spatial_cluster <- kmeans_result$cluster
  
  # Add spatial cluster information to the Seurat object
  seurat_obj@meta.data$spatial_cluster <- data$spatial_cluster[match(colnames(seurat_obj), data$id)]
  
  # Subset the Seurat object based on spatial clusters
  seurat_subsets <- list()
  for (cluster in 1:n_clusters) {
    cells_in_cluster <- data$id[data$spatial_cluster == cluster]
    seurat_subsets[[paste0("cluster_", cluster)]] <- subset(seurat_obj, cells = cells_in_cluster)
  }
  
  # Save Seurat objects as .qs files if save_path is provided AND export_seurat_obj is TRUE
  if (!is.null(save_path) && export_seurat_obj) {
    dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
    for (i in 1:n_clusters) {
      file_path <- file.path(save_path, paste0("Object_", i, object_name, ".qs"))
      qsave(seurat_subsets[[i]], file_path)
      cat("Saved", paste0("Object_", i), "to", file_path, "\n")
    }
  }
  
  # Plot clusters if requested
  if (plot_clusters) {
    plot <- ggplot(result$clustered_data, aes(x = x, y = y, color = factor(spatial_cluster))) +
      geom_point(size = 0.5) +
      #scale_color_brewer(palette = "Set1") +
      labs(title = "Spatial Clusters", color = "Cluster") +
      theme_minimal()
    print(plot)
  }
  
  # Return results
  return(list(
    clustered_data = data,
    seurat_obj = seurat_obj,
    plot = plot
    #kmeans_result = kmeans_result
  ))
}



ImageDimPlot(sobj_list[[3]],
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


result <- FOV_catcheR(
  seurat_obj = sobj_list[[3]],
  n_clusters = 25,
  n_nstart = 1000,
  export_seurat_obj = FALSE,
  plot_clusters = TRUE
)

result$clustered_data %>% dim()
dim(sobj_list[[3]])
table(result$clustered_data$spatial_cluster, useNA = "always")

plot <- ggplot(result$clustered_data, aes(x = x, y = y, color = factor(spatial_cluster))) +
  geom_point(size = 0.5) +
  #scale_color_brewer(palette = "Set1") +
  labs(title = "Spatial Clusters", color = "Cluster") +
  theme_minimal()
print(plot)


## TMA_ID by TMA_core anno. check 
table(sobj_list[[3]]$TMA) %>% length()
df <- data.frame(table(sobj_list[[3]]$TMA, 
      sobj_list[[3]]$fov))

names(df) <- c("TMA_ID", "FOV", "nCells")
df %>% head()

ggplot(df, aes(x = FOV, y = nCells)) +
  geom_point() +
  facet_wrap(~ TMA_ID, scales = "free_x") +
  labs(title = "Cells per FOV for each TMA_ID", x = "FOV", y = "Number of Cells")

# try TMA_2
df <- data.frame(table(sobj_list[[3]]$TMA_2, 
                       sobj_list[[3]]$fov))

names(df) <- c("TMA_ID", "FOV", "nCells")
df %>% head()

ggplot(df, aes(x = FOV, y = nCells)) +
  geom_point() +
  facet_wrap(~ TMA_ID, scales = "free_x") +
  labs(title = "Cells per FOV for each TMA_ID", x = "FOV", y = "Number of Cells")


test <- subset(df, TMA_ID %in% "Ap16-2656")  
n_FOV <- table(df$FOV) %>% names() %>% length()
ggplot(test, aes(x = FOV, y = nCells)) +
  geom_point() +
  #facet_wrap(~ TMA_ID, scales = "free_x") +
  labs(title = "Cells per FOV for each TMA_ID", x = "FOV", y = "Number of Cells")


ImageDimPlot(subset(sobj_list[[3]], TMA_2 %in% "Ap16-2656"),
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



TMA_IDs <- table(df$TMA_ID) %>% names()
test <- subset(df, TMA_ID %in% TMA_IDs[i])  
test$nCells #if there's m




library(dplyr)

# Function to count intervals of nonzero values in a numeric vector
count_intervals <- function(x) {
  # x is a numeric vector (e.g. nCells)
  # Create a logical vector: TRUE if nonzero, FALSE if zero
  nonzero <- x != 0
  
  # Identify runs of TRUE values separated by FALSE
  # We count how many times a TRUE run starts
  # A run start is where nonzero is TRUE and previous is FALSE or NA (start of vector)
  
  starts <- which(nonzero & (c(TRUE, head(!nonzero, -1))))
  
  length(starts)
}

# Get unique TMA_IDs
TMA_IDs <- unique(df$TMA_ID)

# Initialize a vector to store counts
interval_counts <- numeric(length(TMA_IDs))
names(interval_counts) <- TMA_IDs

# Loop over each TMA_ID and count intervals
for (i in seq_along(TMA_IDs)) {
  test <- subset(df, TMA_ID == TMA_IDs[i])
  interval_counts[i] <- count_intervals(test$nCells)
}

# Convert to a dataframe for easier viewing
interval_counts_df <- data.frame(
  TMA_ID = TMA_IDs,
  Interval_Count = interval_counts
)

print(interval_counts_df)


## That's good vis !! 
df <- table(sobj_list[[3]]$TMA_2,
      sobj_list[[3]]$fov) %>% as.matrix()
library(pheatmap)
pheatmap(
  df,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)


df <- table(sobj_list[[3]]$TMA_2,
            sobj_list[[3]]$fov) %>% as.data.frame()
names(df) <- c("TMA_ID", "FOV", "nCells")

df_zero <- filter(df, nCells == 0) 
df_not_zero <- filter(df, nCells > 0) 

setdiff(df_zero$TMA_ID, df_not_zero$TMA_ID)
setdiff(df_not_zero$TMA_ID, df_zero$TMA_ID)




## Matching Niches Across Slides (this automated approach hasn't work out yet----------
Slide_names <- names(sobj_list)
metadata_1 <- sobj_list[[1]]@meta.data
metadata_1$niches_5_Slide_1 <- metadata_1$niches_5

df1 <- metadata_1[, c("Slide", 
             "niches_5",
             "major_celltype")]

metadata_2 <- sobj_list[[2]]@meta.data
metadata_2$niches_5_Slide_2 <- metadata_2$niches_5
df2 <- metadata_2[, c("Slide", 
                      "niches_5",
                      "major_celltype")]

metadata_3 <- sobj_list[[3]]@meta.data
metadata_3$niches_5_Slide_3 <- metadata_3$niches_5
df3 <- metadata_3[, c("Slide", 
                      "niches_5",
                      "major_celltype")]

metadata_4 <- sobj_list[[4]]@meta.data
metadata_4$niches_5_Slide_4 <- metadata_4$niches_5
df4 <- metadata_4[, c("Slide", 
                      "niches_5",
                      "major_celltype")]


metadata <- plyr::rbind.fill(metadata_1, metadata_2,
                 metadata_3, metadata_4)

metadata[, c("Slide", 
             "niches_5_Slide_1", 
             "niches_5_Slide_2",
             "niches_5_Slide_3",
             "niches_5_Slide_4",
             "major_celltype")] %>% head()

niche_metadata <- metadata[, c("Slide", 
                               "Group",
                               "niches_5_Slide_1", 
                               "niches_5_Slide_2",
                               "niches_5_Slide_3",
                               "niches_5_Slide_4",
                               "major_celltype")] 

head(niche_metadata)


### Get celltype proportion by Niche split by Group
library(tidyr)
library(dplyr)

# Assuming your data frame is named df
result <- niche_metadata %>%
  pivot_longer(
    cols = starts_with("niches_5_Slide_"),
    names_to = "niche_slide",
    values_to = "niche"
  ) %>%
  filter(!is.na(niche)) %>%
  group_by(Group, niche, major_celltype) %>%
  dplyr::summarise(count = n(), .groups = "drop_last") %>%
  mutate(
    total = sum(count),
    proportion = count / total
  ) %>%
  ungroup() %>%
  select(Group, niche, major_celltype, proportion) %>% 
  data.frame()

result[result$niche %in% 1, ] %>% head()
 

### Heatmap of 1 niche at a time
niche_1 <- result[result$niche %in% 1 & , ] %>%
  # 1. Reshape to matrix
  unite("group_niche", Group, niche, sep = "|") %>%
  pivot_wider(names_from = major_celltype, values_from = proportion, values_fill = 0) %>%
  column_to_rownames("group_niche") 

niche_1 <- t(niche_1)

# pheatmap(niche_1)

cor_mat <- cor(niche_1, method = "pearson")
pheatmap(cor_mat,
         # clustering_distance_rows = "correlation",
         # clustering_distance_cols = "correlation",
         main = "Correlation Matrix of Niches by Cell Type Proportions Across Slides (Niche 1)")



### Heatmap of all niches 
# Load necessary libraries
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)

# 1. Prepare the matrix of cell type proportions per niche
niche_all <- result %>%
  unite("group_niche", Group, niche, sep = "|") %>%
  pivot_wider(names_from = major_celltype, values_from = proportion, values_fill = 0) %>%
  column_to_rownames("group_niche") %>%
  t()

# 2. Compute correlation matrix between niches
cor_mat <- cor(niche_all, method = "pearson")

# 3. Create annotation dataframe for columns
niche_names <- colnames(cor_mat)
annotation_df <- as.data.frame(do.call(rbind, strsplit(niche_names, "\\|")))
colnames(annotation_df) <- c("Group", "Niche")
rownames(annotation_df) <- niche_names
annotation_df$Group <- NULL
# 4. Assign viridis colors to Groups and Niches
group_levels <- unique(annotation_df$Group)
niche_levels <- unique(annotation_df$Niche)
group_colors <- setNames(viridis(length(group_levels)), group_levels)
niche_colors <- setNames(viridis(length(niche_levels)), niche_levels)
ann_colors <- list(Group = group_colors, Niche = niche_colors)

# 5. Plot correlation heatmap with viridis-labeled annotation bars
pheatmap(
  # cor_mat[c(1, 6:40), c(1, 6:40)]
  cor_mat[c(1, 6:40), c(1, 6:40)],
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  main = "Correlation Matrix of Niches by Cell Type Proportions Across Slides (All 5 Niches)"
)






### Get celltype proportion by Niche split by Slide
library(tidyr)
library(dplyr)

# Assuming your data frame is named df
result <- niche_metadata %>%
  pivot_longer(
    cols = starts_with("niches_5_Slide_"),
    names_to = "niche_slide",
    values_to = "niche"
  ) %>%
  filter(!is.na(niche)) %>%
  group_by(Slide, niche, major_celltype) %>%
  dplyr::summarise(count = n(), .groups = "drop_last") %>%
  mutate(
    total = sum(count),
    proportion = count / total
  ) %>%
  ungroup() %>%
  select(Slide, niche, major_celltype, proportion) %>% 
  data.frame()

library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)

# 1. Prepare the matrix of cell type proportions per niche
niche_all <- result %>%
  unite("group_niche", Slide, niche, sep = "|") %>%
  pivot_wider(names_from = major_celltype, values_from = proportion, values_fill = 0) %>%
  column_to_rownames("group_niche") %>%
  t()

niche_all <- niche_all[, !colnames(niche_all) %in% c("Adult_1|2", 
                                           "Adult_1|3",
                                           "Adult_1|4", 
                                           "Adult_1|5")]

# 2. Compute correlation matrix between niches
cor_mat <- cor(niche_all, method = "pearson")

# 3. Create annotation dataframe for columns
niche_names <- colnames(cor_mat)
annotation_df <- as.data.frame(do.call(rbind, strsplit(niche_names, "\\|")))
colnames(annotation_df) <- c("Slide", "Niche")
rownames(annotation_df) <- niche_names
annotation_df$Slide <- NULL
annotation_df$Test_Gru_Nich <- "other"
annotation_df[rownames(annotation_df) %in% "Adult_1|1", ]$Test_Gru_Nich <- "test niche"

# 4. Assign viridis colors to Groups and Niches
group_levels <- unique(annotation_df$Slide)
group_colors <- setNames(viridis(length(group_levels)), group_levels)

niche_levels <- unique(annotation_df$Niche)
niche_colors <- setNames(viridis(length(niche_levels)), niche_levels)

Test_Gru_Nich_levels <- unique(annotation_df$Test_Gru_Nich)
Test_Gru_Nich_colors <- setNames(c("red", "gray"), Test_Gru_Nich_levels)


ann_colors <- list(Slide = group_colors, Niche = niche_colors,
                   Test_Gru_Nich = Test_Gru_Nich_colors)

# 5. Plot correlation heatmap with viridis-labeled annotation bars
pheatmap(
  cor_mat[, ],
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  main = ""
)


# qsave(sobj_list, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj_list_celltype_niche.qs")

#

#

#

#

# 

#

#

#

#

#


### -------------------------------------------------------
### Ignoring the matching niches across TMAs/slides for now 
### Going down the analysis with what we have 

# Global var --------
sobj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj_list_celltype_niche.qs")

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

# Setting color map 
cols <-
  c(
    '#8DD3C7', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
    '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A',
    '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999'
  )

factor_levels <- sobj_list[[1]]$major_celltype %>% table() %>% names()
celltype_cols <- setNames(cols[1:length(factor_levels)], factor_levels)

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


table(sobj_list[[1]]$Group)
table(sobj_list[[2]]$Group)
table(sobj_list[[3]]$Group)
table(sobj_list[[4]]$Group)

sobj_list[[1]]@meta.data$niches_5  <- as.factor(sobj_list[[1]]@meta.data$niches_5 )
sobj_list[[1]]@meta.data$niches_5 <- factor(sobj_list[[1]]@meta.data$niches_5, levels = as.factor(c(1, 2, 3, 4, 5)))

sobj_list[[2]]@meta.data$niches_5  <- as.factor(sobj_list[[2]]@meta.data$niches_5 )
sobj_list[[2]]@meta.data$niches_5 <- factor(sobj_list[[2]]@meta.data$niches_5, levels = as.factor(c(1, 2, 3, 4, 5)))

sobj_list[[3]]@meta.data$niches_5  <- as.factor(sobj_list[[3]]@meta.data$niches_5 )
sobj_list[[3]]@meta.data$niches_5 <- factor(sobj_list[[3]]@meta.data$niches_5, levels = as.factor(c(1, 2, 3, 4, 5)))

sobj_list[[4]]@meta.data$niches_5  <- as.factor(sobj_list[[4]]@meta.data$niches_5 )
sobj_list[[4]]@meta.data$niches_5 <- factor(sobj_list[[4]]@meta.data$niches_5, levels = as.factor(c(1, 2, 3, 4, 5)))





## Comparison #1 -----

comparisons$`Comparison 1`


plot_spatial_and_bar <- function(
    obj,
    group_filter,
    group_by,
    border_color = "white",
    fov = NULL,
    border_size = 0.1,
    cols = NULL,
    size = 0.8,
    axes = FALSE,
    bar_var,
    bar_group_by,
    bar_colors = NULL,
    bar_main = ""
) {
  # Subset the object
  obj_sub <- subset(obj, Group %in% group_filter)
  
  # ImageDimPlot
  p1 <- ImageDimPlot(
    object = obj_sub,
    group.by = group_by,
    border.color = border_color,
    fov = fov,
    border.size = border_size,
    cols = cols,
    size = size,
    axes = axes
  ) + ggtitle("")
  
  # dittoBarPlot
  p2 <- dittoBarPlot(
    obj_sub,
    bar_var,
    group.by = bar_group_by,
    color.panel = bar_colors,
    main = bar_main
  )
  
  return(list(spatial_plot = p1, bar_plot = p2))
}



Obj_1_plots <- plot_spatial_and_bar(
  obj = sobj_list[[1]],
  group_filter = "Adult RSV infected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_2_plots <- plot_spatial_and_bar(
  obj = sobj_list[[2]],
  group_filter = "Adult RSV infected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_3_plots <- plot_spatial_and_bar(
  obj = sobj_list[[3]],
  group_filter = "Neonate IFN and RSV infected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_4_plots <- plot_spatial_and_bar(
  obj = sobj_list[[4]],
  group_filter = "Neonate IFN and RSV infected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)


Obj_1_plots$spatial_plot
Obj_1_plots$bar_plot

Obj_2_plots$spatial_plot
Obj_2_plots$bar_plot


Obj_3_plots$spatial_plot
Obj_3_plots$bar_plot

Obj_4_plots$spatial_plot
Obj_4_plots$bar_plot

### Matching niches manually 
## Across Slides - same Group
## sobj_list[[3]] is the based for what we're calling "niches 1,2,3,4,and 5"
sobj_list[[4]]@meta.data$niches_5 <- mapvalues(sobj_list[[4]]@meta.data$niches_5, from = as.factor(c(1, 2, 3, 4, 5)), to = as.factor(c(4, 5, 1, 2, 3)))

sobj_list[[4]]@meta.data$niches_5 <- factor(sobj_list[[4]]@meta.data$niches_5, levels = as.factor(c(1, 2, 3, 4, 5)))

## based on sobj_list[[1]]
sobj_list[[2]]@meta.data$niches_5 <- mapvalues(sobj_list[[2]]@meta.data$niches_5, from = as.factor(c(1, 2, 3, 4, 5)), to = as.factor(c(1, 2, 5, 3, 4)))

sobj_list[[2]]@meta.data$niches_5 <- factor(sobj_list[[2]]@meta.data$niches_5, levels = as.factor(c(1, 2, 3, 4, 5)))



## Across Groups 
## based on sobj_list[[1]]
sobj_list[[3]]@meta.data$niches_5 <- mapvalues(sobj_list[[3]]@meta.data$niches_5, from = as.factor(c(1, 2, 3, 4, 5)), to = as.factor(c(1, 4, 3, 5, 2)))

sobj_list[[3]]@meta.data$niches_5 <- factor(sobj_list[[3]]@meta.data$niches_5, levels = as.factor(c(1, 2, 3, 4, 5)))


## based on sobj_list[[2]] == sobj_list[[1]]
sobj_list[[4]]@meta.data$niches_5 <- mapvalues(sobj_list[[4]]@meta.data$niches_5, from = as.factor(c(1, 2, 3, 4, 5)), to = as.factor(c(1, 4, 3, 5, 2)))

sobj_list[[4]]@meta.data$niches_5 <- factor(sobj_list[[4]]@meta.data$niches_5, levels = as.factor(c(1, 2, 3, 4, 5)))


## Adjusting niche 3 from Neonates [[3]] and [[4]] to niche 1 
sobj_list[[3]]@meta.data$niches_5 <- mapvalues(sobj_list[[3]]@meta.data$niches_5, from = as.factor(c(1, 2, 3, 4, 5)), to = as.factor(c(3, 2, 1, 4, 5)))

sobj_list[[3]]@meta.data$niches_5 <- factor(sobj_list[[3]]@meta.data$niches_5, levels = as.factor(c(1, 2, 3, 4, 5)))

sobj_list[[4]]@meta.data$niches_5 <- mapvalues(sobj_list[[4]]@meta.data$niches_5, from = as.factor(c(1, 2, 3, 4, 5)), to = as.factor(c(3, 2, 1, 4, 5)))

sobj_list[[4]]@meta.data$niches_5 <- factor(sobj_list[[4]]@meta.data$niches_5, levels = as.factor(c(1, 2, 3, 4, 5)))




Obj_1_plots <- plot_spatial_and_bar(
  obj = sobj_list[[1]],
  group_filter = "Adult RSV infected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_2_plots <- plot_spatial_and_bar(
  obj = sobj_list[[2]],
  group_filter = "Adult RSV infected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_3_plots <- plot_spatial_and_bar(
  obj = sobj_list[[3]],
  group_filter = "Neonate IFN and RSV infected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_4_plots <- plot_spatial_and_bar(
  obj = sobj_list[[4]],
  group_filter = "Neonate IFN and RSV infected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)


Obj_1_plots$spatial_plot
Obj_1_plots$bar_plot

Obj_2_plots$spatial_plot
Obj_2_plots$bar_plot


Obj_3_plots$spatial_plot
Obj_3_plots$bar_plot

Obj_4_plots$spatial_plot
Obj_4_plots$bar_plot



### Merge all the objects into one for "single cell vis"
# Start with the first object
options(future.globals.maxSize = 99999 * 1024^2)

# All_Sobj <- merge(x = sobj_list[[1]], y = sobj_list[-1])
# All_Sobj <- JoinLayers(All_Sobj)

All_Sobj <- sobj_list[[1]]

# Loop through the remaining objects and merge them one by one
for (i in 2:length(sobj_list)) {
  All_Sobj <- merge(All_Sobj, y = sobj_list[[i]])
}

DefaultAssay(All_Sobj) <- "Nanostring"
All_Sobj <- JoinLayers(All_Sobj)

# qsave(All_Sobj, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_merged_Sobj.qs")


## Plot Dotplot 
library(presto)
All_Sobj@assays$RNA <- All_Sobj@assays$Nanostring
DefaultAssay(All_Sobj) <- "RNA"
All_Sobj <- All_Sobj %>% NormalizeData()

markers <- wilcoxauc(All_Sobj, 'major_celltype')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 5) %>% 
  ungroup() -> top_markers
table(top_markers$group)

library(ggplot2)
library(viridis)
Idents(All_Sobj) <- "major_celltype"
DotPlot(All_Sobj, 
        #features = unique(top_markers$feature)
        features =  unique(c(top_markers$feature))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Cell Typer Markers")


# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 20) %>% 
  ungroup() -> top_markers
table(top_markers$group)

library(ggplot2)
library(viridis)
Idents(All_Sobj) <- "major_celltype"
DotPlot(All_Sobj, 
        #features = unique(top_markers$feature)
        features =  unique(c(top_markers[top_markers$group %in% "Lymphocyte", ]$feature))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Lymphocyte top markers ONLY")


## Plot Cell Prop 
library(dittoSeq)
dittoBarPlot(All_Sobj, 
             "major_celltype", 
             group.by = "Group", 
             color.panel = celltype_cols,
             main = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) + ggtitle("")



library(dittoSeq)
dittoBarPlot(subset(All_Sobj, Group %in% 
                      c(comparisons$`Comparison 1`$group_1,
                        comparisons$`Comparison 1`$group_2)), 
             "niches_5", 
             group.by = "TMA", 
             split.by = "Group",
             color.panel = niche_cols,
             main = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) + ggtitle("")


library(dittoSeq)
dittoBarPlot(subset(All_Sobj, Group %in% 
                      c(comparisons$`Comparison 1`$group_1,
                        comparisons$`Comparison 1`$group_2)), 
             "major_celltype", 
             group.by = "TMA", 
             split.by = "Group",
             color.panel = celltype_cols,
             main = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) + ggtitle("")



### Highlight niche 1 cell types from Cropping 
obj <- subset(sobj_list[[1]], TMA %in% "Ap16-3779 2") #pick just one TMA
DefaultBoundary(obj@images$FOV) <- "centroids"

data <- GetTissueCoordinates(obj)
names(data) <- c("x", "y", "id")
data <- data[!duplicated(data$id), ]

obj@meta.data$cell <- rownames(obj@meta.data)

dim(data) #
data <- merge(data, obj@meta.data[, c("niches_5", "cell")], by.x = "id", by.y = "cell")  # Different column names
dim(data) #

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
    ggplot(data, aes(x = x, y = y, color = factor(niches_5))) +
      geom_point(size = 0.5) +  # Adjust point size as needed
      #scale_color_brewer(palette = "Paired") +  # Use the 'Paired' color palette
      scale_color_manual(values = as.vector(niche_cols)) +  # Custom colors
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

DefaultBoundary(obj@images$FOV) <- "segmentation"
niche_1_cropp_plt <- ImageDimPlot(subset(obj, cells = extracted_coord$id),
                             group.by = 'major_celltype',
                             fov = "FOV",
                             size = 0.5,
                             dark.background = TRUE,
                             cols = celltype_cols,
                             border.size = 0.3) +
  #ggtitle(paste0("pre-treatment")) +
  NoLegend() #+
  #scale_y_reverse() 
print(niche_1_cropp_plt)


DefaultBoundary(obj@images$FOV) <- "segmentation"
niche_1_cropp_plt <- ImageDimPlot(subset(obj, cells = extracted_coord$id),
                                  group.by = 'niches_5',
                                  fov = "FOV",
                                  size = 0.5,
                                  dark.background = TRUE,
                                  cols = niche_cols,
                                  border.size = 0.3) +
  #ggtitle(paste0("pre-treatment")) +
  NoLegend() #+
#scale_y_reverse() 
print(niche_1_cropp_plt)









## Comparison #2 -----


comparisons$`Comparison 2`

table(sobj_list[[1]]$Group)
table(sobj_list[[2]]$Group)
table(sobj_list[[3]]$Group)
table(sobj_list[[4]]$Group)



plot_spatial_and_bar <- function(
    obj,
    group_filter,
    group_by,
    border_color = "white",
    fov = NULL,
    border_size = 0.1,
    cols = NULL,
    size = 0.8,
    axes = FALSE,
    bar_var,
    bar_group_by,
    bar_colors = NULL,
    bar_main = ""
) {
  # Subset the object
  obj_sub <- subset(obj, Group %in% group_filter)
  
  # ImageDimPlot
  p1 <- ImageDimPlot(
    object = obj_sub,
    group.by = group_by,
    border.color = border_color,
    fov = fov,
    border.size = border_size,
    cols = cols,
    size = size,
    axes = axes
  ) + ggtitle("")
  
  # dittoBarPlot
  p2 <- dittoBarPlot(
    obj_sub,
    bar_var,
    group.by = bar_group_by,
    color.panel = bar_colors,
    main = bar_main
  )
  
  return(list(spatial_plot = p1, bar_plot = p2))
}


Obj_1_plots <- plot_spatial_and_bar(
  obj = sobj_list[[1]],
  group_filter = "Adult reinfected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_2_plots <- plot_spatial_and_bar(
  obj = sobj_list[[2]],
  group_filter = "Adult reinfected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_3_plots <- plot_spatial_and_bar(
  obj = sobj_list[[3]],
  group_filter = "Neonate IFN and RSV reinfected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_4_plots <- plot_spatial_and_bar(
  obj = sobj_list[[4]],
  group_filter = "Neonate IFN and RSV reinfected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)


Obj_1_plots$spatial_plot
Obj_1_plots$bar_plot

Obj_2_plots$spatial_plot
Obj_2_plots$bar_plot


Obj_3_plots$spatial_plot
Obj_3_plots$bar_plot

Obj_4_plots$spatial_plot
Obj_4_plots$bar_plot



## Adjusting niche 3 from Neonates [[3]] to niche 1 
sobj_list[[3]]@meta.data$niches_5 <- mapvalues(sobj_list[[3]]@meta.data$niches_5, from = as.factor(c(1, 2, 3, 4, 5)), to = as.factor(c(3, 2, 1, 4, 5)))

sobj_list[[3]]@meta.data$niches_5 <- factor(sobj_list[[3]]@meta.data$niches_5, levels = as.factor(c(1, 2, 3, 4, 5)))



Obj_1_plots <- plot_spatial_and_bar(
  obj = sobj_list[[1]],
  group_filter = "Adult reinfected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_2_plots <- plot_spatial_and_bar(
  obj = sobj_list[[2]],
  group_filter = "Adult reinfected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_3_plots <- plot_spatial_and_bar(
  obj = sobj_list[[3]],
  group_filter = "Neonate IFN and RSV reinfected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_4_plots <- plot_spatial_and_bar(
  obj = sobj_list[[4]],
  group_filter = "Neonate IFN and RSV reinfected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)


Obj_1_plots$spatial_plot
Obj_1_plots$bar_plot

Obj_2_plots$spatial_plot
Obj_2_plots$bar_plot


Obj_3_plots$spatial_plot
Obj_3_plots$bar_plot

Obj_4_plots$spatial_plot
Obj_4_plots$bar_plot


### All the objects into one for "single cell vis"
All_Sobj <- qread('/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_merged_Sobj.qs')


## Plot Dotplot 
library(presto)
All_Sobj@assays$RNA <- All_Sobj@assays$Nanostring
DefaultAssay(All_Sobj) <- "RNA"
All_Sobj <- All_Sobj %>% NormalizeData()

markers <- wilcoxauc(All_Sobj, 'major_celltype')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 5) %>% 
  ungroup() -> top_markers
table(top_markers$group)

library(ggplot2)
library(viridis)
Idents(All_Sobj) <- "major_celltype"
DotPlot(All_Sobj, 
        #features = unique(top_markers$feature)
        features =  unique(c(top_markers$feature))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Cell Typer Markers")


# # Filter DEGs based only in FC
# markers %>%
#   group_by(group) %>%
#   dplyr::filter(logFC > 0 & padj <= 0.05) %>%
#   # slice_head(n = 10) %>%
#   dplyr::top_n(wt = logFC, n = 20) %>% 
#   ungroup() -> top_markers
# table(top_markers$group)
# 
# library(ggplot2)
# library(viridis)
# Idents(All_Sobj) <- "major_celltype"
# DotPlot(All_Sobj, 
#         #features = unique(top_markers$feature)
#         features =  unique(c(top_markers[top_markers$group %in% "Lymphocyte", ]$feature))
# ) +
#   geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
#   scale_colour_viridis(option="magma") +
#   guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Lymphocyte top markers ONLY")


## Plot Cell Prop 
library(dittoSeq)
dittoBarPlot(All_Sobj, 
             "major_celltype", 
             group.by = "Group", 
             color.panel = celltype_cols,
             main = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) + ggtitle("")

library(dittoSeq)
dittoBarPlot(subset(All_Sobj, Group %in% 
                      c(comparisons$`Comparison 2`$group_1,
                        comparisons$`Comparison 2`$group_2)), 
             "niches_5", 
             group.by = "TMA", 
             split.by = "Group",
             color.panel = niche_cols,
             main = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) + ggtitle("")


library(dittoSeq)
dittoBarPlot(subset(All_Sobj, Group %in% 
                      c(comparisons$`Comparison 2`$group_1,
                        comparisons$`Comparison 2`$group_2)), 
             "major_celltype", 
             group.by = "TMA", 
             split.by = "Group",
             color.panel = celltype_cols,
             main = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) + ggtitle("")









## Comparison #3 -----


comparisons$`Comparison 3`

table(sobj_list[[1]]$Group)
table(sobj_list[[2]]$Group)
table(sobj_list[[3]]$Group)
table(sobj_list[[4]]$Group)



plot_spatial_and_bar <- function(
    obj,
    group_filter,
    group_by,
    border_color = "white",
    fov = NULL,
    border_size = 0.1,
    cols = NULL,
    size = 0.8,
    axes = FALSE,
    bar_var,
    bar_group_by,
    bar_colors = NULL,
    bar_main = ""
) {
  # Subset the object
  obj_sub <- subset(obj, Group %in% group_filter)
  
  # ImageDimPlot
  p1 <- ImageDimPlot(
    object = obj_sub,
    group.by = group_by,
    border.color = border_color,
    fov = fov,
    border.size = border_size,
    cols = cols,
    size = size,
    axes = axes
  ) + ggtitle("")
  
  # dittoBarPlot
  p2 <- dittoBarPlot(
    obj_sub,
    bar_var,
    group.by = bar_group_by,
    color.panel = bar_colors,
    main = bar_main
  )
  
  return(list(spatial_plot = p1, bar_plot = p2))
}


Obj_1_plots <- plot_spatial_and_bar(
  obj = sobj_list[[3]],
  group_filter = "Neonate RSV infected (NO IFN)",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_2_plots <- plot_spatial_and_bar(
  obj = sobj_list[[4]],
  group_filter = "Neonate RSV infected (NO IFN)",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_3_plots <- plot_spatial_and_bar(
  obj = sobj_list[[3]],
  group_filter = "Neonate IFN and RSV infected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)

Obj_4_plots <- plot_spatial_and_bar(
  obj = sobj_list[[4]],
  group_filter = "Neonate IFN and RSV infected",
  group_by = "niches_5",
  border_color = "white",
  fov = "FOV",
  border_size = 0.1,
  cols = niche_cols,
  size = 0.8,
  axes = FALSE,
  bar_var = "major_celltype",
  bar_group_by = "niches_5",
  bar_colors = celltype_cols,
  bar_main = ""
)


Obj_1_plots$spatial_plot
Obj_1_plots$bar_plot

Obj_2_plots$spatial_plot
Obj_2_plots$bar_plot


Obj_3_plots$spatial_plot
Obj_3_plots$bar_plot

Obj_4_plots$spatial_plot
Obj_4_plots$bar_plot



## Adjusting niche ... didnt do it. It seemed all right 





### All the objects into one for "single cell vis"
All_Sobj <- qread('/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_merged_Sobj.qs')


## Plot Dotplot 
library(presto)
All_Sobj@assays$RNA <- All_Sobj@assays$Nanostring
DefaultAssay(All_Sobj) <- "RNA"
All_Sobj <- All_Sobj %>% NormalizeData()

markers <- wilcoxauc(All_Sobj, 'major_celltype')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 5) %>% 
  ungroup() -> top_markers
table(top_markers$group)

library(ggplot2)
library(viridis)
Idents(All_Sobj) <- "major_celltype"
DotPlot(All_Sobj, 
        #features = unique(top_markers$feature)
        features =  unique(c(top_markers$feature))
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Cell Typer Markers")


# # Filter DEGs based only in FC
# markers %>%
#   group_by(group) %>%
#   dplyr::filter(logFC > 0 & padj <= 0.05) %>%
#   # slice_head(n = 10) %>%
#   dplyr::top_n(wt = logFC, n = 20) %>% 
#   ungroup() -> top_markers
# table(top_markers$group)
# 
# library(ggplot2)
# library(viridis)
# Idents(All_Sobj) <- "major_celltype"
# DotPlot(All_Sobj, 
#         #features = unique(top_markers$feature)
#         features =  unique(c(top_markers[top_markers$group %in% "Lymphocyte", ]$feature))
# ) +
#   geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
#   scale_colour_viridis(option="magma") +
#   guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Lymphocyte top markers ONLY")


## Plot Cell Prop 
library(dittoSeq)
dittoBarPlot(All_Sobj, 
             "major_celltype", 
             group.by = "Group", 
             color.panel = celltype_cols,
             main = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) + ggtitle("")




library(dittoSeq)
library(ggplot2)

# Generate the bar plot and extract data
result <- dittoBarPlot(
  All_Sobj, 
  "major_celltype", 
  group.by = "Group", 
  color.panel = celltype_cols,
  main = "",
  data.out = TRUE
)

# Extract plot and data
p <- result$p
plot_data <- result$data

# Add percentage labels inside the bars
p + 
  geom_text(
    data = plot_data,
    aes(
      x = grouping, 
      y = percent, 
      label = paste0(round(percent, 2), ""), 
      group = label
    ),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) +
  ggtitle("")





library(dittoSeq)
dittoBarPlot(subset(All_Sobj, Group %in% 
                      c(comparisons$`Comparison 3`$group_1,
                        comparisons$`Comparison 3`$group_2)), 
             "niches_5", 
             group.by = "TMA", 
             split.by = "Group",
             color.panel = niche_cols,
             main = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) + ggtitle("")


library(dittoSeq)
dittoBarPlot(subset(All_Sobj, Group %in% 
                      c(comparisons$`Comparison 3`$group_1,
                        comparisons$`Comparison 3`$group_2)), 
             "major_celltype", 
             group.by = "TMA", 
             split.by = "Group",
             color.panel = celltype_cols,
             main = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) + ggtitle("")


# ## Niches already matched 
# qsave(sobj_list, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj_list_celltype_niche.qs")





