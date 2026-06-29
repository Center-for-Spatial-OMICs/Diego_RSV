## ====== Date: Jun 18, 2025  ===== 

# ========================================================
## Directions: 
# Re-annotate the dataset for major cell types
# Find lymphocytes
# Then separate Bcells from Tcells 
# Go for lym. minor cell types 
# Get DEGs across the three groups

## Workflow 

# ========================================================


### Load and define functions ---------
library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(qs)


call_insitu_cluster <- function(Sobj, n_of_cluster, assay, return_flightpath = FALSE) {
  # Check if assay exists
  if (!assay %in% names(Sobj@assays)) {
    message(sprintf("Assay '%s' not found in the Seurat object. Returning input object unchanged.", assay))
    if (return_flightpath) {
      return(list(Sobj = Sobj, flightpath = NULL))
    } else {
      return(Sobj)
    }
  }
  
  # Store flightpaths if requested
  flightpaths <- list()
  
  for (n in n_of_cluster) {
    # Run insitutypeML
    unsup <- insitutype(
      x = as.matrix(t(as.matrix(Sobj@assays[[assay]]$counts))),
      neg = Sobj@meta.data$nCount_negprobes,
      reference_profiles = NULL,
      bg = NULL,
      n_clusts = n,
      n_phase1 = 200,
      n_phase2 = 500,
      n_phase3 = 2000,
      n_starts = 1,
      max_iters = 5
    )
    
    # Create new column for this cluster number, default 'other'
    cluster_col <- paste0("Inst_Cluster_n", n)
    Sobj@meta.data[[cluster_col]] <- "other"
    if(identical(rownames(Sobj@meta.data), names(unsup$clust))) {
      Sobj@meta.data[[cluster_col]] <- unsup$clust
    } else {
      print("We can't add cluster label into Sobj@meta.data that way, bro...")
    }
    
    # Plot flightpath if requested
    if (return_flightpath) {
      cols = RColorBrewer::brewer.pal(10, 'Paired')
      cols <- cols[seq_along(unique(unsup$clust))]
      names(cols) <- unique(unsup$clust)
      fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup, col = cols[unsup$clust])
      print(fp)
      flightpaths[[cluster_col]] <- fp
    }
  }
  
  if (return_flightpath) {
    return(list(Sobj = Sobj, flightpaths = flightpaths))
  } else {
    return(Sobj)
  }
}


### Load and define variables  ---------
Sobj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj_list_celltype_niche.qs")

All_Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_5/Objects/scObj_processed.qs")

comparisons <- list(
  "Comparison" = list(
    group_1 = "Neonate RSV infected (NO IFN)",
    group_2 = "Neonate IFN and RSV infected",
    group_3 = "Adult RSV infected"
  )
)

load("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/Lung_MCA.RData")
metadata
cellGroups
profile_matrix[1:10, 1:10] # row: mouse genes; col: cell type; content: cell type "fraction" by gene
ref_sc_data <- profile_matrix
head(ref_sc_data)
colnames(ref_sc_data)

major_celltypes <- c(
  "Alveolar", "AT1", "AT2", "B", "Basophil", "Ciliated", "Clara",
  "Dendritic", "Endothelial", "Eosinophil", "IgA", "Interstitial",
  "Monocyte", "Neutrophil", "NK", "Nuocyte", "Plasmacytoid",
  "Stromal", "T"
)

major_celltypes_short <- c(
  "Alveolar", 
  "Alveolar.bipotent.progenitor",  # so we can have a positive control cause it should be found mainly in neonate rather than adults 
  "Ciliated", 
  "Clara-nonCiliated",
  "Dendritic", 
  "Endothelial", 
  "Interstitial", 
  "Stromal" 
)

# Lymphocytes markers from Diego 
lym_markers_Diego <- list("Lymphocytes" = list(
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
  )))


### Unsup. Insitutype clustering  ---------
# set n of clusters to be == expected cell types
All_Sobj <- call_insitu_cluster(Sobj = All_Sobj, 
                    n_of_cluster = length(major_celltypes_short), 
                    assay = "RNA", 
                    return_flightpath = TRUE)
flightpaths_plot <- All_Sobj$flightpaths
All_Sobj <- All_Sobj$Sobj


adu1 <- call_insitu_cluster(Sobj = subset(All_Sobj, Slide %in% "Adult_1"), 
                                n_of_cluster = length(major_celltypes_short), 
                                assay = "RNA", 
                                return_flightpath = TRUE)


adu2 <- call_insitu_cluster(Sobj = subset(All_Sobj, Slide %in% "Adult_2"), 
                                 n_of_cluster = length(major_celltypes_short), 
                                 assay = "RNA", 
                                 return_flightpath = TRUE)

neo1 <- call_insitu_cluster(Sobj = subset(All_Sobj, Slide %in% "Neonate_1"), 
                            n_of_cluster = length(major_celltypes_short), 
                            assay = "RNA", 
                            return_flightpath = TRUE)

neo2 <- call_insitu_cluster(Sobj = subset(All_Sobj, Slide %in% "Neonate_2"), 
                            n_of_cluster = length(major_celltypes_short), 
                            assay = "RNA", 
                            return_flightpath = TRUE)

# Get cluster markers 
library(presto)
markers <- wilcoxauc(adu1$Sobj, 'Inst_Cluster_n8')
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 10) %>% 
  ungroup() %>% 
  data.frame() -> top_markers

library(ggplot2)
library(viridis)
Idents(adu1$Sobj) <- "Inst_Cluster_n8"
DotPlot(adu1$Sobj, 
        features =  unique(c(top_markers$feature))

) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16)) + ggtitle("Inst_Cluster_n8 Adult 1 - cluster markers")

library(ggplot2)
library(viridis)
Idents(adu1$Sobj) <- "Inst_Cluster_n8"
DotPlot(adu1$Sobj, 
        features =  unique(c(unlist(lym_markers_Diego) %>% as.vector()))

) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16)) + ggtitle("Inst_Cluster_n8 Adult 1 - lymphocyte markers")

Sobj <- adu1$Sobj
DefaultAssay(Sobj) <- "RNA"
Sobj <- NormalizeData(Sobj)
Sobj <- FindVariableFeatures(Sobj, nfeatures = length(rownames(Sobj@assays[["RNA"]]$counts)) )
Sobj <- ScaleData(Sobj)
Sobj <- RunPCA(Sobj)
Sobj <- FindNeighbors(Sobj, dims = 1:15)
#Sobj <- FindClusters(Sobj, resolution = c(0.2, 0.5)) 
Sobj <- RunUMAP(Sobj, dims = 1:15, return.model = T)
adu1$Sobj <- Sobj

Idents(adu1$Sobj) <- "Inst_Cluster_n8"
DimPlot(adu1$Sobj)

FeaturePlot(adu1$Sobj, 
            features = c(unlist(lym_markers_Diego) %>% as.vector()),
            order = TRUE
            )



### Testing Felipe's code to check on cell annotation quality -----
# Load ref single cell data
ref_sc_data
obj_list <- SplitObject(object = All_Sobj, split.by = 'Slide')


library(Seurat)
library(dplyr)
library(purrr)
# install.packages('ggpmisc')
library(ggpmisc)

obj_list <- lapply(obj_list, function(x) { x <- NormalizeData(x); return(x)})


# Loop over each Seurat object

plot_list <- list()  
for (i in seq_along(obj_list)) {
  obj <- obj_list[[i]]
  obj_name <- names(obj_list)[i] %||% paste0("Object_", i)
  
  # Make sure insitutype_cluster_Superv_celltypes exists and is a factor
  if (!"insitutype_cluster_Superv_celltypes" %in% colnames(obj@meta.data)) {
    message(paste("Skipping", obj_name, "- insitutype_cluster_Superv_celltypes not found"))
    next
  }
  
  obj$insitutype_cluster_Superv_celltypes <- as.factor(obj$insitutype_cluster_Superv_celltypes)
  
  # Get average expression
  avg_expr <- AverageExpression(obj, group.by = "insitutype_cluster_Superv_celltypes", assays = "RNA", slot = "data")$RNA
  
  # Intersect genes and cell types
  common_celltypes <- intersect(colnames(avg_expr), colnames(ref_sc_data))
  common_genes <- intersect(rownames(avg_expr), rownames(ref_sc_data))
  
  if (length(common_celltypes) == 0 || length(common_genes) == 0) {
    message(paste("Skipping", obj_name, "- no overlap in genes or cell types"))
    next
  }
  
  expr_obj <- avg_expr[common_genes, common_celltypes]
  expr_ref <- as.matrix(ref_sc_data[common_genes, common_celltypes])
  
  # Melt into long format for plotting
  plot_df <- data.frame(
    Gene = rep(common_genes, times = length(common_celltypes)),
    Celltype = rep(common_celltypes, each = length(common_genes)),
    Expr_obj = as.vector(expr_obj),
    Expr_ref = as.vector(expr_ref)
  )
  
  # Plot
  p <- ggplot(plot_df, aes(x = log10(Expr_ref + 0.1), y = log10(Expr_obj + 0.1))) +
    geom_point(alpha = 0.5) +
    facet_wrap(~ Celltype, scales = "free") +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3) +
    theme_bw() +
    labs(
      title = paste("Gene Expression Correlation with Reference -", obj_name),
      x = "log10(Expression in ref_sc_data + 0.1)",
      y = paste("log10(Expression in", obj_name, "+ 0.1)")
    )
  
  plot_list[[obj_name]] <- p  
}

plot_list$Adult_1
plot_list$Adult_2
plot_list$Neonate_1
plot_list$Neonate_2




### Current cell type markers and Inst_Cluster_n8 markers on the sc object
library(presto)
markers <- wilcoxauc(All_Sobj, 'major_celltype') 
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 5) %>% 
  ungroup() %>% 
  data.frame() -> top_markers

library(ggplot2)
library(viridis)
Idents(All_Sobj) <- "major_celltype"
DotPlot(All_Sobj, 
        features =  unique(c(top_markers$feature))
        
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16)) + ggtitle("Current cell type markers")




library(presto)
markers <- wilcoxauc(All_Sobj, 'Inst_Cluster_n8') 
markers %>% head()
dim(markers) 
# Filter DEGs based only in FC
markers %>%
  group_by(group) %>%
  dplyr::filter(logFC > 0 & padj <= 0.05) %>%
  # slice_head(n = 10) %>%
  dplyr::top_n(wt = logFC, n = 5) %>% 
  ungroup() %>% 
  data.frame() -> top_markers

library(ggplot2)
library(viridis)
Idents(All_Sobj) <- "Inst_Cluster_n8"
DotPlot(All_Sobj, 
        features =  unique(c(top_markers$feature))
        
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16)) + ggtitle("Inst_Cluster_n8 markers")



