library(qs)
library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(zellkonverter)

input_qs <- "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_all467158cells_clustered.qs"
output_h5ad <- "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/SeuObj_Nanostring.h5ad"

SeuObj <- qread(input_qs)
DefaultAssay(SeuObj) <- "Nanostring"

counts <- LayerData(SeuObj, assay = "Nanostring", layer = "counts")
data <- LayerData(SeuObj, assay = "Nanostring", layer = "data")
scale_data <- LayerData(SeuObj, assay = "Nanostring", layer = "scale.data")

sce <- SingleCellExperiment(
  assays = list(
    counts = counts,
    logcounts = data,
    scaledata = scale_data
  ),
  colData = DataFrame(SeuObj@meta.data, row.names = rownames(SeuObj@meta.data)),
  rowData = DataFrame(gene = rownames(counts), row.names = rownames(counts))
)

writeH5AD(
  sce,
  file = output_h5ad,
  X_name = "counts"
)

cat("Wrote:", output_h5ad, "\n")
cat("Cells:", ncol(sce), " Genes:", nrow(sce), "\n")
