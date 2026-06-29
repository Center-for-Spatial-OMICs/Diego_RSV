#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(Matrix)
})

# =========================
# INPUT / OUTPUT / CONFIG
# =========================
input_file  <- "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/data/Processed/merged_Sobj.qs"
output_file <- "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/data/Processed/merged_Sobj_UMAP.qs"

# assay from config
assay_name <- "Nanostring"

# QC thresholds from config
min_features <- 15
min_counts   <- 30

cluster_res_02 <- 0.2
cluster_res_05 <- 0.5

dimensions <- 1:30
# =========================
# LOG FILE
# same basename as current script, saved in output directory
# =========================
args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_full, value = TRUE)

script_name <- if (length(script_arg) > 0) {
  normalizePath(sub("^--file=", "", script_arg[1]))
} else {
  "seurat_workflow.R"
}

log_file <- file.path(
  dirname(output_file),
  paste0(tools::file_path_sans_ext(basename(script_name)), ".log")
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message", append = TRUE)

on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

cat("========================================\n")
cat("Pipeline started:", as.character(Sys.time()), "\n")
cat("Input file :", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("Log file   :", log_file, "\n")
cat("Assay      :", assay_name, "\n")
cat("========================================\n\n")

# =========================
# LOAD OBJECT
# =========================
cat("[1] Reading object...\n")
obj <- qread(input_file)

# =========================
# QC METRICS
# get nfeature, ncount from assay counts
# filter out cells with low QC values:
# nFeature < 30, nCount < 60
# =========================
cat("[2] Computing QC metrics from assay counts...\n")
DefaultAssay(obj) <- assay_name

counts_mat <- GetAssayData(obj, assay = assay_name, slot = "counts")

obj$nCount   <- Matrix::colSums(counts_mat)
obj$nFeature <- Matrix::colSums(counts_mat > 0)

cat("Cells before QC:", ncol(obj), "\n")
obj <- subset(
  obj,
  subset = nFeature >= min_features & nCount >= min_counts
)
cat("Cells after QC :", ncol(obj), "\n\n")

# =========================
# NORMALIZATION / FEATURE SELECTION / SCALING
# =========================
cat("[3] NormalizeData...\n")
obj <- NormalizeData(obj)

cat("[4] FindVariableFeatures...\n")
obj <- FindVariableFeatures(
  obj,
  nfeatures = nrow(obj@assays[[assay_name]]$counts)
)

cat("[5] ScaleData...\n")
obj <- ScaleData(obj)

# =========================
# PCA
# =========================
cat("[6] RunPCA...\n")
obj <- Seurat::RunPCA(obj)
qsave(obj, output_file)
cat("Saved after RunPCA ->", output_file, "\n\n")

# =========================
# NEIGHBORS
# =========================
cat("[7] FindNeighbors...\n")
obj <- FindNeighbors(obj, dims = dimensions)
qsave(obj, output_file)
cat("Saved after FindNeighbors ->", output_file, "\n\n")

# =========================
# UMAP
# =========================
cat("[8] RunUMAP...\n")
obj <- RunUMAP(obj, dims = dimensions, return.model = FALSE)
qsave(obj, output_file)
cat("Saved after RunUMAP ->", output_file, "\n\n")


# =========================
# CLUSTERING
# =========================
cat("[9] FindClusters...\n")
obj <- FindClusters(obj, resolution = c(cluster_res_02, cluster_res_05))
qsave(obj, output_file)
cat("Saved after FindClusters ->", output_file, "\n\n")

cat("========================================\n")
cat("Pipeline finished:", as.character(Sys.time()), "\n")
cat("Final object saved to:", output_file, "\n")
cat("========================================\n")


