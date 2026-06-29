#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
})

# =========================================================
# Config
# =========================================================
input_dir <- "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects"
pattern <- "Sobj_list_TMA_ID_fixed.qs"

output_file <- file.path(
  "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/data/Processed",
  "merged_Sobj.qs"
)

log_file <- file.path(
  dirname(output_file),
  "merge_seurat_obj.log"
)

gene_list_file <- NULL

project_name <- "CosMx_1k_RSVinfe"
assay_name <- "Nanostring"

# Name of metadata column to store original cell IDs before merge/RenameCells
cell_id_column <- "cell_id_prior_merge"

# =========================================================
# Setup
# =========================================================
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Logging
# =========================================================
log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

cat("", file = log_file)
log_message("Merge pipeline started")
log_message("input_dir =", input_dir)
log_message("pattern =", pattern)
log_message("output_file =", output_file)
log_message("log_file =", log_file)
log_message("gene_list_file =", ifelse(is.null(gene_list_file), "NULL", gene_list_file))
log_message("project_name =", project_name)
log_message("assay_name =", assay_name)
log_message("cell_id_column =", cell_id_column)

# =========================================================
# Allowed genes logic
# =========================================================
allowed_genes <- NULL

if (!is.null(gene_list_file)) {
  if (!file.exists(gene_list_file)) {
    stop("Gene list file not found: ", gene_list_file)
  }
  
  allowed_genes <- readLines(gene_list_file, warn = FALSE)
  allowed_genes <- unique(allowed_genes)
  allowed_genes <- allowed_genes[nzchar(allowed_genes)]
  
  log_message("Loaded", length(allowed_genes), "allowed genes from gene list")
} else {
  log_message("gene_list_file is NULL -> all genes from each object will be used")
}

# =========================================================
# Helpers
# =========================================================
is_seurat_object <- function(x) {
  inherits(x, "Seurat")
}

get_base_name <- function(file_path) {
  sub("\\.qs$", "", basename(file_path))
}

make_unique_names <- function(x) {
  make.unique(as.character(x), sep = "_")
}

# =========================================================
# Read source(s)
# Supports:
#   1) multiple .qs files, each with one Seurat object
#   2) one .qs file containing a list of Seurat objects
# =========================================================
read_input_objects <- function(input_dir, pattern) {
  qs_files <- sort(list.files(
    input_dir,
    pattern = pattern,
    full.names = TRUE
  ))
  
  log_message("Found", length(qs_files), "matching .qs file(s)")
  
  if (length(qs_files) == 0) {
    stop("No matching .qs files found")
  }
  
  # Case 1: multiple files -> assume each file is one Seurat object
  if (length(qs_files) > 1) {
    log_message("Input mode: multiple .qs files, each expected to contain one Seurat object")
    
    object_list <- vector("list", length(qs_files))
    object_names <- character(length(qs_files))
    
    for (i in seq_along(qs_files)) {
      fp <- qs_files[i]
      nm <- get_base_name(fp)
      
      log_message("Reading file:", fp)
      x <- qread(fp)
      
      if (!is_seurat_object(x)) {
        stop("File does not contain a Seurat object: ", fp)
      }
      
      object_list[[i]] <- x
      object_names[i] <- nm
    }
    
    names(object_list) <- make_unique_names(object_names)
    return(object_list)
  }
  
  # Case 2: one file -> can be Seurat object OR list of Seurat objects
  fp <- qs_files[1]
  log_message("Reading single .qs file:", fp)
  x <- qread(fp)
  
  # Single Seurat object
  if (is_seurat_object(x)) {
    nm <- get_base_name(fp)
    log_message("Input mode: single .qs containing one Seurat object")
    out <- list(x)
    names(out) <- nm
    return(out)
  }
  
  # List of Seurat objects
  if (is.list(x)) {
    log_message("Input mode: single .qs containing a list")
    
    keep <- vapply(x, is_seurat_object, logical(1))
    if (!any(keep)) {
      stop("Single .qs file contains a list, but no Seurat objects were found: ", fp)
    }
    
    x <- x[keep]
    
    obj_names <- names(x)
    if (is.null(obj_names)) {
      obj_names <- paste0(get_base_name(fp), "_obj", seq_along(x))
    } else {
      obj_names[obj_names == ""] <- paste0("obj", which(obj_names == ""))
      obj_names <- paste0(get_base_name(fp), "_", obj_names)
    }
    
    names(x) <- make_unique_names(obj_names)
    
    log_message("Recovered", length(x), "Seurat object(s) from Seurat-object list")
    return(x)
  }
  
  stop("Unsupported content in .qs file: ", fp)
}

# =========================================================
# Prepare object
# =========================================================
prepare_obj <- function(
    obj,
    sample_name,
    allowed_genes = NULL,
    assay_name = "RNA",
    cell_id_column = "cell_id_prior_merge"
) {
  DefaultAssay(obj) <- assay_name
  
  genes_before <- nrow(obj)
  cells_before <- ncol(obj)
  
  if (!is.null(allowed_genes)) {
    keep_genes <- intersect(rownames(obj), allowed_genes)
    
    if (length(keep_genes) == 0) {
      stop("No overlapping genes found for sample/object: ", sample_name)
    }
    
    obj <- subset(obj, features = keep_genes)
    
    log_message(
      "Prepared object:", sample_name,
      "| genes_before =", genes_before,
      "| genes_after =", nrow(obj),
      "| cells =", ncol(obj)
    )
  } else {
    log_message(
      "Prepared object without gene filtering:", sample_name,
      "| genes =", genes_before,
      "| cells =", cells_before
    )
  }
  
  # Preserve original cell IDs before RenameCells()
  if (cell_id_column %in% colnames(obj@meta.data)) {
    log_message(
      "WARNING | Column", cell_id_column,
      "already exists in", sample_name,
      "- it will be overwritten"
    )
  }
  
  obj[[cell_id_column]] <- rownames(obj@meta.data)
  
  # Add sample metadata and rename cells to keep merged object unique
  obj$sample <- sample_name
  obj <- RenameCells(obj, add.cell.id = sample_name)
  
  log_message(
    "Stored original cell IDs in metadata column:", cell_id_column,
    "| sample =", sample_name
  )
  
  obj
}

# =========================================================
# Load all input objects
# =========================================================
obj_list <- read_input_objects(input_dir, pattern)


sample_names <- names(obj_list)

log_message("Total Seurat objects to merge:", length(obj_list))
log_message("Object/sample names:", paste(sample_names, collapse = ", "))

if (length(obj_list) == 0) {
  stop("No Seurat objects available to merge")
}

# =========================================================
# Prepare all objects
# =========================================================
for (i in seq_along(obj_list)) {
  log_message("Preparing object", i, "of", length(obj_list), ":", sample_names[i])
  obj_list[[i]] <- prepare_obj(
    obj = obj_list[[i]],
    sample_name = sample_names[i],
    allowed_genes = allowed_genes,
    assay_name = assay_name,
    cell_id_column = cell_id_column
  )
  gc()
}

# =========================================================
# Single object case
# =========================================================
if (length(obj_list) == 1) {
  log_message("Only one Seurat object available; saving directly")
  obj <- obj_list[[1]]
  obj <- JoinLayers(obj)
  
  qsave(obj, output_file, preset = "high")
  log_message("Saved single object to:", output_file)
  log_message("Cells =", ncol(obj), "| Genes =", nrow(obj))
  log_message("Merge pipeline finished")
  quit(save = "no")
}

# =========================================================
# Merge first two
# =========================================================
log_message("Merging first two objects:", sample_names[1], "+", sample_names[2])
obj <- merge(
  x = obj_list[[1]],
  y = obj_list[[2]],
  project = project_name
)

log_message(
  "Merged object stats after first merge",
  "| cells =", ncol(obj),
  "| genes =", nrow(obj)
)

log_message("Joining layers after first merge")
obj <- JoinLayers(obj)

log_message(
  "Post-JoinLayers stats after first merge",
  "| cells =", ncol(obj),
  "| genes =", nrow(obj)
)

qsave(obj, output_file, preset = "high")
log_message("Saved merged object:", output_file)

rm(obj)
gc()

# =========================================================
# Iterative merge
# =========================================================
if (length(obj_list) > 2) {
  for (i in 3:length(obj_list)) {
    next_sample <- sample_names[i]
    
    log_message("Reading current merged object from disk:", output_file)
    obj <- qread(output_file)
    
    log_message(
      "Current merged object loaded",
      "| cells =", ncol(obj),
      "| genes =", nrow(obj)
    )
    
    log_message("Merging current merged object with:", next_sample)
    obj <- merge(
      x = obj,
      y = obj_list[[i]],
      project = project_name
    )
    
    log_message(
      "Merged object stats after adding", next_sample,
      "| cells =", ncol(obj),
      "| genes =", nrow(obj)
    )
    
    log_message("Joining layers after adding", next_sample)
    obj <- JoinLayers(obj)
    
    log_message(
      "Post-JoinLayers stats after adding", next_sample,
      "| cells =", ncol(obj),
      "| genes =", nrow(obj)
    )
    
    qsave(obj, output_file, preset = "high")
    log_message("Overwrote merged object:", output_file)
    
    rm(obj)
    gc()
  }
}

# =========================================================
# Final check
# =========================================================
log_message("Reading final merged object for verification:", output_file)
final_obj <- qread(output_file)

log_message(
  "Final merged object verified",
  "| cells =", ncol(final_obj),
  "| genes =", nrow(final_obj)
)

log_message("Merge pipeline finished")