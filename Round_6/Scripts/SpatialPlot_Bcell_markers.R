library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(shiny)
library(plotly)
library(ggplot2)
library(data.table)
library(FNN)
library(qs)
library(presto)
library(UpSetR)
library(pheatmap)
library(gridExtra)
library(shiny)
library(ggplot2)
library(pheatmap)



## Load data ---------
obj_list <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Objects/Sobj_list_TMA_ID_fixed.qs")

obj_list$`NA` <- NULL
obj_list$`NA.1` <- NULL



## Spatial Plot Bcell markers ---------
obj_list

# Make sure required packages are loaded
library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)

# Extract dimensions for each object
dims_df <- purrr::imap_dfr(obj_list, function(obj, name) {
  d <- dim(obj)
  tibble::tibble(
    object_name = name,
    n_features = d[1],
    n_cells = d[2]
  )
})

# Print table
print(dims_df %>% data.frame())


obj_list[["Ap16-3894 A1 3"]]@assays$Nanostring$counts %>% rownames()

# chatBOT strong Bcell markers based on this CosMx 1k panel: "Cd19", #Cd19, Ms4a1, Cd79a, Pax5, Cd22, Igkc

Seurat::ImageFeaturePlot(obj_list[["Ap16-3894 A1 3"]],
                         features =  "Cd19", #Cd19, Ms4a1, Cd79a, Pax5, Cd22, Igkc
                         max.cutoff = c("q99"),
                         fov = "FOV") &
  scale_fill_viridis() &
  #coord_flip() &
  ggtitle("") &
  theme(plot.title = element_text(size = 8, hjust = 0))




library(Seurat)
library(ggplot2)
library(viridis)
library(purrr)

features_vec <- c("Cd19", "Ms4a1", "Cd79a", "Pax5", "Cd22", "Igkc")

plot_list <- purrr::map(
  features_vec,
  function(feat) {

    Seurat::ImageFeaturePlot(obj_list[["Ap16-3894 A1 3"]],
                             features =  feat,
                             max.cutoff = c("q99"),
                             fov = "FOV") &
      scale_fill_viridis() &
      #coord_flip() &
      ggtitle("") &
      theme(plot.title = element_text(size = 8, hjust = 0))
    
  }
)

names(plot_list) <- features_vec
plot_list




library(Seurat)
library(purrr)
library(ggplot2)
library(viridis)
library(patchwork)



plot_spatial_by_group_tma_v2 <- function(
    obj_list,
    features = "Cd19",
    fov = "FOV",
    max_cutoff = "q99",
    export_pdf = FALSE,
    pdf_file = "ImageFeaturePlot_by_Group_TMA.pdf",
    pdf_width = 5,
    pdf_height = 5
) {
  
  # helper to check metadata presence
  check_meta <- function(obj, col) {
    if (!col %in% colnames(obj@meta.data)) {
      stop(sprintf("Object is missing metadata column: '%s'", col))
    }
  }
  
  plots_nested <- purrr::imap(obj_list, function(obj, obj_name) {
    
    # ensure expected metadata exists
    check_meta(obj, "Group")
    check_meta(obj, "TMA")
    
    groups <- sort(unique(as.character(obj$Group)))
    out_group <- purrr::map(groups, function(g) {
      
      # TMAs present for this group
      tmas <- sort(unique(as.character(obj$TMA[obj$Group == g])))
      
      out_tma <- purrr::map(tmas, function(t) {
        
        # select cells explicitly (avoids Seurat::subset)
        sel_cells <- colnames(obj)[which(obj$Group == g & obj$TMA == t)]
        
        if (length(sel_cells) == 0) {
          warning(sprintf("No cells for Object '%s' Group '%s' TMA '%s' — returning NULL plot", obj_name, g, t))
          return(NULL)
        }
        
        obj_sub <- obj[, sel_cells, drop = FALSE]  # safe explicit subsetting
        
        
        mat <- obj_sub@assays$Nanostring$counts
        if (feat %in% rownames(mat)) {
          expr_vec <- mat[feat, ]
          cells_use <- colnames(obj_sub)[which(expr_vec >= 1)]
        } else {
          cells_use <- colnames(obj_sub)
        }
        
        # mat <- obj_sub@assays$Nanostring$counts
        # 
        # if (feat %in% rownames(mat)) {
        #   cells_keep <- colnames(obj_sub)[mat[feat, ] >= 1]
        #   obj_sub <- obj_sub[, cells_keep, drop = FALSE]
        # }
        # 
        
        feat_plots <- purrr::map(features, function(feat) {
          
          # your exact plot code, with feat variable
          
          Seurat::ImageFeaturePlot(
            obj_sub,
            #cells = cells_use, 
            features = feat, 
            max.cutoff = c(max_cutoff), 
            border.size = NA,
            #border.color = NA,
            #ize = 0.7,
            fov = fov,
            dark.background = TRUE,
            axes = FALSE
          ) &
            scale_fill_viridis(begin = 0) & # being > 0 to not start from purple 
            ggtitle("") &
            theme(
              panel.grid = element_blank(),     # remove grid
              panel.background = element_rect(fill = "black"),
              plot.background = element_rect(fill = "black"),
              legend.background = element_rect(fill = "black"),
              legend.key = element_rect(fill = "black"),
              plot.title = element_text(size = 8, hjust = 0, color = "white"),
              legend.text = element_text(color = "white"),
              legend.title = element_text(color = "white")
            )
        })
        
        names(feat_plots) <- features
        
        # combine if multiple features
        if (length(features) == 1) {
          p <- feat_plots[[1]]
        } else {
          p <- patchwork::wrap_plots(feat_plots) +
            patchwork::plot_annotation(title = NULL)
        }
        
        # label the page with object/group/tma
        p + patchwork::plot_annotation(
          title = paste0(obj_name, " | Group: ", g, " | TMA: ", t)
        )
      })
      
      names(out_tma) <- tmas
      out_tma
    })
    
    names(out_group) <- groups
    out_group
  })
  
  # export to single multi-page PDF if requested
  if (isTRUE(export_pdf)) {
    grDevices::pdf(pdf_file, width = pdf_width, height = pdf_height, onefile = TRUE)
    
    purrr::iwalk(plots_nested, function(by_group, obj_name) {
      purrr::iwalk(by_group, function(by_tma, g) {
        purrr::iwalk(by_tma, function(p, t) {
          if (!is.null(p)) print(p)
        })
      })
    })
    
    grDevices::dev.off()
  }
  
  invisible(plots_nested)
}

plots <- plot_spatial_by_group_tma_v2(
  obj_list[],
  features = c("Cd19", "Ms4a1", "Cd79a", "Pax5", "Cd22", "Igkc"),
  pdf_file = "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_6/Plots/Bcell_marker_by_Group_TMA.pdf",
  pdf_width = 10,
  pdf_height = 10,
  export_pdf = T
 
)

plots[1]

# Note: I'm able to detect a couple of TMAs presting the "Bcell clamp" very marked on "Igkc" expression
# They're 
# Ap16−3894 A1 3 | Group: Adult reinfected | TMA: Ap16−3894 3
# Ap17−0332 A1 1.1 | Group: Adult Control | TMA: Ap17−0332 A1 1
# Ap16−2655 A2 2 | Group: Neonate (NO IFN) reinfected | TMA: Ap16−2655 A2 1
# Ap16−2657 A2 1 | Group: Neonate IFN and RSV reinfected | TMA: Ap16−2657 A2 2
# and
# Ap17−0333 A1 3 | Group: Adult Control | TMA: Ap17−0333 3


