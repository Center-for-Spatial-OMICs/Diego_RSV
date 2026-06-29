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

source("/mnt/scratch2/Maycon/Utils/R_codabase/Utils.R")

### Load data ------
SeuObj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/Sobj_all467158cells_clustered.qs")

# # columns 
# SeuObj$TMA_4 # TMA ID 
# SeuObj$Group # infection groups 
# SeuObj$Group_2 # adult or neoante labels 
# SeuObj$cluster_regular_res03 # cluster that originated "major_celltype"
# SeuObj$major_celltype # major cell types of Diego's interest (not from insitutype)
# SeuObj$cell_based_final # Tcd4, Tcd8, Bcell, pDC, cDC - no minor for epithelail for lack of markers in the 1k panel 
# SeuObj$cell_id_prior_merge # whenever we need to get those cell anno back to the obj_list with polygons/transcripts


### Number of cells within each comparison group --------
library(dplyr)
library(ggplot2)
library(patchwork)

meta <- SeuObj@meta.data

# comparisons
prim_inf_comp <- list(
  c("Neonate RSV infected (NO IFN)", "Neonate IFN and RSV infected"),
  c("Neonate RSV infected (NO IFN)", "Adult RSV infected"),
  c("Neonate IFN and RSV infected", "Adult RSV infected")
)

reinf_comp <- list(
  c("Neonate (NO IFN) reinfected", "Neonate IFN and RSV reinfected"),
  c("Neonate (NO IFN) reinfected", "Adult reinfected"),
  c("Neonate IFN and RSV reinfected", "Adult reinfected")
)

all_comp <- c(prim_inf_comp, reinf_comp)

names(all_comp) <- c(
  "Primary: Neo RSV no IFN vs Neo IFN+RSV",
  "Primary: Neo RSV no IFN vs Adult RSV",
  "Primary: Neo IFN+RSV vs Adult RSV",
  "Reinfected: Neo no IFN vs Neo IFN+RSV",
  "Reinfected: Neo no IFN vs Adult",
  "Reinfected: Neo IFN+RSV vs Adult"
)

# build comparison dataframe
plot_df <- lapply(names(all_comp), function(comp_name) {
  
  groups_use <- all_comp[[comp_name]]
  
  meta %>%
    dplyr::filter(
      Group %in% groups_use,
      !cell_based_final %in% c(
        "other",
        "Ambiguous_overlap_across_diff_cell_based"
      )
    ) %>%
    dplyr::count(Group, cell_based_final, name = "n_cells") %>%
    dplyr::mutate(comparison = comp_name)
  
}) %>%
  dplyr::bind_rows()

# barplot
ggplot(plot_df, aes(x = cell_based_final, y = n_cells, fill = Group)) +
  geom_col(position = "dodge") +
  facet_wrap(~ comparison, scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 10),
    legend.title = element_blank()
  ) +
  labs(
    x = "cell_based_final",
    y = "Number of cells"
  )



### DEG --------

# # for 1 comparison 
# obj <- subset(SeuObj, cell_based_final %in% "Tcd4")
# plt <- Plot_Volcano_DEG_sc(SeuratObj = obj,
#                            assay = "Nanostring",
#                            Identity = "Group",
#                            pCutoff = 0.05,
#                            FCcutoff = 0.5,
#                            Group_1 = prim_inf_comp[[1]][1],
#                            Group_2 = prim_inf_comp[[1]][2],
#                            logFC_direction = prim_inf_comp[[1]][2],
#                            #markers_to_label = markers,
#                            remove_mt = TRUE,
#                            genes_to_remove_from_vis = c("Hba-a1/2","Hbb"),
#                            title = "Tcd4")
# 
# plt[[1]]
# 
# right_group <- plt[[2]][plt[[2]]$logFC > 0 & plt[[2]]$padj < 0.05, ]$feature
# plt[[2]][plt[[2]]$feature %in% right_group, ]
# left_group <- plt[[2]][plt[[2]]$logFC < 0 & plt[[2]]$padj < 0.05, ]$feature 
# plt[[2]][plt[[2]]$feature %in% left_group, ]
# 
# DEG <- plt[[2]]
# write.csv(DEG[DEG$padj < 0.05, ], "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/DEGs/DEGs_Tcd4_filtered_padj05.csv")


# looping it for all comparisons 
library(Seurat)
library(dplyr)
library(ggplot2)

# =========================
# Settings
# =========================
celltypes_use <- c("Tcd4", "Tcd8", "Bcell")

outdir <- "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/DEGs/Lyn_sub"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prim_inf_comp <- list(
  c("Neonate RSV infected (NO IFN)", "Neonate IFN and RSV infected"),
  c("Neonate RSV infected (NO IFN)", "Adult RSV infected"),
  c("Neonate IFN and RSV infected", "Adult RSV infected")
)

reinf_comp <- list(
  c("Neonate (NO IFN) reinfected", "Neonate IFN and RSV reinfected"),
  c("Neonate (NO IFN) reinfected", "Adult reinfected"),
  c("Neonate IFN and RSV reinfected", "Adult reinfected")
)

all_comp <- c(prim_inf_comp, reinf_comp)

names(all_comp) <- c(
  "primary_NeoRSVnoIFN_vs_NeoIFNRSV",
  "primary_NeoRSVnoIFN_vs_AdultRSV",
  "primary_NeoIFNRSV_vs_AdultRSV",
  "reinf_NeonoIFN_vs_NeoIFNRSV",
  "reinf_NeonoIFN_vs_Adult",
  "reinf_NeoIFNRSV_vs_Adult"
)

# =========================
# DEG loop
# =========================
for (ct in celltypes_use) {
  
  message("\n==============================")
  message("Running cell type: ", ct)
  message("==============================")
  
  obj_ct <- subset(
    SeuObj,
    subset = cell_based_final == ct
  )
  
  for (comp_name in names(all_comp)) {
    
    Group_1 <- all_comp[[comp_name]][1]
    Group_2 <- all_comp[[comp_name]][2]
    
    message("  Comparison: ", comp_name)
    
    group_counts <- table(obj_ct$Group)
    
    if (!(Group_1 %in% names(group_counts)) | !(Group_2 %in% names(group_counts))) {
      message("    Skipped: one group missing")
      next
    }
    
    if (group_counts[Group_1] == 0 | group_counts[Group_2] == 0) {
      message("    Skipped: one group has 0 cells")
      next
    }
    
    plt <- Plot_Volcano_DEG_sc(
      SeuratObj = obj_ct,
      assay = "Nanostring",
      Identity = "Group",
      Group_1 = Group_1,
      Group_2 = Group_2,
      logFC_direction = Group_2,
      pCutoff = 0.05,
      FCcutoff = 0.5,
      remove_mt = TRUE,
      genes_to_remove_from_vis = c("Hba-a1/2", "Hbb"),
      title = paste0(ct, " | ", comp_name)
    )
    
    DEG <- plt[[2]]
    
    deg_bar <- Plot_Top_DEG_Barplot(
      DEG = DEG,
      Group_1 = Group_1,
      Group_2 = Group_2,
      logFC_direction = Group_2,
      pCutoff = 0.05,
      FCcutoff = 0.5,
      top_n = 20,
      title = paste0(ct, " | ", comp_name)
    )
    
    pos_deg <- DEG[
      DEG$logFC >= 0.5 & DEG$padj < 0.05,
    ]$feature
    
    neg_deg <- DEG[
      DEG$logFC <= -0.5 & DEG$padj < 0.05,
    ]$feature
    
    message("    Positive DEG / higher in: ", Group_2, " = ", length(pos_deg), " genes")
    message("    Negative DEG / higher in: ", Group_1, " = ", length(neg_deg), " genes")
    
    base_name <- paste0(
      "DEGs_",
      ct,
      "_",
      comp_name,
      "_",
      gsub("[^A-Za-z0-9]+", "_", Group_2),
      "_direction"
    )
    
    csv_file <- file.path(outdir, paste0(base_name, "_filtered_padj05_logFC05.csv"))
    volcano_pdf_file <- file.path(outdir, paste0(base_name, "_volcano.pdf"))
    barplot_pdf_file <- file.path(outdir, paste0(base_name, "_topDEG_barplot.pdf"))
    
    write.csv(
      DEG[
        DEG$padj < 0.05 &
          abs(DEG$logFC) >= 0.5,
      ],
      csv_file,
      row.names = FALSE
    )
    
    pdf(volcano_pdf_file, width = 14, height = 8)
    print(plt[[1]])
    dev.off()
    
    pdf(barplot_pdf_file, width = 10, height = 8)
    print(deg_bar)
    dev.off()
  }
}


library(Seurat)
library(dplyr)
library(ggplot2)

# =========================
# Settings
# =========================
celltypes_use <- c("cDC (Dendritic)", "Apoe+/C1q+ macrophage/APC (Dendritic)")

outdir <- "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/DEGs/Den_sub"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prim_inf_comp <- list(
  c("Neonate RSV infected (NO IFN)", "Neonate IFN and RSV infected"),
  c("Neonate RSV infected (NO IFN)", "Adult RSV infected"),
  c("Neonate IFN and RSV infected", "Adult RSV infected")
)

reinf_comp <- list(
  c("Neonate (NO IFN) reinfected", "Neonate IFN and RSV reinfected"),
  c("Neonate (NO IFN) reinfected", "Adult reinfected"),
  c("Neonate IFN and RSV reinfected", "Adult reinfected")
)

all_comp <- c(prim_inf_comp, reinf_comp)

names(all_comp) <- c(
  "primary_NeoRSVnoIFN_vs_NeoIFNRSV",
  "primary_NeoRSVnoIFN_vs_AdultRSV",
  "primary_NeoIFNRSV_vs_AdultRSV",
  "reinf_NeonoIFN_vs_NeoIFNRSV",
  "reinf_NeonoIFN_vs_Adult",
  "reinf_NeoIFNRSV_vs_Adult"
)

# =========================
# DEG loop
# =========================
for (ct in celltypes_use) {

  message("\n==============================")
  message("Running cell type: ", ct)
  message("==============================")

  obj_ct <- subset(
    SeuObj,
    subset = minor_celltype == ct
  )

  for (comp_name in names(all_comp)) {

    Group_1 <- all_comp[[comp_name]][1]
    Group_2 <- all_comp[[comp_name]][2]

    message("  Comparison: ", comp_name)

    group_counts <- table(obj_ct$Group)

    if (!(Group_1 %in% names(group_counts)) | !(Group_2 %in% names(group_counts))) {
      message("    Skipped: one group missing")
      next
    }

    if (group_counts[Group_1] == 0 | group_counts[Group_2] == 0) {
      message("    Skipped: one group has 0 cells")
      next
    }

    plt <- Plot_Volcano_DEG_sc(
      SeuratObj = obj_ct,
      assay = "Nanostring",
      Identity = "Group",
      Group_1 = Group_1,
      Group_2 = Group_2,
      logFC_direction = Group_2,
      pCutoff = 0.05,
      FCcutoff = 0.5,
      remove_mt = TRUE,
      genes_to_remove_from_vis = c("Hba-a1/2", "Hbb"),
      title = paste0(ct, " | ", comp_name)
    )

    DEG <- plt[[2]]

    deg_bar <- Plot_Top_DEG_Barplot(
      DEG = DEG,
      Group_1 = Group_1,
      Group_2 = Group_2,
      logFC_direction = Group_2,
      pCutoff = 0.05,
      FCcutoff = 0.5,
      top_n = 20,
      title = paste0(ct, " | ", comp_name)
    )

    pos_deg <- DEG[
      DEG$logFC >= 0.5 & DEG$padj < 0.05,
    ]$feature

    neg_deg <- DEG[
      DEG$logFC <= -0.5 & DEG$padj < 0.05,
    ]$feature

    message("    Positive DEG / higher in: ", Group_2, " = ", length(pos_deg), " genes")
    message("    Negative DEG / higher in: ", Group_1, " = ", length(neg_deg), " genes")

    safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)

    safe_ct <- gsub("_+$", "", gsub("^_+", "", safe_ct))
    
    safe_group2 <- gsub("[^A-Za-z0-9]+", "_", Group_2)
    
    safe_group2 <- gsub("_+$", "", gsub("^_+", "", safe_group2))
    
    base_name <- paste0(
    
      "DEGs_",
    
      safe_ct,
    
      "_",
    
      comp_name,
    
      "_",
    
      safe_group2,
    
      "_direction"
    
    )

    csv_file <- file.path(outdir, paste0(base_name, "_filtered_padj05_logFC05.csv"))
    volcano_pdf_file <- file.path(outdir, paste0(base_name, "_volcano.pdf"))
    barplot_pdf_file <- file.path(outdir, paste0(base_name, "_topDEG_barplot.pdf"))

    write.csv(
      DEG[
        DEG$padj < 0.05 &
          abs(DEG$logFC) >= 0.5,
      ],
      csv_file,
      row.names = FALSE
    )

    pdf(volcano_pdf_file, width = 14, height = 8)
    print(plt[[1]])
    dev.off()

    pdf(barplot_pdf_file, width = 10, height = 8)
    print(deg_bar)
    dev.off()
  }
}


# library(Seurat)
# library(dplyr)
# library(ggplot2)
# 
# # =========================
# # Settings
# # =========================
# celltypes_use <- c("pDC", "cDC")
# 
# outdir <- "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/results/final_celltype/DEGs/Den_sub"
# dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# 
# prim_inf_comp <- list(
#   c("Neonate RSV infected (NO IFN)", "Neonate IFN and RSV infected"),
#   c("Neonate RSV infected (NO IFN)", "Adult RSV infected"),
#   c("Neonate IFN and RSV infected", "Adult RSV infected")
# )
# 
# reinf_comp <- list(
#   c("Neonate (NO IFN) reinfected", "Neonate IFN and RSV reinfected"),
#   c("Neonate (NO IFN) reinfected", "Adult reinfected"),
#   c("Neonate IFN and RSV reinfected", "Adult reinfected")
# )
# 
# all_comp <- c(prim_inf_comp, reinf_comp)
# 
# names(all_comp) <- c(
#   "primary_NeoRSVnoIFN_vs_NeoIFNRSV",
#   "primary_NeoRSVnoIFN_vs_AdultRSV",
#   "primary_NeoIFNRSV_vs_AdultRSV",
#   "reinf_NeonoIFN_vs_NeoIFNRSV",
#   "reinf_NeonoIFN_vs_Adult",
#   "reinf_NeoIFNRSV_vs_Adult"
# )
# 
# # =========================
# # DEG loop
# # =========================
# for (ct in celltypes_use) {
#   
#   message("\n==============================")
#   message("Running cell type: ", ct)
#   message("==============================")
#   
#   obj_ct <- subset(
#     SeuObj,
#     subset = cell_based_final == ct
#   )
#   
#   for (comp_name in names(all_comp)) {
#     
#     Group_1 <- all_comp[[comp_name]][1]
#     Group_2 <- all_comp[[comp_name]][2]
#     
#     message("  Comparison: ", comp_name)
#     
#     group_counts <- table(obj_ct$Group)
#     
#     if (!(Group_1 %in% names(group_counts)) | !(Group_2 %in% names(group_counts))) {
#       message("    Skipped: one group missing")
#       next
#     }
#     
#     if (group_counts[Group_1] == 0 | group_counts[Group_2] == 0) {
#       message("    Skipped: one group has 0 cells")
#       next
#     }
#     
#     plt <- Plot_Volcano_DEG_sc(
#       SeuratObj = obj_ct,
#       assay = "Nanostring",
#       Identity = "Group",
#       Group_1 = Group_1,
#       Group_2 = Group_2,
#       logFC_direction = Group_2,
#       pCutoff = 0.05,
#       FCcutoff = 0.5,
#       remove_mt = TRUE,
#       genes_to_remove_from_vis = c("Hba-a1/2", "Hbb"),
#       title = paste0(ct, " | ", comp_name)
#     )
#     
#     DEG <- plt[[2]]
#     
#     deg_bar <- Plot_Top_DEG_Barplot(
#       DEG = DEG,
#       Group_1 = Group_1,
#       Group_2 = Group_2,
#       logFC_direction = Group_2,
#       pCutoff = 0.05,
#       FCcutoff = 0.5,
#       top_n = 20,
#       title = paste0(ct, " | ", comp_name)
#     )
#     
#     pos_deg <- DEG[
#       DEG$logFC >= 0.5 & DEG$padj < 0.05,
#     ]$feature
#     
#     neg_deg <- DEG[
#       DEG$logFC <= -0.5 & DEG$padj < 0.05,
#     ]$feature
#     
#     message("    Positive DEG / higher in: ", Group_2, " = ", length(pos_deg), " genes")
#     message("    Negative DEG / higher in: ", Group_1, " = ", length(neg_deg), " genes")
#     
#     base_name <- paste0(
#       "DEGs_",
#       ct,
#       "_",
#       comp_name,
#       "_",
#       gsub("[^A-Za-z0-9]+", "_", Group_2),
#       "_direction"
#     )
#     
#     csv_file <- file.path(outdir, paste0(base_name, "_filtered_padj05_logFC05.csv"))
#     volcano_pdf_file <- file.path(outdir, paste0(base_name, "_volcano.pdf"))
#     barplot_pdf_file <- file.path(outdir, paste0(base_name, "_topDEG_barplot.pdf"))
#     
#     write.csv(
#       DEG[
#         DEG$padj < 0.05 &
#           abs(DEG$logFC) >= 0.5,
#       ],
#       csv_file,
#       row.names = FALSE
#     )
#     
#     pdf(volcano_pdf_file, width = 14, height = 8)
#     print(plt[[1]])
#     dev.off()
#     
#     pdf(barplot_pdf_file, width = 10, height = 8)
#     print(deg_bar)
#     dev.off()
#   }
# }



