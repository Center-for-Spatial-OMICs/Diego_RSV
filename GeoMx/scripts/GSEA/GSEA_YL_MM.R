library(limma)
library(clusterProfiler)
library(enrichplot)
library(dplyr)

library(ggplot2)

library(patchwork)

library(ggrepel)

# check all available DEG names
deg_list <- readRDS("/mnt/scratch1/maycon/Diego_RSV_CosMx/GeoMx_from_Yutian/results/DEG/Reinfection_pairwise_DEG/Reinfected_pairwise_DEG_list.rds")
cat("All DEG names:\n")
print(names(deg_list))

# =========================
# Select primary infection comparisons
# =========================

# # keep only:
# # INR_vs_NR
# # AR_vs_NR
# # INR_vs_AR
# 
# primary_idx <- grep(
#   "^(INR_vs_NR|AR_vs_NR|INR_vs_AR)\\s*=",
#   names(deg_list)
# )
# 
# primary_deg <- deg_list[primary_idx]
# 
# # rename to clean short names
# names(primary_deg) <- sub(
#   "\\s*=.*$",
#   "",
#   names(primary_deg)
# )
# 
# cat("Primary infection comparisons selected:\n")
# print(names(primary_deg))



# =========================
# Safer GSEA function
# =========================

run_gsea_one <- function(deg, comp_name, gene_sets, outdir) {
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(deg)) {
    message("Skipping NULL object: ", comp_name)
    return(NULL)
  }
  
  if (missing(gene_sets) || is.null(gene_sets)) {
    stop(
      "gene_sets must be provided for GSEA. ",
      "Use a list of TERM2GENE data.frames."
    )
  }
  
  message("Running GSEA: ", comp_name)
  
  deg <- as.data.frame(deg)
  
  deg <- deg %>%
    dplyr::filter(
      !is.na(gene),
      !is.na(logFC)
    )
  
  if (exists("contrast_defs") && comp_name %in% names(contrast_defs)) {
    
    contrast_groups <- strsplit(contrast_defs[comp_name], " - ")[[1]]
    
    positive_group <- contrast_groups[1]
    negative_group <- contrast_groups[2]
    
    direction_title <- paste0(
      comp_name,
      " | Positive NES/logFC = higher in ",
      positive_group,
      " ; Negative NES/logFC = higher in ",
      negative_group
    )
    
  } else if (grepl(" = ", comp_name) && grepl(" - ", comp_name)) {
    
    contrast_def <- sub(".* = ", "", comp_name)
    contrast_groups <- strsplit(contrast_def, " - ")[[1]]
    
    positive_group <- contrast_groups[1]
    negative_group <- contrast_groups[2]
    
    direction_title <- paste0(
      comp_name,
      " | Positive NES/logFC = higher in ",
      positive_group,
      " ; Negative NES/logFC = higher in ",
      negative_group
    )
    
  } else {
    
    direction_title <- comp_name
  }
  
  gene_list <- deg$logFC
  names(gene_list) <- deg$gene
  
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  cat("Genes used:", length(gene_list), "\n")
  
  if (length(gene_list) < 10) {
    message("Too few genes for: ", comp_name)
    return(NULL)
  }
  
  safe_comp_name <- gsub("[^A-Za-z0-9_]+", "_", comp_name)
  
  for (set_name in names(gene_sets)) {
    
    message("  Database: ", set_name)
    
    gsea_res <- tryCatch({
      
      clusterProfiler::GSEA(
        geneList = gene_list,
        TERM2GENE = gene_sets[[set_name]],
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        eps = 0,
        verbose = FALSE
      )
      
    }, error = function(e) {
      
      message("GSEA failed: ", comp_name, " | ", set_name)
      message(e$message)
      return(NULL)
    })
    
    if (!is.null(gsea_res)) {
      
      gsea_df <- as.data.frame(gsea_res)
      
      if (nrow(gsea_df) == 0) {
        message("No enriched pathways found.")
        next
      }
      
      gsea_df$comparison <- comp_name
      gsea_df$database <- set_name
      gsea_df$direction <- direction_title
      
      write.csv(
        gsea_df,
        file.path(
          outdir,
          paste0(safe_comp_name, "_", set_name, "_GSEA.csv")
        ),
        row.names = FALSE
      )
      
      p <- enrichplot::dotplot(
        gsea_res,
        showCategory = 20
      ) +
        ggplot2::ggtitle(
          paste0(direction_title, "\n", set_name, " GSEA")
        ) +
        ggplot2::theme_bw(base_size = 6) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 7),
          axis.title = ggplot2::element_text(size = 6),
          axis.text.x = ggplot2::element_text(size = 5),
          axis.text.y = ggplot2::element_text(size = 5),
          legend.title = ggplot2::element_text(size = 6),
          legend.text = ggplot2::element_text(size = 5)
        )
      
      ggplot2::ggsave(
        file.path(
          outdir,
          paste0(safe_comp_name, "_", set_name, "_GSEA_dotplot.png")
        ),
        p,
        width = 9,
        height = 7,
        dpi = 300
      )
      
      top_pathway <- gsea_df$ID[1]
      
      p2 <- enrichplot::gseaplot2(
        gsea_res,
        geneSetID = top_pathway,
        title = paste0(
          top_pathway,
          "\n",
          direction_title
        ),
        base_size = 6
      )
      
      png(
        filename = file.path(
          outdir,
          paste0(safe_comp_name, "_", set_name, "_top_pathway.png")
        ),
        width = 9,
        height = 6,
        units = "in",
        res = 300
      )
      
      print(p2)
      
      dev.off()
    }
  }
}

# =========================
# Run primary infection GSEA
# =========================
library(msigdbr)
library(dplyr)

# =========================
# MSigDB gene sets
# =========================

hallmark_t2g <- msigdbr(
  species = "Mus musculus",
  category = "H"
) %>%
  dplyr::select(gs_name, gene_symbol)

reactome_t2g <- msigdbr(
  species = "Mus musculus",
  category = "C2",
  subcategory = "CP:REACTOME"
) %>%
  dplyr::select(gs_name, gene_symbol)

go_bp_t2g <- msigdbr(
  species = "Mus musculus",
  category = "C5",
  subcategory = "GO:BP"
) %>%
  dplyr::select(gs_name, gene_symbol)

immune_t2g <- msigdbr(
  species = "Mus musculus",
  category = "C7"
) %>%
  dplyr::select(gs_name, gene_symbol)

# =========================
# Bundle for run_gsea_one()
# =========================

gene_sets <- list(
  Hallmark = hallmark_t2g,
  Reactome = reactome_t2g,
  GO_BP = go_bp_t2g,
  Immune = immune_t2g
)

# quick summary
sapply(gene_sets, function(x) length(unique(x$gs_name)))

for (nm in names(deg_list)) {
  
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", nm)
  
  run_gsea_one(
    deg = deg_list[[nm]],
    comp_name = nm,
    gene_sets = gene_sets,
    outdir = file.path(
      "/mnt/scratch1/maycon/Diego_RSV_CosMx/GeoMx_from_Yutian/results/GSEA",
      safe_name
    )
  )
}


# Maycon: I stopped here. ----------------









# =========================
# Focused immune/epithelial pathway filtering
# =========================

keywords <- c(
  # B cell
  "B_CELL", "B CELL", "BCR", "B_CELL_RECEPTOR",
  "B_CELL_ACTIVATION", "B_CELL_DIFFERENTIATION",
  "PLASMA_CELL", "GERMINAL_CENTER",
  
  # CD8 / T cell
  "T_CELL", "T CELL", "CD8", "CYTOTOXIC",
  "EFFECTOR", "GRANZYME", "PERFORIN",
  
  # antigen presentation / IFN
  "ANTIGEN", "PRESENTATION", "MHC",
  "INTERFERON", "IFN", "JAK", "STAT",
  
  # epithelial injury / repair / inflammation
  "EPITHELIAL", "EPITHELIUM", "REPAIR", "WOUND",
  "MIGRATION", "INFLAMMATORY", "INFLAMMATION",
  "TNF", "NFKB", "IL6", "CYTOKINE"
)

gsea_files <- list.files(
  gsea_outdir,
  pattern = "_GSEA\\.csv$",
  full.names = TRUE
)

focused_gsea <- lapply(gsea_files, function(f) {
  
  df <- read.csv(f)
  
  df %>%
    filter(
      grepl(
        paste(keywords, collapse = "|"),
        ID,
        ignore.case = TRUE
      )
    ) %>%
    mutate(source_file = basename(f))
  
}) %>%
  bind_rows() %>%
  filter(!is.na(NES)) %>%
  arrange(p.adjust)

write.csv(
  focused_gsea,
  file.path(gsea_outdir, "Primary_infection_focused_immune_epithelial_GSEA.csv"),
  row.names = FALSE
)

top_pathways <- focused_gsea %>%
  group_by(ID) %>%
  summarise(best_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
  arrange(best_padj) %>%
  slice_head(n = 35) %>%
  pull(ID)

plot_df <- focused_gsea %>%
  filter(ID %in% top_pathways) %>%
  mutate(
    comparison = factor(comparison, levels = c("INR_vs_NR", "AR_vs_NR", "INR_vs_AR")),
    sig = ifelse(p.adjust < 0.05, "*", ""),
    neglog10_padj = -log10(p.adjust)
  )

p_heat <- ggplot(plot_df, aes(x = comparison, y = ID, fill = NES)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 5) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Primary Infection Focused GSEA",
    x = "Comparison",
    y = "Immune / epithelial pathway",
    fill = "NES"
  )

ggsave(
  file.path(gsea_outdir, "Primary_infection_focused_GSEA_heatmap.png"),
  p_heat,
  width = 12,
  height = 10,
  dpi = 300
)


library(dplyr)
library(ggplot2)
library(stringr)

gsea_outdir <- "/media/ResearchHome/plummgrp/home/common/CSO-Dir/Files/StJude/GeoMx/DiegoGrp/Data/Diego-RSV-Processed/Manual_Annotation/2026/Main_pairwise_DEG/Primary_infection_GSEA"

gsea_files <- list.files(
  gsea_outdir,
  pattern = "_GSEA\\.csv$",
  full.names = TRUE
)

for (f in gsea_files) {
  
  df <- read.csv(f)
  
  if (nrow(df) == 0) next
  
  plot_df <- df %>%
    filter(!is.na(NES), !is.na(p.adjust)) %>%
    arrange(p.adjust) %>%
    slice_head(n = 12) %>%
    mutate(
      pathway = ifelse(
        !is.na(Description) & Description != "",
        Description,
        ID
      ),
      pathway = str_replace_all(pathway, "^GOBP_|^HALLMARK_|^REACTOME_", ""),
      pathway = str_replace_all(pathway, "_", " "),
      pathway = str_wrap(pathway, width = 45),
      pathway = factor(pathway, levels = rev(pathway)),
      neglog10_fdr = -log10(p.adjust)
    )
  
  if (nrow(plot_df) == 0) next
  
  title_name <- basename(f) %>%
    str_replace("_GSEA.csv", "") %>%
    str_replace_all("_", " ")
  
  p <- ggplot(plot_df, aes(x = NES, y = pathway)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(
      aes(size = setSize, color = neglog10_fdr),
      alpha = 0.9
    ) +
    scale_color_gradient(
      low = "steelblue",
      high = "red"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 8.5),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      plot.margin = margin(10, 20, 10, 10)
    ) +
    labs(
      title = title_name,
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      color = "-log10(FDR)",
      size = "Pathway size"
    )
  
  outfile <- file.path(
    gsea_outdir,
    paste0(
      basename(f) %>% str_replace("_GSEA.csv", ""),
      "_clean_nonoverlap_dotplot.png"
    )
  )
  
  ggsave(
    outfile,
    p,
    width = 12,
    height = 8,
    dpi = 300
  )
}

cat("Done. Clean non-overlap plots saved in:\n")
cat(gsea_outdir, "\n")

p_heat