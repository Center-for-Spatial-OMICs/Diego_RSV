library(limma)

library(dplyr)

library(ggplot2)

library(patchwork)

library(ggrepel)

# # don't exist
# indir <- "/media/ResearchHome/plummgrp/home/common/CSO-Dir/Files/StJude/GeoMx/DiegoGrp/Data/Diego-RSV-Processed/Manual_Annotation/2026"

outdir <- file.path("/mnt/scratch1/maycon/Diego_RSV_CosMx/GeoMx_from_Yutian/results/DEG", "Reinfection_pairwise_DEG")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================

# Read expression matrix

# =========================

# =========================
# Read expression matrix
# =========================
df2 <- read.csv(
  file.path("/media/ResearchHome/plummgrp/home/common/CSO-Dir/Files/StJude/GeoMx/DiegoGrp/Data/Diego-RSV-Processed/Data_RAW/ex_df_diego_rsv.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# first column is gene name, even if not named X
gene_names <- make.unique(as.character(df2[[1]]))

expr <- df2[, -1, drop = FALSE]
expr <- as.data.frame(lapply(expr, function(x) as.numeric(as.character(x))))
expr <- as.matrix(expr)

rownames(expr) <- gene_names

colnames(expr) <- colnames(df2)[-1] |>
  gsub("\\.dcc$", "", x = _) |>
  gsub("\\.", "-", x = _)

cat("Gene check:\n")
print(head(rownames(expr)))

# =========================
# Read metadata
# =========================
adult1 <- read.csv(file.path("/media/ResearchHome/plummgrp/home/common/CSO-Dir/Files/StJude/GeoMx/DiegoGrp/Data/Diego-RSV-Processed/Metadata", "adult1_annotaion.csv"), check.names = FALSE)
adult2 <- read.csv(file.path("/media/ResearchHome/plummgrp/home/common/CSO-Dir/Files/StJude/GeoMx/DiegoGrp/Data/Diego-RSV-Processed/Metadata", "adult2_annotaion.csv"), check.names = FALSE)
neonate1 <- read.csv(file.path("/media/ResearchHome/plummgrp/home/common/CSO-Dir/Files/StJude/GeoMx/DiegoGrp/Data/Diego-RSV-Processed/Metadata", "neonate1_annotaion.csv"), check.names = FALSE)
neonate2 <- read.csv(file.path("/media/ResearchHome/plummgrp/home/common/CSO-Dir/Files/StJude/GeoMx/DiegoGrp/Data/Diego-RSV-Processed/Metadata", "neonate2_annotaion.csv"), check.names = FALSE)

meta <- bind_rows(adult1, adult2, neonate1, neonate2)

meta$Sample_ID_clean <- meta$Sample_ID |>
  gsub("\\.dcc$", "", x = _) |>
  gsub("\\.", "-", x = _)

cat("Original groups:\n")
print(table(meta$Group, useNA = "ifany"))

# =========================
# Match samples
# =========================
common_ids <- intersect(colnames(expr), meta$Sample_ID_clean)

cat("Matched samples:", length(common_ids), "\n")

expr_sub <- expr[, common_ids, drop = FALSE]
meta_sub <- meta[match(common_ids, meta$Sample_ID_clean), ]

stopifnot(all(colnames(expr_sub) == meta_sub$Sample_ID_clean))

expr_log <- log2(expr_sub + 1)

cat("Expression rowname check:\n")
print(head(rownames(expr_log)))

# =========================
# Clean biological groups
# =========================
meta_sub$Group_clean <- case_when(
  
  meta_sub$Group == "Adult Control" ~ "Adult_Control",
  meta_sub$Group == "Adult RSV infected" ~ "Adult_RSV",
  meta_sub$Group == "Adult reinfected" ~ "Adult_Reinfected",
  
  meta_sub$Group == "Neonate control" ~ "Neonate_Control",
  
  grepl("Neonate", meta_sub$Group, ignore.case = TRUE) &
    grepl("RSV infected", meta_sub$Group, ignore.case = TRUE) &
    grepl("NO IFN", meta_sub$Group, ignore.case = TRUE) &
    !grepl("reinfect", meta_sub$Group, ignore.case = TRUE) ~ "Neonate_RSV_NO_IFN",
  
  meta_sub$Group == "Neonate IFN and RSV infected" ~ "Neonate_RSV_IFN",
  
  grepl("Neonate", meta_sub$Group, ignore.case = TRUE) &
    grepl("NO IFN", meta_sub$Group, ignore.case = TRUE) &
    grepl("reinfect", meta_sub$Group, ignore.case = TRUE) ~ "Neonate_Reinfected_NO_IFN",
  
  meta_sub$Group == "Neonate IFN and RSV reinfected" ~ "Neonate_Reinfected_IFN",
  
  TRUE ~ NA_character_
)

cat("Clean groups:\n")
print(table(meta_sub$Group_clean, useNA = "ifany"))

# =========================
# Filter valid groups
# =========================
keep <- !is.na(meta_sub$Group_clean)

expr_log2 <- expr_log[, keep, drop = FALSE]
meta_sub2 <- meta_sub[keep, ]

group <- factor(meta_sub2$Group_clean)

cat("Final groups:\n")
print(table(group))

# =========================
# Design matrix
# =========================
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

cat("Design columns:\n")
print(colnames(design))

# =========================
# Run limma
# =========================
fit <- lmFit(expr_log2, design)

contrast_list <- c(

  # "AR_vs_AC = Adult_RSV - Adult_Control",
  # "ARR_vs_AC = Adult_Reinfected - Adult_Control",
  # "ARR_vs_AR = Adult_Reinfected - Adult_RSV",

  # "NC_vs_AC = Neonate_Control - Adult_Control",
  # "NR_vs_AC = Neonate_RSV_NO_IFN - Adult_Control",
  # "INR_vs_AC = Neonate_RSV_IFN - Adult_Control",
  # "NRR_vs_AC = Neonate_Reinfected_NO_IFN - Adult_Control",
  # "INRR_vs_AC = Neonate_Reinfected_IFN - Adult_Control",

  # "NC_vs_AR = Neonate_Control - Adult_RSV",
  # "NR_vs_AR = Neonate_RSV_NO_IFN - Adult_RSV",
  # "INR_vs_AR = Neonate_RSV_IFN - Adult_RSV",
  # "NRR_vs_AR = Neonate_Reinfected_NO_IFN - Adult_RSV",
  # "INRR_vs_AR = Neonate_Reinfected_IFN - Adult_RSV",

  #"NC_vs_ARR = Neonate_Control - Adult_Reinfected",
  #"NR_vs_ARR = Neonate_RSV_NO_IFN - Adult_Reinfected",
  #"INR_vs_ARR = Neonate_RSV_IFN - Adult_Reinfected",
  "NRR_vs_ARR = Neonate_Reinfected_NO_IFN - Adult_Reinfected",
  "INRR_vs_ARR = Neonate_Reinfected_IFN - Adult_Reinfected",

  # "NR_vs_NC = Neonate_RSV_NO_IFN - Neonate_Control",
  # "INR_vs_NC = Neonate_RSV_IFN - Neonate_Control",
  # "NRR_vs_NC = Neonate_Reinfected_NO_IFN - Neonate_Control",
  # "INRR_vs_NC = Neonate_Reinfected_IFN - Neonate_Control",

  # "INR_vs_NR = Neonate_RSV_IFN - Neonate_RSV_NO_IFN",
  # "NRR_vs_NR = Neonate_Reinfected_NO_IFN - Neonate_RSV_NO_IFN",
  # "INRR_vs_NR = Neonate_Reinfected_IFN - Neonate_RSV_NO_IFN",

  # "NRR_vs_INR = Neonate_Reinfected_NO_IFN - Neonate_RSV_IFN",
  # "INRR_vs_INR = Neonate_Reinfected_IFN - Neonate_RSV_IFN",

  "INRR_vs_NRR = Neonate_Reinfected_IFN - Neonate_Reinfected_NO_IFN"
)

contrast <- makeContrasts(
  contrasts = contrast_list,
  levels = design
)

contrast_defs <- sub(".* = ", "", contrast_list)
names(contrast_defs) <- colnames(contrast)

fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# =========================
# DEG tables
# =========================
deg_list_main <- list()

for (coef_name in colnames(contrast)) {
  
  tt <- topTable(
    fit2,
    coef = coef_name,
    number = Inf,
    sort.by = "none"
  )
  
  tt$gene <- rownames(expr_log2)
  
  contrast_groups <- strsplit(contrast_defs[coef_name], " - ")[[1]]
  
  higher_group <- contrast_groups[1]
  lower_group <- contrast_groups[2]
  
  tt$change <- "Not Sig"
  tt$change[tt$adj.P.Val < 0.05 & tt$logFC > 0.5] <- paste0("Higher in ", higher_group)
  tt$change[tt$adj.P.Val < 0.05 & tt$logFC < -0.5] <- paste0("Higher in ", lower_group)
  
  deg <- tt %>%
    arrange(adj.P.Val)
  
  deg_list_main[[coef_name]] <- deg
  
  write.csv(
    deg,
    file.path(outdir, paste0(coef_name, "_DEG.csv")),
    row.names = FALSE
  )
}

saveRDS(
  deg_list_main,
  file.path(outdir, "Reinfected_pairwise_DEG_list.rds")
)

# =========================
# Volcano plot function
# =========================
plot_volcano <- function(deg, title, top_n = 15) {
  
  deg <- deg %>%
    mutate(
      adj.P.Val = ifelse(is.na(adj.P.Val), 1, adj.P.Val),
      adj.P.Val = ifelse(adj.P.Val <= 0, 1e-300, adj.P.Val),
      neglog10 = -log10(adj.P.Val)
    )
  
  direction_labels <- setdiff(unique(deg$change), "Not Sig")
  
  top_up <- deg %>%
    filter(change == direction_labels[1]) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = top_n)
  
  top_down <- deg %>%
    filter(change == direction_labels[2]) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = top_n)
  
  label_df <- bind_rows(top_up, top_down)
  
  ggplot(
    deg,
    aes(x = logFC, y = neglog10, color = change)
  ) +
    
    geom_point(
      alpha = 0.7,
      size = 1.5
    ) +
    
    geom_vline(
      xintercept = c(-0.5, 0.5),
      linetype = "dashed",
      color = "grey50"
    ) +
    
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey50"
    ) +
    
    geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 3.2,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50"
    ) +
    
    scale_color_manual(
      values = c(
        setNames("#D73027", direction_labels[1]),
        setNames("#4575B4", direction_labels[2]),
        "Not Sig" = "grey80"
      )
    ) +
    
    theme_classic(base_size = 13) +
    
    labs(
      title = title,
      subtitle = "Top significant genes labeled",
      x = "log2 Fold Change",
      y = "-log10 adjusted P value"
    ) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      legend.position = "top"
    )
}

# =========================
# Generate plots
# =========================
plots <- lapply(
  names(deg_list_main),
  function(nm) {
    plot_volcano(
      deg_list_main[[nm]],
      nm,
      top_n = 15
    )
  }
)

# =========================
# Save combined PDF
# =========================
pdf(
  file.path(outdir, "Reinfected_pairwise_volcano_plots_labeled.pdf"),
  width = 16,
  height = 12
)

print(
  wrap_plots(plots, ncol = 2)
)

dev.off()

# =========================
# Save individual PNGs
# =========================
for (nm in names(deg_list_main)) {
  
  p <- plot_volcano(
    deg_list_main[[nm]],
    nm,
    top_n = 15
  )
  
  ggsave(
    file.path(outdir, paste0(nm, "_volcano_labeled.png")),
    p,
    width = 8,
    height = 6.5,
    dpi = 300
  )
}

cat("Done!\n")
cat("Results saved to:\n")
cat(outdir, "\n")