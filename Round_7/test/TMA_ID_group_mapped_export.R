merged_obj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/data/Processed/merged_Sobj_UMAP.qs")


meta <- merged_obj@meta.data[, c("TMA_4", "TMA_3", "TMA_2", "TMA", "Group")]


library(dplyr)
meta_unique <- meta %>%
  distinct(TMA_4, .keep_all = TRUE)

write.csv(meta_unique, "/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_7/test/TMA_ID_group_mapped.csv")
