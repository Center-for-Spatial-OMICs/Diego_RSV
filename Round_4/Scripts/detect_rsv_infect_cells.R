library(Seurat); library(tidyverse); library(data.table); library(dplyr)
library(RColorBrewer); library(InSituType)
library(ggplot2)
library(dplyr)
library(viridis)
library(dittoSeq)
library(shiny)
library(plotly)
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(qs)



## Load data
Sobj <- qread("/mnt/scratch1/maycon/Diego_RSV_CosMx/Round_4/Objects/current_Sobj.qs")



## Load list of custom genes 
# The "RSV related genes" ought to be among those
custom_genes <- c("Acp5", "Aicda", "Atp5b", "Bcl6", "Birc3", "Braf", "Ccdc80", "Ccl17", 
           "Ccl22", "Cc26", "Ccr9", "Cd1d1", "Cd22", "Cd5", "Cd93", "Ccr2", 
           "Ctsd", "Cxcl13", "Dpt", "Fap", "Flt3", "Gsk3b", "Havcr1", "Hhex",
           "Il21", "Il21r", "Isg15", "Itm2b", "Map2k1", "Marcksl1", "Ms4a4a",
           "Ncam1", "Ndufa4l2", "Nhej1", "Pax5", "Pgk1", "Plcg1", "Pld3",
           "Prdm1", "Prox1", "Pten", "Rbm47", "Sdc1", "Sla", "Spn",
           "Tgfbi", "Tnfrsf13c", "Tnfsf13", "Tpi1")


## Plot heatmap with custom genes only
rownames(Sobj) %in% custom_genes %>% table() 
setdiff(custom_genes , rownames(Sobj) ) # we just dont have "Cc26" on this CosMx panel 

































