library(Seurat)
library(dplyr)
library(future)
options(future.globals.maxSize = 1000 * 1024^2)


#setwd("./bgmp/shiny-apps-main/smrtseq_vs_10X_scRNAseq/data/seurat_obj_input/")
setwd("/Volumes/easystore/SIMR_2019/shiny_backup/shiny-apps-main/smrtseq_vs_10X_scRNAseq")
seurat_obj <- readRDS("./data/scaled.filtered_adj_fpkm_1828_smartseq_integ.RDS")

# ======================================= ===========================
Idents(seurat_obj) <- "seq.method"
smartseq <- WhichCells(seurat_obj, idents = "smartseq2")
tenX <- WhichCells(seurat_obj, idents = "10X", downsample = round(length(colnames(seurat_obj)) * .30))
sampled_seurat_obj <- subset(seurat_obj, cells = c(smartseq, tenX))
print(object.size(sampled_seurat_obj), units = "MB")

levels(sampled_seurat_obj$cell.type.ident)



