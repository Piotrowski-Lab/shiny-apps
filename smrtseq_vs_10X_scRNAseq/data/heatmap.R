library(Seurat)
library(ggplot2)

# ======================================================== function ===================================
obj_integrated@meta.data$data.set <- plyr::revalue(
  obj_integrated@meta.data$data.set, c("homeo-2047" = "homeo-10X-2047",
                                       "homeo-2410-7" = "homeo-10X-2410-7",
                                       "homeo-2410-8" = "homeo-10X-2410-8"))

saveRDS(obj_integrated, "./scaled.filtered_adj_fpkm_1828_smartseq_integ.RDS")

#ids <- as.list(levels(obj_integrated$data.set))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

obj_integrated <- readRDS("./subsampled_30_smrtseq_10X.RDS")

smartseq_tomatch <- c("1hr-smrtseq", "homeo-smrtseq")

tenX_tomatch <- c("homeo-10X-isl1", "homeo-10X-2047", "homeo-10X-2410-7", "homeo-10X-2410-8")

split_heatmap <- function(seurat_obj, method, tomatch){
  split_obj <- subset(seurat_obj, subset = seq.method == method)
 
  split_obj$data.set <- droplevels(split_obj$data.set)
  
  split_obj$cell.type.ident <- droplevels(split_obj$cell.type.ident)
  
  
  return(split_obj)
}

smartseq <- split_heatmap(obj_integrated, method = "smartseq2", tomatch = smartseq_tomatch)

 s <- DoHeatmap(smartseq, features = features, group.by = "data.set")
 
tenX <- split_heatmap(obj_integrated, method = "10X", tomatch = tenX_tomatch)

t <- DoHeatmap(tenX, features = features, group.by = "data.set")
 
s + t

cowplot::plot_grid(s, t, nrow = 2)
