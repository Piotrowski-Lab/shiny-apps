#seurat_obj <- readRDS("./data/filtered_adj_fpkm_1828_smartseq_integ.RDS")


files <- list.files("./data", pattern = ".RDS", full.names = TRUE)
file_list <- list()

print("Loading Seurat objects...")
for (i in 1:length(files)) {
  file_list[[i]] <- readRDS(files[i])
  print(object.size(file_list[[i]]), units = "MB")

  DefaultAssay(file_list[[i]]) <- "RNA"
  file_list[[i]] <- FindVariableFeatures(file_list[[i]], selection.method = "vst", nfeatures = 2000)
  file_list[[i]] <- ScaleData(file_list[[i]])
}

#delete integrated matrix
file_list[[1]][["integrated"]]@scale.data <- matrix() 

print(object.size(file_list[[1]]), units = "MB")


saveRDS(file_list[[1]], "./data/seurat_obj_input/scaled.filtered_adj_fpkm_1828_smartseq_integ.RDS")

# print(object.size(seurat_obj), units = "MB")
# #10238.8 mb
# 
# class(seurat_obj[["RNA"]]@scale.data)
# 
# x <- as.sparse(seurat_obj[["RNA"]]@scale.data)
# head(x)
# #seurat_obj[["RNA"]]@scale.data <- as.sparse(seurat_obj[["RNA"]]@scale.data)
# 
# as.sparse(seurat_obj[["RNA"]]@scale.data)
# 
# head(seurat_obj[["RNA"]]@scale.data)
# 
# seurat_obj[["integrated"]]@scale.data <- matrix() 
# print(object.size(seurat_obj), units = "MB")
# #9802.9 Mb
# 
# seurat_obj[["RNA"]]@var.features
# 
# seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
# 
# seurat_obj <- ScaleData(seurat_obj,verbose = TRUE)
# print(object.size(seurat_obj), units = "MB")
# #3079.5 Mb
# 
