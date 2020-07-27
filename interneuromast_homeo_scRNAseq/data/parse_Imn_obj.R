library(Seurat)
library(ggplot2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

app_name <- "interneuromast_homeo_scRNAseq"

obj_integrated <- readRDS("./scaled.filtered_adj_fpkm_1828_smartseq_integ.RDS")

files <- list.files(".", pattern = "scaled", full.names = TRUE)

file_list <- list()

print("Loading Seurat objects...")
for (i in 1:length(files)) {
  file_list[[i]] <- readRDS(files[i])
  print(object.size(file_list[[i]]), units = "MB")
  
  DefaultAssay(file_list[[i]]) <- "integrated"
  file_list[[i]] <- subset(file_list[[i]], idents = c("Inm", "mantle-cells") , subset = seq.method == "10X")
  file_list[[i]] <- ScaleData(file_list[[i]], features = rownames(file_list[[i]]))
  file_list[[i]] <- RunPCA(file_list[[i]])
  file_list[[i]] <- FindNeighbors(file_list[[i]], dims = 1:10)
  file_list[[i]] <- FindClusters(file_list[[i]], resolution = 0.6)
  file_list[[i]] <- RunUMAP(file_list[[i]], dims = 1:10)
  #omit unfixed sample
  file_list[[i]] <- subset(file_list[[i]], subset = data.set == c("homeo-10X-isl1", "homeo-10X-2410-8"))
  
  file_list[[i]]$data.set <- droplevels(file_list[[i]]$data.set) #drop unused data.sets (smartseq cells)
  file_list[[i]]$cell.type.ident <- droplevels(file_list[[i]]$cell.type.ident) #drop unused cell type idents
  print(object.size(file_list[[i]]), units = "MB")
  }

DimPlot(file_list[[1]], reduction = "umap")

head(file_list[[1]][["RNA"]]@var.features)

ElbowPlot(file_list[[1]])

saveRDS(file_list[[1]], paste0(app_name, ".RDS"))

# ===================================== Trim Obj ====================================

# Processing for individual data set that is recently added


files <- list.files(".", pattern = "*.RDS", full.names = TRUE)
new_obj <- readRDS(files[1])
print(format(object.size(new_obj), units = "Mb"))


new_obj@assays$RNA@counts <- matrix()
new_obj@assays$RNA@scale.data <- matrix()
new_obj@assays$integrated@counts <- matrix()
new_obj@assays$integrated@scale.data <- matrix()
saveRDS(new_obj, paste0("TRIMMED_", app_name, ".RDS"))
