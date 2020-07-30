library(Seurat)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
files <- list.files(".", pattern = "*.RDS", full.names = TRUE)
combined_obj <- readRDS(files[1]) #combined seurat list

#all_obj <- list("all_data_sets" = combined_obj)

print(object.size(combined_obj), units = "MB") #before

for (i in 1:length(combined_obj)) {
  combined_obj[[i]]@assays$RNA@counts <- matrix()
  combined_obj[[i]]@assays$RNA@scale.data <- matrix()
  combined_obj[[i]]@assays$integrated@scale.data <- matrix()
  print(paste0("saving: ", names(combined_obj)[i]))
  saveRDS(combined_obj[[i]], paste0("TRIMMED_", names(combined_obj)[i], ".RDS"))
}

print(object.size(combined_obj), units = "MB") #after
