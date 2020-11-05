library(Seurat)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

files <- list.files(".", pattern = "*.RDS", full.names = TRUE)

# =============================================== single object input

combined_obj <- readRDS(files[1]) #combined seurat list

combined_obj <- list("SeurObj_anchored_cell_type_update_nmast_skin_additional_timepoints_seurat3_v1.0_" = combined_obj)

print(object.size(combined_obj), units = "MB") #before

for (i in 1:length(combined_obj)) {
  combined_obj[[i]]@assays$RNA@counts <- matrix()
  combined_obj[[i]]@assays$RNA@scale.data <- matrix()
  combined_obj[[i]]@assays$integrated@scale.data <- matrix()
  print(paste0("saving: ", names(combined_obj)[i]))
  saveRDS(combined_obj[[i]], paste0("TRIMMED_", names(combined_obj)[i], ".RDS"))
}

print(object.size(combined_obj), units = "MB") #after


# =================================================== add new objects with list
files <- files[2]
combined_obj <- list()

for (i in 1:length(files)){
  print(files[i])
  combined_obj[[i]] <- readRDS(files[i]) #combined seurat list
}
combined_obj <- list("SeurObj_ptime_subset_experiments_seurat3_v1.4.1_" = combined_obj[[1]])

print(object.size(combined_obj), units = "MB") #before

for (i in 1:length(combined_obj)) {
  combined_obj[[i]]@assays$RNA@counts <- matrix()
  combined_obj[[i]]@assays$RNA@scale.data <- matrix()
  combined_obj[[i]]@assays$integrated@scale.data <- matrix()
  print(paste0("saving: ", names(combined_obj)[i]))
  saveRDS(combined_obj[[i]], paste0("TRIMMED_", names(combined_obj)[i], ".RDS"))
}

print(object.size(combined_obj), units = "MB") #after
