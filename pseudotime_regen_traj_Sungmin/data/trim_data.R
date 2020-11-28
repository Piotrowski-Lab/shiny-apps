library(Seurat)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# =============================================================================
# Processing for individual seurat obj that is added

files <- list.files(".", pattern = "*.RDS", full.names = TRUE)
new_obj <- readRDS(files[4])
print(format(object.size(new_obj), units = "Mb"))

new_obj@assays$RNA@counts <- matrix()
new_obj@assays$RNA@scale.data <- matrix()
new_obj@assays$integrated@counts <- matrix()
new_obj@assays$integrated@scale.data <- matrix()
saveRDS(new_obj, paste0("TRIMMED_", 
  "SeurObj_ptime_subset_experiments_seurat3_10hr_central_v1.0_dim_6_", ".RDS"))

print(format(object.size(new_obj), units = "Mb"))
