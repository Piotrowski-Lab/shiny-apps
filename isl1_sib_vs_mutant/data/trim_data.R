library(Seurat)
files <- list.files(".", pattern = "*.RDS", full.names = TRUE)
combined_obj <- readRDS(files[1])
seurat_obj <- readRDS(files[2])

central_obj <- combined_obj[[1]]
HCprog_obj <- combined_obj[[2]]
AP_obj <- combined_obj[[3]]
mantle_obj <- combined_obj[[4]]

all_obj <- list("all_data_sets" = seurat_obj, "central_obj" = central_obj,
  "HCprog_obj" = HCprog_obj,"AP_obj" = AP_obj, "mantle_obj" = mantle_obj)

for (i in 2:5) {
  cat('\n')
  all_obj[[i]]@meta.data$tree.ident <- 
    factor(all_obj[[i]]@meta.data$tree.ident, ordered = TRUE)
}

for (i in 1:5) {
  all_obj[[i]]@assays$RNA@counts <- matrix()
  all_obj[[i]]@assays$RNA@scale.data <- matrix()
  saveRDS(all_obj[[i]], paste0("TRIMMED_", names(all_obj)[i], ".RDS"))
}


# =============================================================================
# Processing for individual data set that is recently added


files <- list.files(".", pattern = "*.RDS", full.names = TRUE)
new_obj <- readRDS(files[1])

new_obj@meta.data$tree.ident <- 
  factor(new_obj@meta.data$tree.ident, ordered = TRUE)

new_obj@assays$RNA@counts <- matrix()
new_obj@assays$RNA@scale.data <- matrix()
new_obj@assays$integrated@counts <- matrix()
new_obj@assays$integrated@scale.data <- matrix()
saveRDS(new_obj, paste0("TRIMMED_", "all_LL_cells_regen_anchored_v1.2_", ".RDS"))

for (i in seq_along(slotNames(new_obj))) {
  print(format(object.size(slot(new_obj, slotNames(new_obj)[i])), units = "Mb"))
}

print(format(object.size(new_obj@assays$RNA), units = "Mb"))
print(format(object.size(new_obj@assays$integrated), units = "Mb"))


# =============================================================================


# Changing resolution back to tree ident
files <- list.files("./data", pattern = "TRIMMED", full.names = TRUE)
all_data <- readRDS(files[1])
colnames(all_data@meta.data)
all_data@meta.data$tree.ident <- all_data@meta.data$tree.ident.res0.4
saveRDS(all_data, files[1])