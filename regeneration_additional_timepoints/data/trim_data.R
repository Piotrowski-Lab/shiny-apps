library(Seurat)
files <- list.files(".", pattern = "*.RDS", full.names = TRUE)
combined_obj <- readRDS(files[1])

all_obj <- list("all_data_sets" = combined_obj)

for (i in 1) {
  all_obj[[i]]@assays$RNA@counts <- matrix()
  all_obj[[i]]@assays$RNA@scale.data <- matrix()
  saveRDS(all_obj[[i]], paste0("TRIMMED_", names(all_obj)[i], ".RDS"))
}
