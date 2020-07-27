library(Seurat)
library(dplyr)

if (TRUE) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

  script_name <- "Imn_mantle_homeo"
  
  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}
}

seurat_obj <- readRDS("./Imn_mantle_homeo.RDS")

DefaultAssay(seurat_obj) <- "RNA"

diff_table <- FindAllMarkers(seurat_obj)
diff_table["Gene.name.uniq"] <- row.names(diff_table)

gene_table <- read.table("./Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)
diff_results <- merge(diff_table,gene_table, by = "Gene.name.uniq", sort = FALSE)

x <- diff_results %>% group_by(as.factor(cluster)) %>% top_n(n = 100, wt = avg_logFC)

idents <- as.character(unique(Idents(seurat_obj)))
cluster_list <- vector(mode = "list", length = length(idents))
top100_list <- vector(mode = "list", length = length(idents))
seq_nums <- seq(1, 100, by = 20)

for (x in 1:length(idents)){
  print(paste0("finding markers for cluster: ", idents[x]))
  cluster_list[[x]] <- FindMarkers(seurat_obj, ident.1 = idents[x], only.pos = TRUE, verbose = TRUE)
  #print(cluster_list[[x]])
  cluster_list[[x]]["Gene.name.uniq"] <- row.names(cluster_list[[x]])
  cluster_list[[x]] <- merge(cluster_list[[x]], gene_table, by = "Gene.name.uniq", sort = FALSE)
  #sort by desc order of avg-logFC
  cluster_list[[x]] <- cluster_list[[x]][order( cluster_list[[x]][,3], decreasing = TRUE),]
  print(paste0("finding top 100 genes for cluster: ", idents[x]))
  top100_list[[x]] <- cluster_list[[x]] %>% top_n(n = 100, wt = avg_logFC)
  top100_list[[x]] <- top100_list[[x]]$Gene.name.uniq 
  for(j in 1:length(seq_nums)) {
    to_plot <- top100_list[[x]][seq_nums[j]:(seq_nums[j]+19)]
    if (NA %in% to_plot) {break}
    f <- FeaturePlot(seurat_obj, to_plot,
                     reduction = "umap", pt.size = 0.25, combine = FALSE)
    
    v <- VlnPlot(seurat_obj, features = to_plot) 
    for(k in 1:length(f)) {
      f[[k]] <- f[[k]] + NoLegend() + NoAxes()
      
    }
  
    path <- figurePath(paste0(
      "cluster-", "top-markers/", "cluster_",idents[x], "_top_",
      seq_nums[j], "-", (seq_nums[j] + 19), "_features.png"))
    png(path, width = 30, height = 25, units = "in", res = 200)
    print(cowplot::plot_grid(plotlist = f))
    dev.off()
    path <- figurePath(paste0(
      "VlnPlot-cluster-", "top-markers/", "cluster_",idents[x], "_top_",
      seq_nums[j], "-", (seq_nums[j] + 19), "_features.png"))
    png(path, width = 30, height = 25, units = "in", res = 200)
    print(v)
    dev.off()
  }
  
}


names(cluster_list) <- paste("cluster_",idents)
names(top100_list) <- paste("cluster_",idents)

# ============================== heatmap =======================

top100 <- c(top100_list$`cluster_ 0`, top100_list$`cluster_ 1`, top100_list$`cluster_ 2`, top100_list$`cluster_ 3`,
            top100_list$`cluster_ 4`, top100_list$`cluster_ 5`, top100_list$`cluster_ 6`,
            top100_list$`cluster_ 7`, top100_list$`cluster_ 8`, top100_list$`cluster_ 9`)
png(figurePath(paste0("./heatmap-cluster-top-markers/heatmap.allmarkers.png"))
    ,width = 30, height = 100, units = "in", res = 300)
print(DoHeatmap(seurat_obj, features = top100) + NoLegend())
dev.off()


# ========================= excel sheet ======================
library(writexl)

write_xlsx(list("cluster 0" = cluster_list$`cluster_ 0`,
                "cluster 1" = cluster_list$`cluster_ 1`,
                "cluster 2" = cluster_list$`cluster_ 2`,
                "cluster 3" = cluster_list$`cluster_ 3`,
                "cluster 4" = cluster_list$`cluster_ 4`,
                "cluster 5" = cluster_list$`cluster_ 5`,
                "cluster 6" = cluster_list$`cluster_ 6`,
                "cluster 7" = cluster_list$`cluster_ 7`,
                "cluster 8" = cluster_list$`cluster_ 8`,
                "cluster 9" = cluster_list$`cluster_ 9`),
           path = "./Imn_mantle_homeo_figures/gene_tabl_by_cluster/gene_table_by_cluster.xlsx")









