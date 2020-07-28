library(Seurat)
library(dplyr)
library(kable)
library(kableExtra)
library(writexl)


if (TRUE) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  script_name <- "process_seurObj"
  app_name <- "smrtseq_vs_10X_scRNAseq"
  
  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}
}

files <- list.files("./seurat_obj_input", pattern = "*scaled", full.names = TRUE)

seurat_obj <- readRDS(file = files[1])

# ======================================= Modify RNA assay ===========================
print(object.size(seurat_obj), units = "MB")

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
print(object.size(seurat_obj), units = "MB")

#change cell.type.idents
seurat_obj@meta.data$cell.type.ident <- plyr::revalue(
  seurat_obj@meta.data$cell.type.ident, c("early-HCs" = "young-HCs"))

#reoder data.set
#temp <- seurat_obj
Idents(seurat_obj) <- "data.set"
my_levels <- c("1hr-smrtseq", "homeo-smrtseq", "homeo-10X-isl1", "homeo-10X-2047", "homeo-10X-2410-7", "homeo-10X-2410-8")
seurat_obj$data.set <- factor(Idents(seurat_obj), levels= my_levels)

seurat_obj@assays$RNA@counts <- matrix()
seurat_obj@assays$RNA@scale.data <- matrix()
seurat_obj[["integrated"]]@scale.data <- matrix() 

print(object.size(seurat_obj), units = "MB")

saveRDS(seurat_obj, paste0("TRIMMED_", app_name, ".RDS"))

# =======================================Downsample ===========================
downsample <- FALSE
if (downsample){
Idents(seurat_obj) <- "seq.method"
smartseq <- WhichCells(seurat_obj, idents = "smartseq2")
tenX <- WhichCells(seurat_obj, idents = "10X", downsample = round(length(colnames(seurat_obj)) * .30))
seurat_obj <- subset(seurat_obj, cells = c(smartseq, tenX))
print(object.size(seurat_obj), units = "MB")

#revert back to default idents
Idents(seurat_obj) <- "cell.type.ident"

addmargins(table(Idents(seurat_obj), seurat_obj$data.set))

saveRDS(seurat_obj, "./subsampled_30_smrtseq_10X.RDS")


seurat_obj <- readRDS("./subsampled_30_smrtseq_10X.RDS")
seurat_obj
library(kableExtra)
kable(table(seurat_obj$data.set),col.names = c("data.set", "freq")) %>%kable_styling(font_size = 10) %>% row_spec(1:2, color = "red")
}
# ======================================= Generate Figures ===========================

diff_results <- FindAllMarkers(seurat_obj)
diff_results["Gene.name.uniq"] <- row.names(diff_results)

gene_table <- read.table("./Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)
diff_results <- merge(diff_results,gene_table, by = "Gene.name.uniq", sort = FALSE)

mcPlotMaker <- function(object, path ){
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
      
      ifelse(!dir.exists(file.path(path, "cluster-top-markers/")), dir.create(file.path(path, "cluster-top-markers/")), FALSE)
      ifelse(!dir.exists(file.path(path, "VlnPlot-cluster-top-markers/")), dir.create(file.path(path, "VlnPlot-cluster-top-markers/")), FALSE)
      
      img_path <- paste0(path, "cluster-top-markers/", "cluster_",idents[x], "_top_",
        seq_nums[j], "-", (seq_nums[j] + 19), "_features.png")
      png(img_path, width = 30, height = 25, units = "in", res = 200)
      print(cowplot::plot_grid(plotlist = f))
      dev.off()
      img_path <- paste0(path, "VlnPlot-cluster-top-markers/", "cluster_",idents[x], "_top_",
        seq_nums[j], "-", (seq_nums[j] + 19), "_features.png")
      png(img_path, width = 30, height = 25, units = "in", res = 200)
      print(v)
      dev.off()
    }
    
  }
  names(cluster_list) <- paste("cluster_",idents)
  names(top100_list) <- paste("cluster_",idents)
  return(cluster_list)
}
cluster_list <- capture.output(mcPlotMaker(object = seurat_obj))

# ======================================= Generate Figures For Smartseq2 cells ===========================
smartseq2 <- TRUE
tenX <- TRUE
if (smartseq2){
  print("subsetting")
  seurat_obj <- subset(seurat_obj, subset = seq.method == "smartseq2")
  print("running mcPlotMaker")
  cluster_list <- mcPlotMaker(object = seurat_obj, path = "./process_seurObj_figures/smartseq2_analysis/")
  # ============================== smartseq2 heatmap =======================
  
  top100 <- c(cluster_list$`cluster_ central-cells`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ amp-SCs`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ AP-cells`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ DV-cells`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ HC-prog`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ mantle-cells`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ mature-HCs`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ Inm`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ young-HCs`$Gene.name.uniq[1:100])
  png(figurePath(paste0("smartseq2_analysis/heatmap-cluster-top-markers/heatmap.allmarkers.png"))
      ,width = 30, height = 100, units = "in", res = 300)
  print(DoHeatmap(seurat_obj, features = top100) + NoLegend())
  dev.off()
  
  ifelse(!dir.exists(file.path("./process_seurObj_figures/smartseq2_analysis/", "gene_tabl_by_cluster/")), dir.create(file.path("./process_seurObj_figures/smartseq2_analysis/", "gene_tabl_by_cluster/")), FALSE)
  write_xlsx(list("central-cells" = cluster_list$`cluster_ central-cells`,
                    "amp-SCs" = cluster_list$`cluster_ amp-SCs`,
                  "AP-cells" = cluster_list$`cluster_ AP-cells`,
                  "DV-cells" = cluster_list$`cluster_ DV-cells`,
                  "HC-prog" = cluster_list$`cluster_ HC-prog`,
                  "mantle-cells" = cluster_list$`cluster_ mantle-cells`,
                  "mature-HCs" = cluster_list$`cluster_ mature-HCs`,
                  "Inm" = cluster_list$`cluster_ Inm`,
                  "early-HCs" = cluster_list$`cluster_ young-HCs`),
             path = "./process_seurObj_figures/smartseq2_analysis/gene_tabl_by_cluster/gene_table_by_cluster.xlsx")
  
}

#reset seurat_obj
seurat_obj <- readRDS("./subsampled_30_smrtseq_10X.RDS")

if(tenX){
  print("subsetting")
  seurat_obj <- subset(seurat_obj, subset = seq.method == "10X")
  print("running mcPlotMaker")
  cluster_list <- mcPlotMaker(object = seurat_obj, path = "./process_seurObj_figures/10X_analysis/")
  # ============================== 10X heatmap =======================
  
  top100 <- c(cluster_list$`cluster_ central-cells`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ amp-SCs`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ AP-cells`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ DV-cells`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ HC-prog`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ mantle-cells`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ mature-HCs`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ Inm`$Gene.name.uniq[1:100],
              cluster_list$`cluster_ young-HCs`$Gene.name.uniq[1:100])
  png(figurePath(paste0("10X_analysis/heatmap.allmarkers.png"))
      ,width = 30, height = 100, units = "in", res = 300)
  print(DoHeatmap(seurat_obj, features = top100) + NoLegend())
  dev.off()
  
  ifelse(!dir.exists(file.path("./process_seurObj_figures/10X_analysis/", "gene_tabl_by_cluster/")), dir.create(file.path("./process_seurObj_figures/10X_analysis/", "gene_tabl_by_cluster/")), FALSE)
  write_xlsx(list("central-cells" = cluster_list$`cluster_ central-cells`,
                  "amp-SCs" = cluster_list$`cluster_ amp-SCs`,
                  "AP-cells" = cluster_list$`cluster_ AP-cells`,
                  "DV-cells" = cluster_list$`cluster_ DV-cells`,
                  "HC-prog" = cluster_list$`cluster_ HC-prog`,
                  "mantle-cells" = cluster_list$`cluster_ mantle-cells`,
                  "mature-HCs" = cluster_list$`cluster_ mature-HCs`,
                  "Inm" = cluster_list$`cluster_ Inm`,
                  "early-HCs" = cluster_list$`cluster_ young-HCs`),
             path = "./process_seurObj_figures/10X_analysis/gene_tabl_by_cluster/gene_table_by_cluster.xlsx")
  
}

