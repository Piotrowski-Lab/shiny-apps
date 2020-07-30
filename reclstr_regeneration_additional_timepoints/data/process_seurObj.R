library(devtools)
#dev_mode(on=T)
#devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(hrbrthemes)
library(tidyr)
library(reshape2)

if (TRUE) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  script_name <- "process_seurObj"

  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}
  
}

files <- list.files(".", pattern = "TRIMMED", full.names = TRUE)
file_list <- list()

cell.type <- c("mature-HCs","young-HCs","HC-prog" ,"central-cells", "DV-cells","AP-cells",  
               "amp-SCs","mantle-cells", "Inm")
treatments <- c("homeo" ,"0min" , "30min", "1hr", "1.5hr", "2hr", "3hr","5hr", "10hr")

readSeuratObj <- TRUE
modifySeuratObj <-FALSE

for (i in 1:length(files)) {
  if (readSeuratObj){
    print("Loading Seurat objects...")
    file_list[[i]] <- readRDS(files[i])
    DefaultAssay(file_list[[i]]) <- "RNA"
  }
  if (modifySeuratObj){
    print('hi')
    #create new column in meta.data 
    #applied to all-she-pos and neuromast analysis
    if ("cell.type.ident" %in% colnames(file_list[[i]]@meta.data)){
      print('hello')
      file_list[[i]]@meta.data$cell.type.ident.by.data.set <- factor(paste(file_list[[i]]@meta.data$cell.type.ident,
                                                                           file_list[[i]]@meta.data$data.set,sep="_"))
      if ("early-HCs" %in% file_list[[i]]$cell.type.ident){
        file_list[[i]]@meta.data$cell.type.ident <- plyr::revalue(
          file_list[[i]]@meta.data$cell.type.ident, c("early-HCs" = "young-HCs"))
      }
    }
    if ("cell_type" %in% colnames(file_list[[i]]@meta.data)){
      print(files[i])
      file_list[[i]]@meta.data$cell.type.ident.by.data.set <- factor(paste(file_list[[i]]@meta.data$cell_type,
                                                                           file_list[[i]]@meta.data$data.set,sep="_"))
      if ("early-HCs" %in% file_list[[i]]$cell_type){
        file_list[[i]]@meta.data$cell_type <- plyr::revalue(
          file_list[[i]]@meta.data$cell_type, c("early-HCs" = "young-HCs"))
      }
    }
    file_list[[i]]$cell.type.ident.by.data.set <- factor(file_list[[i]]$cell.type.ident.by.data.set)
    
    # print("saving object...")
    # saveRDS(file_list[[i]], file = files[i])
  }
}

# ==================== recluster subsetted ident object list (omit neuromast analysis)
combined_obj <- file_list[2:6]

combined_obj <- readRDS("./CombinedObj_central_AP_mantle_seurat3_v1.0_.RDS")

combined_obj <- list(combined_obj$`AP-cells`, combined_obj$`central-cells`,
                     combined_obj$`DV-cells`, combined_obj$`HC-prog`,
                     combined_obj$`mantle-cells`)

combined_obj <- list("AP-cells" = combined_obj[[1]],
                "central-cells" = combined_obj[[2]],
                "DV-cells" = combined_obj[[3]],
                "HC-prog" = combined_obj[[4]],
                "mantle-cells" = combined_obj[[5]]) #rename list elements

for ( i in 1:length(combined_obj)){
  print(paste0("before: ", object.size(combined_obj[[i]]), units = "MB"))
  DefaultAssay(file_list[[i]]) <- "RNA"
  combined_obj[[i]] <- ScaleData(combined_obj[[i]], features = rownames(combined_obj[[i]]))
  combined_obj[[i]] <- RunPCA(combined_obj[[i]])
  combined_obj[[i]] <- FindNeighbors(combined_obj[[i]], dims = 1:10)
  combined_obj[[i]] <- FindClusters(combined_obj[[i]], resolution = 0.6)
  combined_obj[[i]] <- RunUMAP(combined_obj[[i]], dims = 1:10)
  print(paste0("after: ", object.size(combined_obj[[i]]), units = "MB"))
  
}

DimPlot(combined_obj$`central-cells`, group.by = "data.set")


#reorder cell.type.ident.by.data.set in all-she-pos analysis & neuromast analysis
for (i in 1:1){
  print(file_list[[i]])
  Idents(file_list[[i]]) <- "cell.type.ident.by.data.set"
  my_levels <- as.vector(t(outer(cell.type, treatments, paste, sep="_"))) 
  file_list[[i]]@active.ident <- factor(file_list[[i]]@active.ident, levels= my_levels)
  file_list[[i]]$cell.type.ident.by.data.set <- factor(file_list[[i]]$cell.type.ident.by.data.set, 
                                                       levels= my_levels)
  file_list[[i]]@active.ident <- droplevels(file_list[[i]]@active.ident)
  file_list[[i]]$cell.type.ident.by.data.set <- droplevels(file_list[[i]]$cell.type.ident.by.data.set)
  Idents(file_list[[i]]) <- "cell.type.ident"
  saveRDS(file_list[[i]], file = files[i])
  }

