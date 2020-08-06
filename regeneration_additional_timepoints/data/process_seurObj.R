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
library(stringr)

if (TRUE) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  script_name <- "process_seurObj"
}

files <- list.files(".", pattern = "TRIMMED", full.names = TRUE)
file_list <- list()

cell.type <- c("mature-HCs","young-HCs","HC-prog" ,"central-cells", "DV-cells","AP-cells",  
               "amp-SCs","mantle-cells", "Inm")
treatments <- c("homeo" ,"0min" , "30min", "1hr", "1.5hr", "2hr", "3hr","5hr", "10hr")

readSeuratObj <- TRUE
modifySeuratObj <-TRUE

for (i in 1:length(files)) {
  if (readSeuratObj){
    print("Loading Seurat objects...")
    file_list[[i]] <- readRDS(files[i])
    DefaultAssay(file_list[[i]]) <- "RNA"
  }
  if (modifySeuratObj){
    print('hi')
      if ("early-HCs" %in% file_list[[i]]$cell.type.ident){
        print('changing early-HC to young-HCs')
        file_list[[i]]@meta.data$cell.type.ident <- plyr::revalue(
          file_list[[i]]@meta.data$cell.type.ident, c("early-HCs" = "young-HCs"))
      }
    #create new column in meta.data 
    #applied to all-she-pos and neuromast analysis
        # if ("cell.type.ident" %in% colnames(file_list[[i]]@meta.data)){
        #    print('changing metadata colnames')
        #     file_list[[i]]@meta.data$cell.type.ident.by.data.set <- factor(paste(file_list[[i]]@meta.data$cell.type.ident,
        #                                                                    file_list[[i]]@meta.data$data.set,sep="_"))
        # }
      
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
    
    #print("saving object...")
   # saveRDS(file_list[[i]], file = files[i])
  }
}

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

# =============================== modify ptime regen objects
readSeuratObj <- TRUE
modifySeuratObj <-FALSE

files <- files[7:9]
file_list <- list()

for (i in 1:length(files)) {
  if (readSeuratObj){
    print(paste0("reading: ", files[i]))
    file_list[[i]] <- readRDS(files[i])
    DefaultAssay(file_list[[i]]) <- "RNA"
  }
  if(modifySeuratObj){
    
    if ("early-HCs" %in% file_list[[i]]$cell.type.ident){
      print('changing early-HC to young-HCs')
      file_list[[i]]@meta.data$cell.type.ident <- plyr::revalue(
        file_list[[i]]@meta.data$cell.type.ident, c("early-HCs" = "young-HCs"))
    }
    print('changing column name')
    #change cell.type.ident.by.data.set
    names(file_list[[i]]@meta.data)[names(file_list[[i]]@meta.data) == "cell.type.and.trt"] <- "cell.type.ident.by.data.set"
    file_list[[i]]@meta.data$cell.type.ident.by.data.set <- gsub(pattern = "early-HCs", "young-HCs", 
    file_list[[i]]@meta.data$cell.type.ident.by.data.set)   
    file_list[[i]]$cell.type.ident.by.data.set <- factor(file_list[[i]]$cell.type.ident.by.data.set)
    }
  # print(paste0("saving: ", files[i]))
  # saveRDS(file_list[[i]], file = files[i])
}

# =============================== reorder cell.type.ident, ptime objects
cell.type <- c("central-cells", "HC-prog","young-HCs","mature-HCs")
#reorder cell.type.ident.by.data.set in all-she-pos analysis & neuromast analysis
for (i in 1:3){
  print(file_list[[i]])
  Idents(file_list[[i]]) <- "cell.type.ident"
  file_list[[i]]@active.ident <- factor(file_list[[i]]@active.ident, levels= cell.type)
  file_list[[i]]$cell.type.ident <- factor(file_list[[i]]$cell.type.ident,
                                                       levels= cell.type)
  file_list[[i]]@active.ident <- droplevels(file_list[[i]]@active.ident)
  file_list[[i]]$cell.type.ident <- droplevels(file_list[[i]]$cell.type.ident)
  #drop levels for data.set
  file_list[[i]]$data.set <- droplevels(file_list[[i]]$data.set)
  
  print(paste0("saving: ", files[i]))
  saveRDS(file_list[[i]], file = files[i])
}

# ======= manual reorder of 

ptime_home_regen <- file_list[[1]]
type_treat_levels <- c( "homeo.central-cells", "0min.central-cells" , "30min.central-cells",
                        "1hr.central-cells" ,  "1.5hr.central-cells" ,"2hr.central-cells",
                        "3hr.central-cells" ,  "5hr.central-cells",   "1hr.HC-prog",
                        "1.5hr.HC-prog"  ,     "2hr.HC-prog"   ,      "3hr.HC-prog",
                        "5hr.HC-prog"    ,     "10hr.HC-prog"    ,    "homeo.HC-prog",
                        "1.5hr.young-HCs"   ,  "2hr.young-HCs"   ,    "3hr.young-HCs",
                        "5hr.young-HCs"   ,    "10hr.young-HCs"    ,  "2hr.mature-HCs", "3hr.mature-HCs"   ,
                        "5hr.mature-HCs"   ,   "10hr.mature-HCs",
                        "homeo.mature-HCs")
ptime_home_regen$cell.type.ident.by.data.set <- factor(ptime_home_regen$cell.type.ident.by.data.set, levels = type_treat_levels)
levels(ptime_home_regen$cell.type.ident.by.data.set)
levels(ptime_home_regen$cell.type.ident)
levels(ptime_home_regen$data.set)
saveRDS(ptime_home_regen, file = files[1])

ptime_regen <- file_list[[2]]
type_treat_levels <- c("0min.central-cells" , "30min.central-cells", "1hr.central-cells",
                       "1.5hr.central-cells" ,"2hr.central-cells" ,  "3hr.central-cells",
                       "5hr.central-cells"  , "1hr.HC-prog"      ,   "1.5hr.HC-prog",
                       "2hr.HC-prog" ,        "3hr.HC-prog"  ,       "5hr.HC-prog",
                       "10hr.HC-prog"   ,     "1.5hr.young-HCs"   ,  "2hr.young-HCs",
                       "3hr.young-HCs"   ,    "5hr.young-HCs"   ,    "10hr.young-HCs",
                       "2hr.mature-HCs"   ,   "3hr.mature-HCs"  ,    "5hr.mature-HCs",
                       "10hr.mature-HCs")
ptime_regen$cell.type.ident.by.data.set <- factor(ptime_regen$cell.type.ident.by.data.set, levels = type_treat_levels)
levels(ptime_regen$cell.type.ident.by.data.set)
levels(ptime_regen$cell.type.ident)
levels(ptime_regen$data.set)
saveRDS(ptime_regen, file = files[2])


ptime_homeo <- file_list[[3]]
type_treat_levels <- c( "homeo.central-cells" ,"homeo.HC-prog"    ,   "homeo.young-HCs",
                        "homeo.mature-HCs")
ptime_homeo$cell.type.ident.by.data.set <- factor(ptime_homeo$cell.type.ident.by.data.set, levels = type_treat_levels)
levels(ptime_homeo$cell.type.ident.by.data.set)
levels(ptime_homeo$cell.type.ident)
saveRDS(ptime_homeo, file = files[3])

