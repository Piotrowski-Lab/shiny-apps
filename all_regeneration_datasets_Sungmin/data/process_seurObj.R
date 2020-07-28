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
  
  script_name <- "modify_seuratObj"
  
  date <-Sys.Date()
  
  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}
  
  }

files <- list.files(".", pattern = "TRIMMED", full.names = TRUE)
file_list <- list()

cell.type <- c("mature-HCs","young-HCs","HC-prog" ,"central-cells", "DV-cells","AP-cells",  
               "amp-SCs","mantle-cells", "Inm","blood", "spi1b-pos","krt17-pos","twist3-pos","tm4sf4-pos")
treatments <- c("homeo" ,"0min" , "30min", "1hr", "3hr","5hr", "10hr")

readSeuratObj <- TRUE
modifySeuratObj <-FALSE

for (i in 1:length(files)) {
  if (readSeuratObj){
  print("Loading Seurat objects...")
  file_list[[i]] <- readRDS(files[i])
  DefaultAssay(file_list[[i]]) <- "RNA"
  }
  if (modifySeuratObj){
  #create new column in meta.data 
    #applied to all-she-pos and neuromast analysis
  if ("cell.type.ident" %in% colnames(file_list[[i]]@meta.data)){
    file_list[[i]]@meta.data$cell.type.ident <- plyr::revalue(
      file_list[[i]]@meta.data$cell.type.ident, c("early-HCs" = "young-HCs"))
  file_list[[i]]@meta.data$cell.type.ident.by.data.set <- factor(paste(file_list[[i]]@meta.data$cell.type.ident,
                                                          file_list[[i]]@meta.data$data.set,sep="_"))
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
    
  saveRDS(file_list[[i]], file = files[i])
  }
}

  #reorder cell.type.ident.by.data.set in all-she-pos analysis & neuromast analysis
for (i in 6:5){
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

# ================== test new heatmap ======================
features <- c("atoh1a", "her4.1", "hes2.2", "dld", "sox4a*1", "myclb", "gadd45gb.1",
              "insm1a", "wnt2", "sost", "sfrp1a", "pcna", "mki67", "isl1", "slc1a3a", "glula", "lfng", "cbln20", "ebf3a",
              "znf185", "si:ch211-229d2.5", "si:ch73-261i21.5", "spaca4l", "foxp4", "crip1")

seurat_obj <- file_list[[6]]

dotplot <- DotPlot(seurat_obj, features = features,
                   group.by = "cell.type.ident.by.data.set")

dotplot$data$groupIdent <- gsub("(.+?)(\\_.*)", "\\1",dotplot$data$id)
dotplot$data$groupIdent <- factor(dotplot$data$groupIdent,levels=levels(seurat_obj$cell.type.ident))

g <- ggplot(dotplot$data, aes(id, features.plot,fill= avg.exp.scaled, width = 1, height = 1)) + 
  geom_tile() +
  scale_fill_distiller(
    palette = "RdYlBu") +
  theme_ipsum()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5,size = 13),
        axis.title.y.right = element_text(size=13),panel.spacing = unit(.35, "lines")) + facet_grid( ~ groupIdent, scales='free_x')

g

# =================== Individual Cell Heatmap =====================

#from DoHeatMap

"%||%" <- devtools:::`%||%`

group.by <- "cell.type.ident.by.data.set"
cells <- NULL
col.min = -2.5
col.max = 2.5

#object <- suppressMessages(expr = StashIdent(object = seurat_obj, save.name = 'ident'))

cells <- cells %||% colnames(x = seurat_obj)

data <- as.data.frame(x = t(x = as.matrix(x = GetAssayData(
  object = seurat_obj, slot = "data")[features, cells, drop = FALSE])))


data <- scale(data)
data <- as.data.frame(MinMax(data = data, min = col.min, max = col.max))

data$id <- if (is.null(x = group.by)) {
  Idents(object = seurat_obj)[cells, drop = TRUE]
} else {
  seurat_obj[[group.by, drop = TRUE]][cells, drop = TRUE]
}
if (!is.factor(x = data$id)) {
  data$id <- factor(x = data$id)
}
data$id <- as.vector(x = data$id)

data$Cell <- rownames(data)
data <- melt(data, variable.name  = "Feature")
# data$groupIdent <- gsub("(.+?)(\\_.*)", "\\1",data$id)
# data$groupIdent <- factor(data$groupIdent,levels=cell.type)
#preserve identity order
#group.by.f <- factor(group.by)
if (group.by == "cell.type.ident.by.data.set"){
  print('hi')
data$id <- factor(data$id, levels = levels(seurat_obj$cell.type.ident.by.data.set))
}else if (group.by == "data.set"){
  data$id <- factor(data$id, levels = levels(seurat_obj$data.set))
}else{
  data$id <- factor(data$id, levels = levels(seurat_obj$cell.type.ident))
}

indv.hmap <- ggplot(data, aes(Cell, Feature,fill= value, width = 1, height = 1)) + 
  geom_tile() +
  scale_fill_distiller(
    palette = "RdYlBu") +
  theme_ipsum()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y.right = element_text(size=13),panel.spacing = unit(.25, "lines"),
        strip.text.x  = element_text(angle = 90, vjust = 0.5, hjust=.5,size = 8)) + facet_grid( ~ id, space = 'free', scales = 'free')

png("IndvHeatmap.png"
    ,width = 30, height = 10, units = "in", res = 300)
print(indv.hmap)
dev.off()

View(data)

# ======================================= Downsample ===========================
percentage <- as.numeric(c(.75,.50,.25))
object_list <- list()[1:length(percentage)]
for(i in 1:length(object_list)){
  object_list[[i]] = subset(seurat_obj, cells = sample(Cells(seurat_obj), round(percentage[i]*length(colnames(seurat_obj)))))
  #object_list[[i]] <- seurat_obj[, sample(Cells(seurat_obj), size = round(percentage[i]*length(colnames(seurat_obj))), replace=F)]
  
}
table(seurat_obj$cell.type.ident)
table(object_list[[1]]$cell.type.ident)
