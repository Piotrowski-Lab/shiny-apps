{library(devtools)
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
library(patchwork)

if (TRUE) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  script_name <- "process_seurObj"
}
files <- list.files(".", pattern = "TRIMMED", full.names = TRUE)
file_list <- list()
}

readSeuratObj <- TRUE

for (i in 1:length(files)) {
  if (readSeuratObj){
    print("Loading Seurat objects...")
    file_list[[i]] <- readRDS(files[i])
    DefaultAssay(file_list[[i]]) <- "RNA"
  }
}
seurat_obj <- file_list[[9]] #ptime homeo and regen combined

seurat_obj$data.set <- droplevels(seurat_obj$data.set)
ids <- as.list(levels(seurat_obj$data.set))

selected <- c("atoh1a", "her4.1")


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


obj_trt_list <- list()[1:length(ids)]

for (i in 1:length(ids)) {
  print(ids[[i]])
  obj_trt_list[[i]] <- seurat_obj[,seurat_obj[["data.set"]] == ids[[i]]]
}
names(obj_trt_list) <- ids



stacked_violin_plot <- function(goi, obj_trt_list){
  trt_plot_list <- list()[1:length(ids)]
  names(trt_plot_list) <- ids
  for (i in 1:length(ids)) {
    vln_obj <- VlnPlot(
      obj_trt_list[[i]], features = goi, pt.size = 0.0) +
      xlab("") + ylab(ids[i]) + ggtitle("") +
      theme(legend.position = "none", axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = rel(1), angle = 0),
            axis.text.y = element_text(size = rel(1)),
            plot.margin = unit(c(-1.0, 0.5, -1.0, 0.5), "cm")) 
            # #placeholders
    
    trt_plot_list[[i]] <- vln_obj
  }
  
  trt_plot_list[[length(trt_plot_list)]]<- trt_plot_list[[length(trt_plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line()) +
    scale_x_discrete(limits = c(factor(levels(seurat_obj$seurat_clusters))))
  # change the y-axis tick to only max value, treats ymax from each obj_trt_list independently
  ymaxs <- purrr::map_dbl(trt_plot_list, extract_max)
  #finds highest ymax, normalize
  ymaxs<- max(sapply(ymaxs, max))
  trt_plot_list <- purrr::map2(trt_plot_list, ymaxs, function(x, y) x +
                                 scale_y_continuous(breaks = c(y)) + expand_limits(y = y))
  grid_obj <- cowplot::plot_grid(plotlist = trt_plot_list,
                                 nrow = length(ids), ncol = 1, axis = "l", align = "hv") + #rel_heights  = c(1,1,1,1,1,1)) +
    #theme(plot.margin = margin(2.0, 2.0, 2.0, 2.0, unit = "in")) 
    ggtitle(goi) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  
  return(grid_obj)
  
}
grid_obj <- list()[1:length(selected)]
for (i in 1:length(selected)){
  print(selected[[i]])
  #options(repr.plot.width = 12, repr.plot.height = 3)
  
  grid_obj[[i]] <- stacked_violin_plot(goi = selected[[i]], obj_trt_list = obj_trt_list)
  
}
names(grid_obj) <- selected

final_grid <- cowplot::plot_grid(plotlist = grid_obj, nrow = length(grid_obj), axis = "l", align = "hv", scale = 0.9) +
  theme(plot.margin = margin(.2, .2, .2, .2, unit = "in"))
final_grid

# ======================== grouped heat map

features <- c("atoh1a", "her4.1", "hes2.2", "dld", "sox4a*1", "myclb", "gadd45gb.1",
              "insm1a", "wnt2", "sost", "sfrp1a", "pcna", "mki67", "isl1", "slc1a3a", "glula", "lfng", "cbln20", "ebf3a",
              "znf185", "si:ch211-229d2.5", "si:ch73-261i21.5", "spaca4l", "foxp4", "crip1")


dotplot <- DotPlot(seurat_obj, features = features,
                   group.by = "cell.type.ident.by.data.set")

# dotplot$data$groupIdent <- gsub("(.+?)(\\_.*)", "\\1",dotplot$data$id)
# dotplot$data$groupIdent <- factor(dotplot$data$groupIdent,levels=levels(seurat_obj$cell.type.ident))

if(TRUE){
  dotplot$data$groupIdent <- gsub("^.*\\.", "",dotplot$data$id)
  dotplot$data$groupIdent <- factor(dotplot$data$groupIdent,levels=levels(seurat_obj$cell.type.ident))
}

g <- ggplot(dotplot$data, aes(id, features.plot,fill= avg.exp.scaled, width = 1, height = 1)) + 
  geom_tile() +
  scale_fill_distiller(
    palette = "RdYlBu") +
  theme_ipsum()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5,size = 13),
        axis.title.y.right = element_text(size=13),panel.spacing = unit(.35, "lines")) + 
  facet_grid( ~ groupIdent, scales='free_x') +
  ylim(rev(levels(dotplot$data$features.plot))) 

g

# ========================== stkd violin plot by gene
group.by <- "data.set"

p_size <- as.numeric(.5)

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = p_size, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          grouped = group.by) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = grouped )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, grouped = group.by, pt.size = p_size))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  #finds highest ymax, normalize
  ymaxs<- max(sapply(ymaxs, max))
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

StackedVlnPlot(obj = seurat_obj, features = features)
