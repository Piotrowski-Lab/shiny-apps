library(Seurat)
library(ggplot2)

if (TRUE) {
  #setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  script_name <- "stacked_violin"
  
  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}
  
}
# ======================================== ==============================================


seurat_obj <- readRDS("./subsampled_30_smrtseq_10X.RDS")

ids <- as.list(levels(seurat_obj$data.set))


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
#})

trt_plot_list <- list()[1:length(ids)]
names(trt_plot_list) <- ids
stacked_violin_plot <- function(goi, obj_trt_list){
  trt_plot_list <- list()[1:length(ids)]
  names(trt_plot_list) <- ids
  for (i in 1:length(ids)) {
    vln_obj <- VlnPlot(
      obj_trt_list[[i]], features = goi, pt.size = 0) +
      xlab("") + ylab(ids[i]) + ggtitle("") +
      theme(legend.position = "none", axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = rel(1), angle = 0),
            axis.text.y = element_text(size = rel(1)),
            plot.margin = unit(c(-1, 0.5, -1, 0.5), "cm"))
    trt_plot_list[[i]] <- vln_obj
  }
  
  trt_plot_list[[length(trt_plot_list)]]<- trt_plot_list[[length(trt_plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  # change the y-axis tick to only max value
  ymaxs <- purrr::map_dbl(trt_plot_list, extract_max)
  ymaxs<- max(sapply(ymaxs, max))
  trt_plot_list <- purrr::map2(trt_plot_list, ymaxs, function(x, y) x +
                                 scale_y_continuous(breaks = c(y)) + expand_limits(y = y))
  grid_obj <- cowplot::plot_grid(plotlist = trt_plot_list,
                                 nrow = length(ids), ncol = 1, axis = "l", align = "hv", rel_heights  = c(3,3,3,3,3,3)) +
    #theme(plot.margin = margin(1.0, 1.0, 1.0, 1.0, unit = "in")) +
    ggtitle(goi) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  

  return(grid_obj)
  
}

selected <- c("atoh1a","wnt2")
selected <- unlist(strsplit(selected, " "))

grid_obj <- list()[1:length(selected)]
for (i in 1:length(selected)){
  print(selected[[i]])
  options(repr.plot.width = 12, repr.plot.height = 3)
  grid_obj[[i]] <- stacked_violin_plot(goi = selected[[i]], obj_trt_list = obj_trt_list)
  
}
names(grid_obj) <- selected

final_grid <- cowplot::plot_grid(plotlist = grid_obj, nrow = length(grid_obj), axis = "l", align = "hv", scale = 0.9) +
  theme(plot.margin = margin(.2, .2, .2, .2, unit = "in"))

final_grid

png(figurePath(paste0("stacked_violin.png"))
    ,width = 11, height = 9, units = "in", res = 300)
print(final_grid)
dev.off()
