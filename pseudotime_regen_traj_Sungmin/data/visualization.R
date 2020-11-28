library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
devtools::load_all("/n/projects/nt2473/Analysis/Scripts/SeuratExtensions")
devtools::load_all("/n/projects/nt2473/Analysis/Scripts/CellTrajectoryExtensions")


dataPath <- function(filename){paste0("/n/projects/",
                                      "nt2473/Analysis/Data/sb2191-regen/", 
                                      filename)}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

seurat_obj <- readRDS(dataPath(paste0(
  "SeurObj_ptime_subset_experiments_seurat3_10hr_central_v1.0_dim_6_.RDS")))

cds <- readRDS(dataPath(paste0("saved_obj_", 
              "HC_lineage_homeo_10hr_monocle3_10_central_v1.0", "_.RDS")))

cell_type_trt <- levels(cds$cell.type.and.trt)
type_trt_cols <- gg_color_hue(length(cell_type_trt))

traject_plot <- plot_cells(cds, color_cells_by = "cell.type.and.trt",
                           label_leaves = TRUE, 
                           show_trajectory_graph = TRUE, 
                           cell_size = 0.70,
                           trajectory_graph_segment_size = 0.55, 
                           label_cell_groups = FALSE,
                           group_label_size = 8) + 
  theme(legend.position="bottom")

traject_plot <- cleanUMAP(traject_plot)
traject_plot <- traject_plot + scale_color_manual(values = type_trt_cols)

traject_plot

ptime_plot <- plot_cells(cds, color_cells_by = "pseudotime",
                         label_leaves = TRUE, 
                         show_trajectory_graph = TRUE, 
                         cell_size = 0.80,
                         trajectory_graph_segment_size = 0.60, 
                         label_cell_groups = FALSE)

ptime_plot <- ptime_plot +
  guides(fill = guide_legend(ncol = 1))
ptime_plot <- cleanUMAP(ptime_plot)
ptime_plot
# ============================================ overlap traj line with UMAP coordinate
x <- 1
y <- 2
ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>%
  as.data.frame() %>%
  dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
  dplyr::mutate(sample_name = rownames(.),
                sample_state = rownames(.))

dp_mst <- cds@principal_graph[["UMAP"]]

edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select_(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       source="sample_name",
                       source_prin_graph_dim_1="prin_graph_dim_1",
                       source_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       target="sample_name",
                       target_prin_graph_dim_1="prin_graph_dim_1",
                       target_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "target")

mst_branch_nodes <- branch_nodes(cds, reduction_method ="UMAP")
branch_point_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
  dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))

mst_leaf_nodes <- leaf_nodes(cds, reduction_method = "UMAP")
leaf_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
  dplyr::mutate(leaf_idx = seq_len(dplyr::n()))

mst_root_nodes <- root_nodes(cds, reduction_method = "UMAP")
root_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
  dplyr::mutate(root_idx = seq_len(dplyr::n()))



features <- c("atoh1a", "her4.1")

trajectory_graph_segment_size=0.75
graph_label_size=2

g <- FeaturePlot(seurat_obj, features = features, slot = "data")
g <- cleanUMAP(g)
g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1",
                                 y="source_prin_graph_dim_2",
                                 xend="target_prin_graph_dim_1",
                                 yend="target_prin_graph_dim_2"),
                      linetype="solid",
                      na.rm=TRUE,
                      data=edge_df) 
g <- g + #plot branching
  geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
             shape = 21, stroke=I(0.75),
             color="white",
             fill="black",
             size=I(2 * 1.5),
             na.rm=TRUE, branch_point_df) +
  geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                       label="branch_point_idx"),
            size=I(2), color="white", na.rm=TRUE,
            branch_point_df)

g <- g + #plot leaves
  geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
             shape = 21, stroke=I(0.75),
             color="black",
             fill="lightgray",
             size=I(2 * 1.5),
             na.rm=TRUE,
             leaf_df) +
  geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                       label="leaf_idx"),
            size=I(2), color="black", na.rm=TRUE, leaf_df)

g <- g + #plot root node
  geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
             shape = 21, stroke=I(0.75),
             color="black",
             fill="white",
             size=I(2 * 1.5),
             na.rm=TRUE,
             root_df) +
  geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                       label="root_idx"),
            size=I(2), color="black", na.rm=TRUE, root_df)
g

# ======================================================== line dynamic individual genes
gene <- c("atoh1a", "her4.1", "sox4a", "isl1")

plot_cells(ptime_central_traj)
plot_cells(ptime_main_traj)

make_plot_df <- function(cds_sub, gene){
  cds_sub <- cds_sub[rownames(cds_sub) %in% gene,]
  cds_exprs <- as.matrix(SingleCellExperiment::counts(cds_sub))
  
  count_mtx_sub <- cds_exprs
  #count_mtx_sub <- count_mtx_sub[!apply(count_mtx_sub,1,sum)==0,]
  count_mtx_sub <- t(count_mtx_sub)
  count_mtx_sub <- scale(log1p(count_mtx_sub))
  count_mtx_sub <- as.data.frame(Seurat::MinMax(count_mtx_sub , min = -2.5, 
                                                max = 2.5))
  mod1_df <- as.data.frame(count_mtx_sub)
  mod1_df$Cell <- rownames(mod1_df)
  
  mod1_df <- reshape2::melt(mod1_df)
  colnames(mod1_df)[2:3] <- c("Gene.name.uniq","expression")
  
  
  ptime_df <- data.frame(pseudotime = t(pseudotime(cds))[1,],
                         Cell = names(pseudotime(cds)),
                         cell_group = colData(cds)$cell.type.and.trt)
  
  #order cells in correspondence to ptime
  mod1_df$Cell <- factor(mod1_df$Cell, levels = ptime_df$Cell)
  
  plot_dt <- inner_join(mod1_df, ptime_df)
  
  
  return(plot_dt)
}

central_traj_df <- make_plot_df(cds_sub = ptime_central_traj,gene = gene)
nrow(central_traj_df)
main_HC_traj_df <- make_plot_df(cds_sub = ptime_main_traj, gene = gene)
nrow(main_HC_traj_df)

central_traj_geom_smooth_col <- "#00BE67"

main_traj_geom_smooth_col <- "#F8766D"

# color_indx <- which(cell_type_trt %in% levels(cds@colData$cell.type.and.trt))
# select_colors <- type_trt_cols[color_indx]
cell_type_trt <- levels(seurat_obj$cell.type.and.trt)
type_trt_cols <- gg_color_hue(length(cell_type_trt))

test_lg <- ggplot(central_traj_df) +
  geom_smooth(
    colour=central_traj_geom_smooth_col,
    span=0.2,
    method='loess',
    #fill = central_traj_geom_smooth_col,
    aes( x=pseudotime, y=expression,fill="Central Cell Lineage"),
    fullrange = TRUE)

test_lg <- test_lg +
  geom_smooth(data = main_HC_traj_df,
              colour=main_traj_geom_smooth_col,
              span=0.2,
              method='loess',
              #fill = main_traj_geom_smooth_col,
              aes( x=pseudotime, y=expression, fill="HC Lineage"))

test_lg <- test_lg + scale_fill_manual(name="legend",guide = 'legend',
                                       values = c("Central Cell Lineage" = central_traj_geom_smooth_col,
                                                  "HC Lineage" = main_traj_geom_smooth_col))
test_lg
             # labels = c("test1", "test2"))
              #values=c(central_traj_geom_smooth_col, main_traj_geom_smooth_col))

test_lg <- test_lg + geom_rug(data=main_HC_traj_df, sides='b', 
                              alpha=.10, aes(x=pseudotime,
                                             color = cell_group) )

test_lg <- test_lg + geom_rug(data=central_traj_df, sides='t', alpha=.10, 
                              aes(x=pseudotime,
                                 color = cell_group) )

test_lg <- test_lg + scale_color_manual(values = type_trt_cols)
test_lg <-   test_lg + theme_bw() +
  theme(
  ) +
  labs(
    x     = 'pseudotime'
    ,y    = 'z-scored gene expression'
  )
test_lg <- test_lg +theme(legend.position="bottom", legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

#add facet title to 
test_lg$data$title <- unique(test_lg$data$Gene.name.uniq)
test_lg <- test_lg + facet_wrap(~title)
test_lg <- test_lg+ 
  theme(strip.text.x = element_text(size = 18))
test_lg

#==========================

test_lg <- ggplot(central_traj_df) +
  geom_smooth(
    #colour=central_traj_geom_smooth_col,
    span=0.2,
    method='loess',
    #fill = central_traj_geom_smooth_col,
    aes( x=pseudotime, y=expression,color="Central Cell Lineage", fill = "Central Cell Lineage"),
    fullrange = TRUE)

test_lg <- test_lg +
  geom_smooth(data = main_HC_traj_df,
              #colour=main_traj_geom_smooth_col,
              span=0.2,
              method='loess',
              #fill = main_traj_geom_smooth_col,
              aes( x=pseudotime, y=expression, color="HC Lineage", fill = "HC Lineage"))

test_lg <- test_lg + scale_color_manual(name="Branching Trajectories",guide = 'legend',
                                       values = c("Central Cell Lineage" = central_traj_geom_smooth_col,
                                                  "HC Lineage" = main_traj_geom_smooth_col)) + guides(guide_legend(override.aes = list(linetype = c("black","black"))))   
test_lg <- test_lg + scale_fill_manual(name="Branching Trajectories",guide = 'legend',
                                        values = c("Central Cell Lineage" = central_traj_geom_smooth_col,
                                                   "HC Lineage" = main_traj_geom_smooth_col))
test_lg <- test_lg + theme(legend.direction = "horizontal",legend.position = "bottom")

test_lg <- test_lg +
  guides(colour = guide_legend(title.position = "top",title.hjust = 0.5))



test_lg

p1 <- gtable::gtable_filter(ggplot_gtable(ggplot_build(test_lg)), "guide-box") 

test_lg2 <- ggplot(central_traj_df) +
  geom_smooth(
    colour=central_traj_geom_smooth_col,
    span=0.2,
    method='loess',
    fill = central_traj_geom_smooth_col,
    aes( x=pseudotime, y=expression,color="Central Cell Lineage", fill = "Central Cell Lineage"),
    fullrange = TRUE)

test_lg2 <- test_lg2 +
  geom_smooth(data = main_HC_traj_df,
              colour=main_traj_geom_smooth_col,
              span=0.2,
              method='loess',
              fill = main_traj_geom_smooth_col,
              aes( x=pseudotime, y=expression, color="HC Lineage", fill = "HC Lineage"))


test_lg2 <- test_lg2 + geom_rug(data=main_HC_traj_df, sides='b', 
                              alpha=.10, aes(x=pseudotime,
                                             color = cell_group) )

test_lg2 <- test_lg2 + geom_rug(data=central_traj_df, sides='t', alpha=.10, 
                              aes(x=pseudotime,
                                  color = cell_group) )

test_lg2 <- test_lg2 + scale_color_manual(values = type_trt_cols)
test_lg2 <-   test_lg2 + theme_bw() +
  theme(
  ) +
  labs(
    x     = 'pseudotime'
    ,y    = 'z-scored gene expression'
  )
test_lg2 <- test_lg2 +theme(legend.position="bottom", legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

#add facet title to 
test_lg2$data$title <- unique(test_lg2$data$Gene.name.uniq)
test_lg2 <- test_lg2 + facet_wrap(~title)
test_lg2 <- test_lg2+ 
  theme(strip.text.x = element_text(size = 18))

test_lg2

p2 <- gtable::gtable_filter(ggplot_gtable(ggplot_build(test_lg2)), "guide-box") 
p2

#no legend plot
p3 <- ggplot(central_traj_df) +
  geom_smooth(
    colour=central_traj_geom_smooth_col,
    span=0.2,
    method='loess',
    fill = central_traj_geom_smooth_col,
    aes( x=pseudotime, y=expression,fill="Central Cell Lineage"),
    fullrange = TRUE)

p3 <- p3 +
  geom_smooth(data = main_HC_traj_df,
              colour=main_traj_geom_smooth_col,
              span=0.2,
              method='loess',
              fill = main_traj_geom_smooth_col,
              aes( x=pseudotime, y=expression, fill="HC Lineage"))

p3 <- p3 + scale_fill_manual(name="legend",guide = 'legend',
                                       values = c("Central Cell Lineage" = central_traj_geom_smooth_col,
                                                  "HC Lineage" = main_traj_geom_smooth_col))
# labels = c("test1", "test2"))
#values=c(central_traj_geom_smooth_col, main_traj_geom_smooth_col))

p3 <- p3 + geom_rug(data=main_HC_traj_df, sides='b', 
                              alpha=.10, aes(x=pseudotime,
                                             color = cell_group) )

p3 <- p3 + geom_rug(data=central_traj_df, sides='t', alpha=.10, 
                              aes(x=pseudotime,
                                  color = cell_group) )

p3 <- p3 + scale_color_manual(values = type_trt_cols)
p3 <-   p3 + theme_bw() +
  theme(
  ) +
  labs(
    x     = 'pseudotime'
    ,y    = 'z-scored gene expression'
  )
p3 <- p3 +theme(legend.position="none", legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

p3 <- p3 + facet_grid( .~ Gene.name.uniq , scales='fixed',cols = 1)

#add facet title to 
p3$data$title <- unique(p3$data$Gene.name.uniq)
p3 <- p3 + facet_wrap(~title)
p3 <- p3+ 
  theme(strip.text.x = element_text(size = 18))


p3


# Arrange the three components (plot, leg1, leg2)
# The two legends are positioned outside the plot: 
# one at the top and the other to the side.
plotNew <- arrangeGrob(p1, p3, 
                       heights = unit.c(p1$height, unit(.80, "npc") - p1$height), ncol = 1)

plotNew <- arrangeGrob(plotNew, p2,
                       widths = unit.c(unit(.5, "npc") - p2$width, p2$heights),
                       heights = unit.c(p2$height, unit(.5, "npc") - p2$height), nrow = 2)

grid.newpage()
grid.draw(plotNew)
lg
#===================== method 2

# 3.3 extract "legends only" from ggplot object
legend1 <- get_legend(test_lg) 
legend2 <- get_legend(test_lg2) 

# 4.1 setup legends grid
legend1_grid <- cowplot::plot_grid(legend1, align = "v", nrow = 2)

# 4.2 add second legend to grid, specifying its location
legends <- legend1_grid +
  ggplot2::annotation_custom(
    grob = legend2,
    xmin = .5, xmax = .5, ymin = .3, ymax = .85
  )

lg <- cowplot::plot_grid(p3, legends,
                         nrow = 2)
lg

# # 5. plot "plots" + "legends" (with legends in between plots)
# cowplot::plot_grid(p3, legends,
#                    nrow = 2,
#                    rel_widths = .5
# )

png(filename = "data/test_line_plot.png", width = 15,height = 18,
    units = 'in',res = 300)
cowplot::plot_grid(p3, legends,
                   nrow = 2
)
dev.off()
  

