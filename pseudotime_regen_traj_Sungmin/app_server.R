# =========================================================================================== server 

server <- function(input, output) {

# ======== Dataset selection ======== #
# SelectDataset <- reactive({
# 	seurat_obj <- file_list[[input$Analysis]]
# 	print(names(file_list[input$Analysis]))
# 
# 	cluster_clrs <<- gg_color_hue(
# 				      length(levels(seurat_obj@active.ident)))
# 	return(seurat_obj)
# })

# Asks if multiple conditions are present
	whichDataset <- function() {
#seurat_obj <- SelectDataset()
		if ("data.set" %in% colnames(seurat_obj@meta.data)) {
			"data.set"
		} else if ("cell.type.ident" %in% colnames(seurat_obj@meta.data)) {
			"cell.type.ident"
		} else {
			"tree.ident"}
	}

	printTreats <- reactive({
#seurat_obj <- SelectDataset()
			print(seurat_obj)
			if (whichDataset() == "data.set") {
			sort(unique(seurat_obj@meta.data$data.set))
			} else {
			NULL # single data set
			}
			})

	printIdents <- reactive({
#seurat_obj <- SelectDataset()
			print(seurat_obj)
			if(input$Analysis %in% multiple_idents_seurObj){
			sort(unique(seurat_obj@meta.data$cell.type.ident))
			} else {
			sort(unique(seurat_obj@meta.data$seurat_clusters))
			}
			})

# returns the correct ID class for cell subset
	IDtype <- function() {
#seurat_obj <- SelectDataset()
		if(input$Analysis %in% multiple_idents_seurObj){
			seurat_obj@meta.data$cell.type.ident
		} else {
			seurat_obj@meta.data$seurat_clusters
		}
	}


# ======== Gene Database ======== #
	GeneDB <- function() {
#seurat_obj <- SelectDataset()
		selected <- unique(unlist(strsplit(input$dbGenes, " ")))

			present <- gene_df$Gene.name.uniq %in% rownames(seurat_obj)
			gene_df <- cbind(in_dataset = present, gene_df)

			ifelse(length(selected) < 2,
					ind <- grep(selected, gene_df$Gene.name.uniq),

					ifelse(selected %in% gene_df$Gene.name.uniq,
						ind <- multiGrep2(selected, gene_df$Gene.name.uniq),

						ifelse(selected %in% gene_df$Gene.stable.ID,
							ind <- multiGrep2(selected, gene_df$Gene.stable.ID),
							"Gene not in dataset")
					      )
			      )
			gene_df[ind,]
	}
	output$GeneDB <- renderTable({GeneDB()})


# ======== Cluster/Data UMAP ======== #
		DatFeatPlotF <- function() {

# ============== traj umap by dataset
			traject_plot <- plot_cells(cds, color_cells_by = "cell.type.and.trt",
					label_leaves = TRUE, 
					show_trajectory_graph = TRUE, 
					cell_size = 0.70,
					trajectory_graph_segment_size = 0.55, 
					label_cell_groups = FALSE,
					group_label_size = 8) + 
				theme(legend.position="bottom", legend.title=element_blank())

# #preserve dataset colors 
# cell_type_trt <- levels(cds$cell.type.and.trt)
# type_trt_cols <- gg_color_hue(length(cell_type_trt))

				traject_plot <- cleanUMAP(traject_plot)
				traject_plot <- traject_plot + scale_color_manual(values = type_trt_cols)

# ================================= traj umap by pseudotime
				ptime_plot <- plot_cells(cds, color_cells_by = "pseudotime",
						label_leaves = TRUE, 
						show_trajectory_graph = TRUE, 
						cell_size = 0.70,
						trajectory_graph_segment_size = 0.55, 
						label_cell_groups = FALSE)

				ptime_plot <- ptime_plot +
				guides(fill = guide_legend(ncol = 1))
				ptime_plot <- cleanUMAP(ptime_plot)


				datfeat_list <- list(traject_plot, ptime_plot)
				plot_h <- plot_grid(plotlist = datfeat_list, ncol = 2)

				plot_v <- plot_grid(plotlist = datfeat_list, ncol = 1)
				plot_v <- plot_grid(plot_v) + theme(
						plot.background = element_rect(size = 2, color = "#DCDCDC"))

				datfeat_list <- list(plot_h, plot_v)
				datfeat_list
		}

	output$myDatFeatPlotH1 <- renderPlot({DatFeatPlotF()[[1]]})
		output$plot.uiDatFeatPlotH1 <- renderUI({
				plotOutput("myDatFeatPlotH1", width = "1050px", height = "525px")
				})

	n_panels <- 1:3

		lapply(n_panels, function(i) {
				output[[paste0("myDatFeatPlotV", i)]] <- 
				renderPlot({DatFeatPlotF()[[2]]})
				})

	lapply(n_panels, function(i) {
			output[[paste0("plot.uiDatFeatPlotV", i)]] <- 
			renderUI({plotOutput(paste0("myDatFeatPlotV", i),
						width = "425px", height = "880px")})
			})


# ======== Feature Plot ======== #
	FeaturePlotF <- reactive({
#seurat_obj <- SelectDataset()
			selected <- unique(unlist(strsplit(input$featureGenes, " ")))

			ifelse(selected %in% com_name,
					selected <- selected[selected %in% com_name],

					ifelse(selected %in% ens_id,
						selected <- gene_df[ens_id %in% selected, 3],"")
			      )

#remove not found in scaled matrix obj
			selected <- selected[selected %in% rownames(seurat_obj[["RNA"]]@data)]

# ============ import principle graph df, root node, branch df, 
# ============ leaf nodes df 
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

			feat <- FeaturePlot(seurat_obj, selected,
					reduction = "umap", cols = c(input$CellBackCol, input$CellForeCol),
					combine = FALSE, pt.size = input$ptSizeFeature)

			for(k in 1:length(feat)) {
				feat[[k]] <- feat[[k]] + labs(x = "UMAP 1", y = "UMAP 2") + 
					theme(axis.text.x = element_blank(), legend.position="none",
							axis.ticks.x = element_blank(), 
							axis.line.x = element_blank(),
							axis.text.y = element_blank(), 
							axis.ticks.y = element_blank(),
							axis.line.y = element_blank(), 
							axis.title = element_text(size = 18),
							panel.border = element_rect(colour = "#FFFFFF",
								fill = NA, size = 1))

#import, plot traj line
					feat[[k]] <- feat[[k]] + geom_segment(
							aes_string(x="source_prin_graph_dim_1",
								y="source_prin_graph_dim_2",
								xend="target_prin_graph_dim_1",
								yend="target_prin_graph_dim_2"),
							linetype="solid",
							na.rm=TRUE,
							data=edge_df) 

#plot branching
					feat[[k]] <- feat[[k]] + 
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

#plot leaves
					feat[[k]] <- feat[[k]]+ 
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

					feat[[k]] <- feat[[k]] + #plot root node
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

			}
# return(plot_grid(plotlist = feat, ncol = 1))

			pg <- plot_grid(plotlist = feat, ncol = 1) +
				labs(title = paste("Selected analysis:",
							as.character("HC Lineage Trajectory")), 
						subtitle = "", caption = "") +
				theme(plot.title = element_text(face = "bold", size = 15, hjust = 0))

				return(pg)
	})

	mismatchFeat <- function(seurat_obj) {
#seurat_obj <- SelectDataset()

		selected <- unique(unlist(strsplit(input$featureGenes, " ")))

			mismatch <- ifelse(!selected %in% c(com_name, ens_id),
					selected[!selected %in% c(com_name, ens_id,
						rownames(seurat_obj[["RNA"]]@data))],"")
			return(mismatch)
	}

	output$notInFeat <- renderText({input$runFeatPlot
			isolate({mismatchFeat()})
			})

# output$SelectedDataFeat <- renderText({input$runFeatPlot
# 	isolate({input$Analysis})
# 			   })

	output$myFeaturePlotF <- renderPlot({input$runFeatPlot
			isolate({withProgress({p <- FeaturePlotF(); print(p)},
						message = "Rendering plot..", min = 0, max = 10, value = 10)})
			})

	getHeightFeat <- function() {
		l <- getLenInput(input$featureGenes)
			if (l == 1) {h <- "1000px"
			} else {
				h <- as.character(ceiling(l) * 800)
					h <- paste0(h, "px")
			}
		return(h)
	}

	output$plot.uiFeaturePlotF <- renderUI({input$runFeatPlot
			isolate({h <- getHeightFeat()
					plotOutput("myFeaturePlotF", width = "840px", height = h)
					})
			})

# output$downloadFeaturePlotF <- downloadHandler(
# 					       filename = "Feature_plot.png", content = function(file) {
# 						       png(file, units = "in", res = as.numeric(input$featDPI),
# 							   width = 12, height = 12 * getLenInput(input$featureGenes))
# 						       print(FeaturePlotF())
# 						       dev.off()
# 					       }
# 					       )

	output$downloadFeaturePlotF <- downloadHandler(
			filename = "Feature_plot.svg", content = function(file) {
			svg(file,
					width = 12, height = 12 * getLenInput(input$featureGenes))
			print(FeaturePlotF())
			dev.off()
			}
			)

# # ======== Indv Gene Ptime Line Dynamics ======== #
		PtimeLinePlotF <- reactive({
#seurat_obj <- SelectDataset()
				selected <- unique(unlist(strsplit(input$ptimeLinePlotGenes, " ")))

				ifelse(selected %in% com_name,
						selected <- selected[selected %in% com_name],

						ifelse(selected %in% ens_id,
							selected <- gene_df[ens_id %in% selected, 3],"")
				      )

#remove not found in norm matrix obj
				selected <- selected[selected %in% rownames(seurat_obj[["RNA"]]@data)]

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

				central_traj_df <- make_plot_df(cds_sub = ptime_central_traj,
						gene = selected)
					main_HC_traj_df <- make_plot_df(cds_sub = ptime_main_traj,
							gene = selected)

					central_traj_geom_smooth_col <- "#00BE67"

					main_traj_geom_smooth_col <- "#F8766D"

					if (input$selectGrpPtimeLinePlot == "NoLegend"){
						lg <- ggplot(central_traj_df) +
							geom_smooth(
									colour=central_traj_geom_smooth_col,
									span=0.2,
									method='loess',
									fill = central_traj_geom_smooth_col,
									aes( x=pseudotime, y=expression,fill="Central Cell Lineage"),
									fullrange = TRUE)

							lg <- lg +
							geom_smooth(data = main_HC_traj_df,
									colour=main_traj_geom_smooth_col,
									span=0.2,
									method='loess',
									fill = main_traj_geom_smooth_col,
									aes( x=pseudotime, y=expression, fill="HC Lineage"))

							lg <- lg + scale_fill_manual(name="legend",guide = 'legend',
									values = c("Central Cell Lineage" = central_traj_geom_smooth_col,
										"HC Lineage" = main_traj_geom_smooth_col))

							lg <- lg + geom_rug(data=main_HC_traj_df, sides='b',
									alpha=.10, aes(x=pseudotime,
										color = cell_group) )

							lg <- lg + geom_rug(data=central_traj_df, sides='t', alpha=.10,
									aes(x=pseudotime,
										color = cell_group) )

							lg <- lg + scale_color_manual(values = type_trt_cols)
							lg <-   lg + theme_bw() +
							theme(
							     ) +
							labs(
									x     = 'pseudotime'
									,y    = 'z-scored gene expression'
							    )
							lg <- lg +theme(legend.position="none", legend.title=element_blank()) +
							guides(colour = guide_legend(override.aes = list(alpha = 1)))

#add facet title to
							lg$data$title <- unique(lg$data$Gene.name.uniq)
							lg <- lg + facet_wrap(~title)
							lg <- lg+
							theme(strip.text.x = element_text(size = 18))

					}else{
						plt_legend1 <- ggplot(central_traj_df) +
							geom_smooth(
#colour=central_traj_geom_smooth_col,
									span=0.2,
									method='loess',
#fill = central_traj_geom_smooth_col,
									aes( x=pseudotime, y=expression,color="Central Cell Lineage", fill = "Central Cell Lineage"),
									fullrange = TRUE)

							plt_legend1 <- plt_legend1 +
							geom_smooth(data = main_HC_traj_df,
#colour=main_traj_geom_smooth_col,
									span=0.2,
									method='loess',
#fill = main_traj_geom_smooth_col,
									aes( x=pseudotime, y=expression, color="HC Lineage", fill = "HC Lineage"))

							plt_legend1 <- plt_legend1 + scale_color_manual(name="Branching Trajectories",guide = 'legend',
									values = c("Central Cell Lineage" = central_traj_geom_smooth_col,
										"HC Lineage" = main_traj_geom_smooth_col)) + guides(guide_legend(override.aes = list(linetype = c("black","black"))))
							plt_legend1 <- plt_legend1 + scale_fill_manual(name="Branching Trajectories",guide = 'legend',
									values = c("Central Cell Lineage" = central_traj_geom_smooth_col,
										"HC Lineage" = main_traj_geom_smooth_col))
							plt_legend1 <- plt_legend1 + theme(legend.direction = "horizontal",legend.position = "bottom")

							plt_legend1 <- plt_legend1 +
							guides(colour = guide_legend(title.position = "top",title.hjust = 0.5))



							plt_legend2 <- ggplot(central_traj_df) +
							geom_smooth(
									colour=central_traj_geom_smooth_col,
									span=0.2,
									method='loess',
									fill = central_traj_geom_smooth_col,
									aes( x=pseudotime, y=expression,color="Central Cell Lineage", fill = "Central Cell Lineage"),
									fullrange = TRUE)

							plt_legend2 <- plt_legend2 +
							geom_smooth(data = main_HC_traj_df,
									colour=main_traj_geom_smooth_col,
									span=0.2,
									method='loess',
									fill = main_traj_geom_smooth_col,
									aes( x=pseudotime, y=expression, color="HC Lineage", fill = "HC Lineage"))


							plt_legend2 <- plt_legend2 + geom_rug(data=main_HC_traj_df, sides='b',
									alpha=.10, aes(x=pseudotime,
										color = cell_group) )

							plt_legend2 <- plt_legend2 + geom_rug(data=central_traj_df, sides='t', alpha=.10,
									aes(x=pseudotime,
										color = cell_group) )

							plt_legend2 <- plt_legend2 + scale_color_manual(values = type_trt_cols)
							plt_legend2 <-   plt_legend2 + theme_bw() +
							theme(
							     ) +
							labs(
									x     = 'pseudotime'
									,y    = 'z-scored gene expression'
							    )
							plt_legend2 <- plt_legend2 +theme(legend.position="bottom", legend.title=element_blank()) +
							guides(colour = guide_legend(override.aes = list(alpha = 1)))

#add facet title to
							plt_legend2$data$title <- unique(plt_legend2$data$Gene.name.uniq)
							plt_legend2 <- plt_legend2 + facet_wrap(~title)
							plt_legend2 <- plt_legend2+
							theme(strip.text.x = element_text(size = 18))

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

#add facet title to
							p3$data$title <- unique(p3$data$Gene.name.uniq)
							p3 <- p3 + facet_wrap(~title)
							p3 <- p3+
							theme(strip.text.x = element_text(size = 18))

# 3.3 extract "legends only" from ggplot object
							legend1 <- get_legend(plt_legend1)
							legend2 <- get_legend(plt_legend2)

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
					}

				return(lg)
		})

	mismatchPtimeLinePlot <- function(seurat_obj) {
#seurat_obj <- SelectDataset()

		selected <- unique(unlist(strsplit(input$ptimeLinePlotGenes, " ")))

			mismatch <- ifelse(!selected %in% c(com_name, ens_id),
					selected[!selected %in% c(com_name, ens_id,
						rownames(seurat_obj[["RNA"]]@data))],"")
			return(mismatch)
	}

	output$notInPtimeLinePlot <- renderText({input$runPtimeLinePlot
			isolate({mismatchPtimeLinePlot()})
			})

# output$SelectedDataFeat <- renderText({input$runPtimeLinePlot
# 	isolate({input$Analysis})
# 			   })

	output$myPtimeLinePlotF <- renderPlot({input$runPtimeLinePlot
			isolate({withProgress({p <- PtimeLinePlotF(); print(p)},
						message = "Rendering plot..", min = 0, max = 10, value = 10)})
			})

	getHeightPtimeLinePlot <- function() {
		l <- getLenInput(input$ptimeLinePlotGenes)
			if (input$selectGrpPtimeLinePlot == "NoLegend"){

				if (l == 1) {h <- "300px"
				} else {
					h <- as.character(ceiling(l) * 300)
						h <- paste0(h, "px")
				}
			}else{ #WithLegend

				if (l == 1) {h <- "600px"
				} else {
					h <- as.character(ceiling(l) * 600)
						h <- paste0(h, "px")
				}
			}
		return(h)
	}

	output$plot.uiPtimeLinePlotF <- renderUI({input$runPtimeLinePlot
			isolate({h <- getHeightPtimeLinePlot()
					plotOutput("myPtimeLinePlotF", width = "740px", height = h)
					})
			})

# output$downloadPtimeLinePlotF <- downloadHandler(
#   filename = "Feature_plot.png", content = function(file) {
#     png(file, units = "in", res = as.numeric(input$featDPI),
#         width = 12, height = 12 * getLenInput(input$ptimeLinePlotGenes))
#     print(PtimeLinePlotF())
#     dev.off()
#   }
# )

	output$downloadPtimeLinePlotF <- downloadHandler(
			filename = "Ptime_line_dynamic_plot.svg", content = function(file) {
			if (input$selectGrpPtimeLinePlot == "NoLegend"){
			svg(file, 
					width = 12, height = 5 * getLenInput(input$ptimeLinePlotGenes))
			print(PtimeLinePlotF())
			dev.off()
			}else{
			svg(file, 
					width = 12, height = 9 * getLenInput(input$ptimeLinePlotGenes))
			print(PtimeLinePlotF())
			dev.off()
			}
			}
			)



} # Server close
