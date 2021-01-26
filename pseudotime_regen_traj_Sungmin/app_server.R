# =========================================================================================== server 

server <- function(input, output) {


# Asks if multiple conditions are present
	whichDataset <- function() {
		if ("data.set" %in% colnames(seurat_obj@meta.data)) {
			"data.set"
		} else if ("cell.type.ident" %in% colnames(seurat_obj@meta.data)) {
			"cell.type.ident"
		} else {
			"tree.ident"}
	}

	printTreats <- reactive({
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

		n_panels <- 1:4
		
		lapply(n_panels, function(i) {
		  output[[paste0("myDatFeatPlotV", i)]] <- 
		    renderPlot({DatFeatPlotF()[[2]]})
		})
		
		lapply(n_panels, function(i) {
		  output[[paste0("plot.uiDatFeatPlotV", i)]] <- 
		    renderUI({plotOutput(paste0("myDatFeatPlotV", i),
		                         width = "425px", height = "880px")})
		})
		
#sidepanel FeaturePlot for Heatmap tab
ClustNumFeatPlotF <- function(){
  clustUMAP <- DimPlot(seurat_obj, group.by = "tree.ident")
  clustUMAP <- cleanUMAP(clustUMAP) + 
    theme(legend.position="bottom")
  return(clustUMAP)
}

# generate sidepanel FeaturePlot for Heatmap tab 
output$myDatFeatPlotV5 <- renderPlot({ClustNumFeatPlotF()[[1]]})
output$plot.uiDatFeatPlotV5 <- renderUI({
  plotOutput("myDatFeatPlotV5", width = "425px", height = "525px")
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
							size=I(as.numeric(input$nodeSize) * 1.5),
							na.rm=TRUE, branch_point_df) +
					geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
								label="branch_point_idx"),
							size=I(as.numeric(input$nodeSize)), color="white", na.rm=TRUE,
							branch_point_df)

#plot leaves
					feat[[k]] <- feat[[k]]+ 
					geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
							shape = 21, stroke=I(0.75),
							color="black",
							fill="lightgray",
							size=I(as.numeric(input$nodeSize) * 1.5),
							na.rm=TRUE,
							leaf_df) +
					geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
								label="leaf_idx"),
							size=I(as.numeric(input$nodeSize)), color="black", na.rm=TRUE, 
							leaf_df)

					feat[[k]] <- feat[[k]] + #plot root node
					geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
							shape = 21, stroke=I(0.75),
							color="black",
							fill="white",
							size=I(as.numeric(input$nodeSize) * 1.5),
							na.rm=TRUE,
							root_df) +
					geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
								label="root_idx"),
							size=I(as.numeric(input$nodeSize)), color="black", na.rm=TRUE, root_df)

			}

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
			if (l == 1) {h <- "800px"
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
  
	#download SVG button
	output$downloadSVGFeaturePlotF <- downloadHandler(
			filename = "Feature_plot.svg", content = function(file) {
			svg(file,
					width = 8, height = 8 * getLenInput(input$featureGenes))
			print(FeaturePlotF())
			dev.off()
			}
			)
	
	#download PDF button
	output$downloadPDFFeaturePlotF <- downloadHandler(
			filename = "Feature_plot.pdf", content = function(file) {
			pdf(file, 
					width = 8, height = 8 * getLenInput(input$featureGenes))
			print(FeaturePlotF())
			dev.off()
				}
			)
	
	#download pnf button
	output$downloadPNGFeaturePlotF <- downloadHandler(
	  filename = "Feature_plot.png", content = function(file) {
	    png(file, 
	        width = 8, height = 8 * getLenInput(input$featureGenes),
	        units = "in", res = 300)
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

				# central_traj_df <- make_plot_df(cds_sub = ptime_central_traj,
				# 		gene = selected)
				# 	main_HC_traj_df <- make_plot_df(cds_sub = ptime_main_traj,
				# 			gene = selected)

					central_traj_geom_smooth_col <- "#00BE67"

					main_traj_geom_smooth_col <- "#F8766D"
					
					# = add vertical line at branching point
					branching_ptime <- get_branching_point(seurat_obj = seurat_obj,
					                                       cds = cds)
					
					#if there is no cutoff on graph (use complete graph)
					if (input$selectCutOffPtimeLinePlot == "NoCutOff"){
					  central_traj_df <- make_plot_df(cds_sub = ptime_central_traj,
					                                  gene = selected)
					  main_HC_traj_df <- make_plot_df(cds_sub = ptime_main_traj,
					                                  gene = selected)
					  #if user selects for a cutoff right before the branching point
					} else {
					  central_traj_df <- make_plot_df(cds_sub = ptime_central_traj,
					                                  gene = selected)
					  central_traj_df <- central_traj_df[central_traj_df$pseudotime > 
					                                       as.numeric(branching_ptime$pseudotime - 2),]
					  main_HC_traj_df <- make_plot_df(cds_sub = ptime_main_traj,
					                                  gene = selected)
					  main_HC_traj_df <- main_HC_traj_df[main_HC_traj_df$pseudotime > 
					                                       as.numeric(branching_ptime$pseudotime - 2),]
					  
					}

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
									,y    = 'scaled gene expression'
							    )
							lg <- lg +theme(legend.position="none", legend.title=element_blank()) +
							guides(colour = guide_legend(override.aes = list(alpha = 1))) +
							scale_x_continuous( expand = c(0, 0)) #remove gap below 0
							

#add facet title to
							lg$data$title <- unique(lg$data$Gene.name.uniq)
							lg <- lg + facet_wrap(~title)
							lg <- lg+
							theme(strip.text.x = element_text(size = 18))
							
# add branching vertical line
							lg <- lg +
							  new_scale_color() + # add new scale color for geom_rug
							  new_scale_fill() +# add new scale color for geom_rug
							  geom_vline(data = branching_ptime, 
							             aes(xintercept = pseudotime, color = "cell_group"),
							             linetype=2) +
							  scale_color_manual(name = " ",  values = "red",
							                     labels = c("branching point"))
							
							if(selected >1){
							  lg <- lg + facet_wrap(~Gene.name.uniq, ncol=1)
							}

					}else{ #withLegend
					  lg <- ggplot(central_traj_df) +
					    geom_smooth(
					      #colour=central_traj_geom_smooth_col,
					      span=0.2,
					      method='loess',
					      #fill = central_traj_geom_smooth_col,
					      aes( x=pseudotime, y=expression,color="Central Cell Lineage", fill = "Central Cell Lineage"),
					      fullrange = TRUE)
					  
					  lg <- lg +
					    geom_smooth(data = main_HC_traj_df,
					                #colour=main_traj_geom_smooth_col,
					                span=0.2,
					                method='loess',
					                #fill = main_traj_geom_smooth_col,
					                aes( x=pseudotime, y=expression, color="HC Lineage", fill = "HC Lineage"))
					  
					  lg <- lg + scale_color_manual(name="Branching Trajectories",
					                                guide = 'legend',
					                                values = c("Central Cell Lineage" = central_traj_geom_smooth_col,
					                                           "HC Lineage" = main_traj_geom_smooth_col)) + guides(guide_legend(override.aes = list(linetype = c("black","black"))))
					  lg <- lg + scale_fill_manual(name="Branching Trajectories",guide = 'legend',
					                               values = c("Central Cell Lineage" = central_traj_geom_smooth_col,
					                                          "HC Lineage" = main_traj_geom_smooth_col))
					  lg <- lg  +
					    new_scale_color() + # add new scale color for branching point
					    new_scale_fill() # add new scale color for branching point
					  
					  # add branching vertical line
					  lg <- lg +
					    new_scale_color() + # add new scale color for geom_rug
					    new_scale_fill() +# add new scale color for geom_rug
					    geom_vline(data = branching_ptime, 
					               aes(xintercept = pseudotime, color = "cell_group"),
					               linetype=2) +
					    scale_color_manual(name = " ",  values = "red",
					                       labels = c("branching point"))
					  
					  lg <- lg  +
					    new_scale_color() + # add new scale color for geom_rug
					    new_scale_fill() # add new scale color for geom_rug
					  
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
					      ,y    = 'scaled gene expression'
					    )
					  lg <- lg +theme(legend.position="bottom",
					                  legend.box = "vertical",
					                  legend.title=element_blank()) +
					    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
					    scale_x_continuous( expand = c(0, 0)) #remove gap below 0
					  
					  
					  #add facet title to
					  lg$data$title <- unique(lg$data$Gene.name.uniq)
					  lg <- lg + facet_wrap(~title)
					  lg <- lg+
					    theme(strip.text.x = element_text(size = 18))
					    
					    #handle multiple genes
							if(selected >1){
							  lg <- lg + facet_wrap(~Gene.name.uniq, ncol=1)
							}

					}

				return(lg)
		})

	mismatchPtimeLinePlot <- function(seurat_obj) {

		selected <- unique(unlist(strsplit(input$ptimeLinePlotGenes, " ")))

			mismatch <- ifelse(!selected %in% c(com_name, ens_id),
					selected[!selected %in% c(com_name, ens_id,
						rownames(seurat_obj[["RNA"]]@data))],"")
			return(mismatch)
	}

	output$notInPtimeLinePlot <- renderText({input$runPtimeLinePlot
			isolate({mismatchPtimeLinePlot()})
			})


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

				if (l == 1) {h <- "400px"
				} else {
					h <- as.character(ceiling(l) * 400)
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

	#Download SVG
	output$downloadSVGPtimeLinePlotF <- downloadHandler(
	  filename = "Ptime_line_dynamic_plot.svg", content = function(file) {
	    if (input$selectGrpPtimeLinePlot == "NoLegend"){
	      svg(file, 
	          width = 12, height = 5 * getLenInput(input$ptimeLinePlotGenes))
	      print(PtimeLinePlotF())
	      dev.off()
	    }else{
	      svg(file, 
	          width = 12, height = 7 * getLenInput(input$ptimeLinePlotGenes))
	      print(PtimeLinePlotF())
	      dev.off()
	    }
	  }
	)
	
	#Download PDF
	output$downloadPDFPtimeLinePlotF <- downloadHandler(
	  filename = "Ptime_line_dynamic_plot.pdf", content = function(file) {
	    if (input$selectGrpPtimeLinePlot == "NoLegend"){
	      pdf(file, 
	          width = 12, height = 5 * getLenInput(input$ptimeLinePlotGenes))
	      print(PtimeLinePlotF())
	      dev.off()
	    }else{
	      pdf(file, 
	          width = 12, height = 7 * getLenInput(input$ptimeLinePlotGenes))
	      print(PtimeLinePlotF())
	      dev.off()
	    }
	  }
	)
	
	#Download PNG
	output$downloadPNGPtimeLinePlotF <- downloadHandler(
	  filename = "Ptime_line_dynamic_plot.png", content = function(file) {
	    if (input$selectGrpPtimeLinePlot == "NoLegend"){
	      png(file, 
	          width = 12, height = 5 * getLenInput(input$ptimeLinePlotGenes),
	          units = "in", res = 300)
	      print(PtimeLinePlotF())
	      dev.off()
	    }else{
	      png(file, 
	          width = 12, height = 7 * getLenInput(input$ptimeLinePlotGenes),
	          units = "in", res = 300)
	      print(PtimeLinePlotF())
	      dev.off()
	    }
	  }
	)

	# # ======== Mulitple Gene Ptime Line Dynamics ======== #
	MultiGPtimeLinePlotF <- reactive({
	  selected <- unique(unlist(strsplit(input$multptimeLinePlotGenes, " ")))
	  
	  ifelse(selected %in% com_name,
	         selected <- selected[selected %in% com_name],
	         
	         ifelse(selected %in% ens_id,
	                selected <- gene_df[ens_id %in% selected, 3],"")
	  )
	  
	  #remove not found in norm matrix obj
	  selected <- selected[selected %in% rownames(seurat_obj[["RNA"]]@data)]
	  
	  # = add vertical line at branching point
	  branching_ptime <- get_branching_point(seurat_obj = seurat_obj, cds = cds)
	  
	  #if there is no cutoff on graph (use complete graph)
	  if (input$selectCutOffMultiGPtimeLinePlot == "NoCutOff"){
	  central_traj_df <- make_plot_df(cds_sub = ptime_central_traj,
	                                  gene = selected)
	  main_HC_traj_df <- make_plot_df(cds_sub = ptime_main_traj,
	                                  gene = selected)
	  #if user selects for a cutoff right before the branching point
	  } else {
	    central_traj_df <- make_plot_df(cds_sub = ptime_central_traj,
	                                    gene = selected)
	    central_traj_df <- central_traj_df[central_traj_df$pseudotime > 
	                                 as.numeric(branching_ptime$pseudotime - 2),]
	    main_HC_traj_df <- make_plot_df(cds_sub = ptime_main_traj,
	                                    gene = selected)
	    main_HC_traj_df <- main_HC_traj_df[main_HC_traj_df$pseudotime > 
	                                 as.numeric(branching_ptime$pseudotime - 2),]
	    
	  }
	  #order genes in order of input sequence 
	  central_traj_df$Gene.name.uniq <- factor(central_traj_df$Gene.name.uniq, 
	                                           levels = selected)
	  main_HC_traj_df$Gene.name.uniq <- factor(main_HC_traj_df$Gene.name.uniq, 
	                                           levels = selected)
	  
	  #specify trajectory trail
	  central_traj_df$trajectory <- "Central Cell Lineage"
	  main_HC_traj_df$trajectory <- "HC Lineage"
	  
	  #cancatenate two traj df 
	  plot_dt <- rbind(central_traj_df,main_HC_traj_df)
	  
	  central_traj_geom_smooth_col <- "#00BE67"
	  
	  main_traj_geom_smooth_col <- "#F8766D"
	  
	  #add smoothing gene line
	  p <- ggplot(data = plot_dt,
	              mapping = aes(x = pseudotime, y = expression, color = Gene.name.uniq,
	                            fill = Gene.name.uniq)) +
	    geom_smooth(method ="loess", span = 0.2, fullrange = TRUE) 
	  
	  #add vertical line at branching point
	  p <- p +
	    new_scale_color() + # add new scale color for geom_rug
	    new_scale_fill() +# add new scale color for geom_rug
	    geom_vline(data = branching_ptime, aes(xintercept = pseudotime, color = "cell_group"),
	               linetype=2) +
	    scale_color_manual(name = " ",  values = "red",
	                       labels = c("branching point"))
	  
	  #add geom_rug for dashed cell types
	  p <- p+
	    new_scale_color() + # add new scale color for geom_rug
	    new_scale_fill() +# add new scale color for geom_rug
	    geom_rug(data=plot_dt[plot_dt$trajectory == "Central Cell Lineage",], sides='b', 
	             alpha=.10, aes(x=pseudotime, color = cell_group) ) +
	    scale_color_manual(values = type_trt_cols, guide = "none") +
	    new_scale_color() + # add new scale color for geom_rug
	    new_scale_fill() +# add new scale color for geom_rug
	    geom_rug(data=plot_dt[plot_dt$trajectory == "HC Lineage",], sides='b', 
	             alpha=.10, aes(x=pseudotime, color = cell_group) ) +
	    scale_color_manual(values = type_trt_cols) 
	  
	  #legend themes
	  p <- p + theme_bw() +
	    guides(colour = guide_legend(override.aes = list(alpha = 1)))+
	    theme(
	    ) +
	    labs(
	      x     = 'pseudotime'
	      ,y    = 'scaled gene expression'
	    ) + 
	    theme(legend.position="bottom", 
	          legend.box = "vertical",
	          legend.title=element_blank()) +
	    facet_wrap(~ trajectory) + #split plot by trajectory path
	    theme(strip.text.x = element_text(size = 18, face = "bold")) + #+#specify facet title size 
	    scale_x_continuous(breaks= scales::pretty_breaks()) + #create automatic x axis labels
	    scale_x_continuous( expand = c(0, 0)) #remove gap below 0
	  
	    #dynamically adjust based on high and low geom smooth values y axis limits without dropping data
	   # p <- p + coord_cartesian(ylim = c(min(ggplot_build(p)$data[[1]]$ymin,na.rm = TRUE), #dynamically plot min  and max y lim
	   #                                   max(ggplot_build(p)$data[[1]]$ymax,na.rm = TRUE)))
	   # warning message: remove n rows means there was NA's in ggplot_build(p)$data[[1]]$ymin
	 
	  # # manually apply facet panel colors 
	  g <- ggplot_gtable(ggplot_build(p))
	  strips <- which(grepl('strip-', g$layout$name))
	  
	  pal <- c(central_traj_geom_smooth_col, main_traj_geom_smooth_col)
	  
	  for (i in seq_along(strips)) {
	    k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	    l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	    g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i] #change facet  background label colors
	    g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- "white" #change facet letter colors
	  }
	  
	  return(g)
	  #return(g)
	})
	
	mismatchMultiPtimeLinePlot <- function(seurat_obj) {
	  #seurat_obj <- SelectDataset()
	  
	  selected <- unique(unlist(strsplit(input$multptimeLinePlotGenes, " ")))
	  
	  mismatch <- ifelse(!selected %in% c(com_name, ens_id),
	                     selected[!selected %in% c(com_name, ens_id,
	                                               rownames(seurat_obj[["RNA"]]@data))],"")
	  return(mismatch)
	}
	
	output$notInMultiPtimeLinePlot <- renderText({input$runMultiPtimeLinePlot
	  isolate({mismatchMultiPtimeLinePlot()})
	})
	
	
	output$myMultiGPtimeLinePlotF <- renderPlot({input$runMultiPtimeLinePlot
	  isolate({withProgress({p <- MultiGPtimeLinePlotF(); print(plot(p))},
	                        message = "Rendering plot..", min = 0, max = 10, value = 10)})
	})
	
	#do not need getHeightMultiPtimeLinePlot() for this tab
	# getHeightMultiPtimeLinePlot <- function() {
	#   l <- getLenInput(input$multptimeLinePlotGenes)
	#     if (l == 1) {h <- "300px"
	#     } else {
	#       h <- as.character(ceiling(l) * 300)
	#       h <- paste0(h, "px")
	#     }
	#     if (l == 1) {h <- "600px"
	#     } else {
	#       h <- as.character(ceiling(l) * 600)
	#       h <- paste0(h, "px")
	#     }
	#   }
	#   return(h)
	# }
	
	output$plot.uiMultiGPtimeLinePlotF <- renderUI({input$runMultiPtimeLinePlot
	  #isolate({h <- getHeightMultiPtimeLinePlot()
	  plotOutput("myMultiGPtimeLinePlotF", width = "1000px", height = "700px")
	  })
	#})
	
	#Download SVG
	output$downloadSVGMultiGPtimeLinePlotF <- downloadHandler(
	  filename = "Multi_Ptime_line_dynamic_plot.svg", content = function(file) {
	    # if (input$selectGrpMultiPtimeLinePlot == "NoLegend"){
	    #   svg(file, 
	    #       width = 12, height = 5)
	    #   print(MultiGPtimeLinePlotF())
	    #   dev.off()
	    # }else{
	      svg(file, 
	          width = 14, height = 7.5)
	      print(plot(MultiGPtimeLinePlotF()))
	      dev.off()
	   # }
	  }
	)
	
	
	#Download PDF
	output$downloadPDFMultiGPtimeLinePlotF <- downloadHandler(
	  filename = "Multi_Ptime_line_dynamic_plot.pdf", content = function(file) {
	      pdf(file, 
	          width = 14, height = 7.5)
	    print(plot(MultiGPtimeLinePlotF()))
	    dev.off()
	  }
	)
	
	#Download PNG
	output$downloadPNGMultiGPtimeLinePlotF <- downloadHandler(
	  filename = "Multi_Ptime_line_dynamic_plot.png", content = function(file) {
	      png(file, 
	          width = 14, height = 7.5,
	          units = "in", res = 300)
	      print(plot(MultiGPtimeLinePlotF()))
	      dev.off()
	  }
	)
	
	# # ======== Cluster heatmap ======== #
	HeatmapF <- reactive({
	  selected <- unique(unlist(strsplit(input$hmapGenes, " ")))
	  
	  ifelse(selected %in% com_name,
	         selected <- selected[selected %in% com_name],
	         
	         ifelse(selected %in% ens_id,
	                selected <- gene_df[ens_id %in% selected, 3],"")
	  )
	  
	  #remove not found in norm matrix obj
	  selected <- selected[selected %in% rownames(seurat_obj[["RNA"]]@data)]
	  
	  # === build ggplot2 heatmap object
	  # change current idents to heirarachical clusters num
	  seurat_obj$tree.ident <- factor(seurat_obj$tree.ident)
	  levels(seurat_obj$tree.ident) <- paste0("cluster",levels(seurat_obj$tree.ident))
	  Idents(seurat_obj) <- "tree.ident"
	  
	  dotplot <- DotPlot(seurat_obj, features = selected,
	                     group.by = "tree.ident")
	  
	  # specify trajectory group based on cluster number
	  clusters <- as.character(dotplot$data$id)
	  #specify treatment and sequencing technique (data.set)
	  dotplot$data <- dotplot$data %>% mutate(groupIdent=case_when(
	    str_detect(dotplot$data$id, "cluster3|cluster2|cluster7|cluster6") ~ "initial injury response", 
	    str_detect(dotplot$data$id, "cluster8|cluster5|^cluster1$") ~ "HC Lineage",
	    str_detect(dotplot$data$id, "cluster9|cluster4|^cluster10$") ~ "Central Cell Lineage",
	    TRUE ~ as.character(dotplot$data$id)))
	  
	  #reorder levels/ build plotting df
	  dotplot$data$groupIdent <- factor(dotplot$data$groupIdent,
	                                    levels = c("initial injury response",
	                                               "HC Lineage",
	                                               "Central Cell Lineage"))
	  
	  dotplot$data$id <- factor(dotplot$data$id, 
	                            levels = paste0("cluster",c(3,2,7,6,8,5,1,9,10,4)))
	  
	  plot_df <- dotplot$data
	  
	  # build ggplot obj
	  p <- ggplot(plot_df, aes(id, features.plot,fill= avg.exp.scaled, width = 1, height = 1)) +
	    geom_tile(color = "gray", size = .75) +
	    scale_x_discrete(expand=c(0,0)) + #expand tile 
	    scale_y_discrete(expand=c(0,0))+ #expand tile
	    scale_fill_distiller(
	      palette = "RdYlBu") +
	    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5,size = 13),
	          axis.title.y.right = element_blank(),
	          panel.spacing = unit(.25, "lines"), #change facet spacing
	          strip.text.x  = element_text(size = 14,face = "bold"), #change facet label 
	          panel.background = element_blank(),
	          axis.title.y  = element_blank()) +
	    facet_grid( ~ groupIdent, scales='free_x') + #split heatmap by meta trajectory
	    ylim(rev(levels(dotplot$data$features.plot))) #put first element in gene list on top of heatmap
	  #scale_y_discrete(position = "right") # <---- uncomment if you want to put gene labels on the right side
	  
	  #= manually apply facet panel colors 
	  g <- ggplot_gtable(ggplot_build(p))
	  strips <- which(grepl('strip-', g$layout$name))
	  
	  pal <- c("#A3A500", "#F8766D", "#39B600")
	  
	  for (i in seq_along(strips)) {
	    k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	    l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	    g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i] #change facet  background label colors
	    g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- "white" #change facet letter colors
	  }
	  return(g)
	})
	
	mismatchHmap <- function(seurat_obj) {
	  selected <- unique(unlist(strsplit(input$hmapGenes, " ")))
	  
	  mismatch <- ifelse(!selected %in% c(com_name, ens_id),
	                     selected[!selected %in% c(com_name, ens_id,
	                                               rownames(seurat_obj[["RNA"]]@data))],"")
	  return(mismatch)
	}
	
	output$notInHmap <- renderText({input$runHmap
	  isolate({mismatchHmap()})
	})
	
	
	output$myHeatmapF <- renderPlot({input$runHmap
	  isolate({withProgress({p <- HeatmapF(); print(plot(p))},
	                        message = "Rendering plot..", 
	                        min = 0, 
	                        max = 10, 
	                        value = 10)})
	})
	
	# getHeightHmap <- function() {
	#   l <- getLenInput(input$hmapGenes)
	#     if (l == 1) {h <- "300px"
	#     } else {
	#       h <- as.character(ceiling(l) * 15)
	#       h <- paste0(h, "px")
	#     }
	#   }
	#   return(h)
	# }
	
	output$plot.uiHeatmapF <- renderUI({input$runHmap
	  #isolate({h <- getHeightHmap()
	  plotOutput("myHeatmapF",
	             width = paste0(input$manAdjustHmapW, "in"),
	             height = paste0(input$manAdjustHmapH, "in"))})

	
	#Download SVG
	output$downloadSVGHeatmapF <- downloadHandler(
	  filename = "Heatmap.svg", content = function(file) {
	    svg(file, 
	        width = as.numeric(input$manAdjustHmapW), 
	        height = as.numeric(input$manAdjustHmapH))
	    print(plot(HeatmapF()))
	    dev.off()
	  }
	)
	
	
	#Download PDF
	output$downloadPDFHeatmapF <- downloadHandler(
	  filename = "Heatmap.pdf", content = function(file) {
	    pdf(file, 
	        width = as.numeric(input$manAdjustHmapW), 
	        height = as.numeric(input$manAdjustHmapH))
	    print(plot(HeatmapF()))
	    dev.off()
	  }
	)
	
	#Download PNG
	output$downloadPNGHeatmapF <- downloadHandler(
	  filename = "Heatmap.png", content = function(file) {
	    png(file, 
	        width = as.numeric(input$manAdjustHmapW), 
	        height = as.numeric(input$manAdjustHmapH),
	        units = "in", res = 300)
	    print(plot(HeatmapF()))
	    dev.off()
	  }
	)

} # Server close
