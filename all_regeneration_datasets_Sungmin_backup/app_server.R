server <- function(input, output) {

  # ======== Dataset selection ======== #
  SelectDataset <- reactive({
    seurat_obj <- file_list[[input$Analysis]]
    print(names(file_list[input$Analysis]))

    cluster_clrs <<- gg_color_hue(
      length(levels(seurat_obj@active.ident)))
    return(seurat_obj)
  })

  # Asks if multiple conditions are present
  whichDataset <- function() {
    seurat_obj <- SelectDataset()
    if ("data.set" %in% colnames(seurat_obj@meta.data)) {
      "data.set"
    } else if ("cell.type.ident" %in% colnames(seurat_obj@meta.data)) {
      "cell.type.ident"
    } else {
      "tree.ident"}
  }

  printTreats <- reactive({
    seurat_obj <- SelectDataset()
    print(seurat_obj)
      if (whichDataset() == "data.set") {
        sort(unique(seurat_obj@meta.data$data.set))
      } else {
        NULL # single data set
    }
  })

  printIdents <- reactive({
    seurat_obj <- SelectDataset()
    print(seurat_obj)
    if ("cell.type.ident" %in% colnames(seurat_obj@meta.data)) {
      sort(unique(seurat_obj@meta.data$cell.type.ident))
    } else {
      sort(unique(seurat_obj@meta.data$tree.ident))
    }
  })

  # returns the correct ID class for cell subset
  IDtype <- function() {
    seurat_obj <- SelectDataset()
    if ("cell.type.ident" %in% colnames(seurat_obj@meta.data)) {
      seurat_obj@meta.data$cell.type.ident
    } else {
      seurat_obj@meta.data$tree.ident
    }
  }


  # ======== Gene Database ======== #
  GeneDB <- function() {
    seurat_obj <- SelectDataset()
    selected <- unlist(strsplit(input$dbGenes, " "))
    
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


  # ======== UMAP Cluster plot ======== #
  ClusterPlotF <- function() {
    seurat_obj <- SelectDataset()

    if ("cell.type.ident" %in% colnames(seurat_obj@meta.data)) {
      umap_clusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.10,
        label = TRUE, label.size = 0, group.by = "cell.type.ident",
        cols = cluster_clrs)
    } else {
      umap_clusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.10,
        label = TRUE, label.size = 0)
    }

    umap_clusters <- umap_clusters + labs(x = "UMAP 1", y = "UMAP 2") + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.title = element_text(size = 12))
    umap_clusters
  }
  
  output$myClusterPlotF <- renderPlot({ClusterPlotF()})
  output$plot.uiClusterPlotF <- renderUI({plotOutput("myClusterPlotF",
    width = "600px", height = "500px")})

  DatasetPlotF <- function() {
    seurat_obj <- SelectDataset()

    umap_dataset <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.10,
      label = TRUE, label.size = 0, group.by = "data.set", cols = trt_colors)

    umap_dataset <- umap_dataset + labs(x = "UMAP 1", y = "UMAP 2") + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.title = element_text(size = 12))
    umap_dataset
  }

  output$myDatasetPlotF <- renderPlot({DatasetPlotF()})
  output$plot.uiDatasetPlotF <- renderUI({plotOutput("myDatasetPlotF",
    width = "600px", height = "500px")})


  # ======== Cluster/Data UMAP ======== #
  DatFeatPlotF <- function() {

    seurat_obj <- SelectDataset()

    if ("cell.type.ident" %in% colnames(seurat_obj@meta.data)) {
      umap_clusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.10,
        label = TRUE, label.size = 0, group.by = "cell.type.ident",
        cols = cluster_clrs)
    } else {
      umap_clusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.10,
        label = TRUE, label.size = 0)
    }

    umap_clusters <- umap_clusters + labs(x = "UMAP 1", y = "UMAP 2") + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.title = element_text(size = 12), legend.position="bottom")

    umap_dataset <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.10,
      label = TRUE, label.size = 0, group.by = "data.set", cols = trt_colors)

    umap_dataset <- umap_dataset + labs(x = "UMAP 1", y = "UMAP 2") + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.title = element_text(size = 12), legend.position = "bottom")

    datfeat_list <- list(umap_clusters, umap_dataset)
    plot_h <- plot_grid(plotlist = datfeat_list, ncol = 2)

    plot_v <- plot_grid(plotlist = datfeat_list, ncol = 1)
    plot_v <- plot_grid(plot_v) + theme(
      plot.background = element_rect(size = 2, color = "#DCDCDC"))
    
    datfeat_list <- list(plot_h, plot_v)
    datfeat_list
  }

  output$myDatFeatPlotH1 <- renderPlot({DatFeatPlotF()[[1]]})
  output$plot.uiDatFeatPlotH1 <- renderUI({
    plotOutput("myDatFeatPlotH1", width = "850px", height = "450px")
  })

  n_panels <- 1:6

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
    seurat_obj <- SelectDataset()
    selected <- unlist(strsplit(input$featureGenes, " "))
    
    ifelse(selected %in% com_name,
      selected <- selected[selected %in% com_name],
    
      ifelse(selected %in% ens_id,
        selected <- gene_df[ens_id %in% selected, 3],"")
    )

    cells_to_plt <- rownames(seurat_obj@meta.data[
      seurat_obj@meta.data$data.set %in% input$cellIdentsFeat,])
    
    feat <- FeaturePlot(seurat_obj[,cells_to_plt], selected,
      reduction = "umap", cols = c(input$CellBackCol, input$CellForeCol),
      combine = FALSE, pt.size = input$ptSizeFeature)

    for(k in 1:length(feat)) {
      feat[[k]] <- feat[[k]] + labs(x = "UMAP 1", y = "UMAP 2") + 
      theme(axis.text.x = element_blank(), legend.position="none",
      axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.y = element_blank(), axis.title = element_text(size = 18),
      panel.border = element_rect(colour = "#FFFFFF", fill = NA, size = 1))
    }
  # return(plot_grid(plotlist = feat, ncol = 1))

  pg <- plot_grid(plotlist = feat, ncol = 1) +
      labs(title = paste("Selected analysis:",
        as.character(input$Analysis)), subtitle = "", caption = "") +
        theme(plot.title = element_text(face = "bold", size = 15, hjust = 0))
    
  return(pg)
  })

  output$cellSelectFeat <- renderUI({
    pickerInput("cellIdentsFeat", "Add or remove treatments from plot:",
      choices = as.character(printTreats()), multiple = TRUE,
      selected = as.character(printTreats()), options = list(
       `actions-box` = TRUE), width = "85%")
  })

  mismatchFeat <- function() {
    selected <- unlist(strsplit(input$featureGenes, " "))

    mismatch <- ifelse(!selected %in% c(com_name, ens_id),
      selected[!selected %in% c(com_name, ens_id)],"")
    return(mismatch)
  }

  output$notInFeat <- renderText({input$runFeatPlot
    isolate({mismatchFeat()})
  })

  output$SelectedDataFeat <- renderText({input$runFeatPlot
    isolate({input$Analysis})
  })

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

  output$downloadFeaturePlotF <- downloadHandler(
    filename = "Feature_plot.png", content = function(file) {
      png(file, units = "in", res = as.numeric(input$featDPI),
        width = 12, height = 12 * getLenInput(input$featureGenes))
      print(FeaturePlotF())
      dev.off()
    }
  )

  # ======== Violin Plot ======== #
  VlnPlotF <- reactive({
    seurat_obj <- SelectDataset()
    selected <- unlist(strsplit(input$vlnGenes, " "))
    
    ifelse(selected %in% com_name,
      selected <- selected[selected %in% com_name],
    
      ifelse(selected %in% ens_id,
        selected <- gene_df[ens_id %in% selected, 3],"")
    )

    seurat_obj <- seurat_obj[,IDtype() %in% input$cellIdentsVln]

    g <- VlnPlot(seurat_obj, selected,
      pt.size = input$ptSizeVln, combine = FALSE,
      group.by = input$selectGrpVln, cols = cluster_clrs)

    for(k in 1:length(g)) {
      g[[k]] <- g[[k]] + theme(legend.position = "none")
    }

    pg <- plot_grid(plotlist = g, ncol = 1) +
      labs(title = paste("Selected analysis:",
        as.character(input$Analysis)), subtitle = "", caption = "") +
        theme(plot.title = element_text(face = "bold", size = 15, hjust = 0))
    
    return(pg)
  })

  output$cellSelectVln <- renderUI({ # New cell type select
    pickerInput("cellIdentsVln", "Add or remove clusters:",
      choices = as.character(printIdents()), multiple = TRUE,
      selected = as.character(printIdents()), options = list(
       `actions-box` = TRUE), width = "85%")
  })

  mismatchVln <- function() {
    selected <- unlist(strsplit(input$vlnGenes, " "))

    mismatch <- ifelse(!selected %in% c(com_name,ens_id),
      selected[!selected %in% c(com_name,ens_id)],"")
    return(mismatch)
  }

  output$notInVln <- renderText({input$runVlnPlot
    isolate({mismatchVln()})
  })

  output$SelectedDataVln <- renderText({input$runVlnPlot
      isolate({input$Analysis})
  })

  output$myVlnPlotF <- renderPlot({input$runVlnPlot
    isolate({withProgress({p <- VlnPlotF(); print(p)},
      message = "Rendering plot..",
      min = 0, max = 10, value = 10)
    })
  })

  getHeightVln <- function() {
    l <- getLenInput(input$vlnGenes)
    if (l == 1) {h <- "600px"
    } else {
      h <- as.character(ceiling(l) * 600)
      h <- paste0(h, "px")
    }
    return(h)
  }

  output$plot.uiVlnPlotF <- renderUI({input$runVlnPlot
    isolate({h <- getHeightVln(); plotOutput("myVlnPlotF",
      width = "800px", height = h)})
  })

  output$downloadVlnPlot <- downloadHandler(
    filename = "Violin_plot.pdf", content = function(file) {
      pdf(file, onefile = FALSE,
        width = 12,
        height = 10 * getLenInput(input$vlnGenes))
      print(VlnPlotF())
      dev.off()
    }
  )


  # ======== Ridge Plot ======== #
  RdgPlotF <- reactive({
    seurat_obj <- SelectDataset()
    selected <- unlist(strsplit(input$rdgGenes, " "))
    
    ifelse(selected %in% com_name,
      selected <- selected[selected %in% com_name],
    
      ifelse(selected %in% ens_id,
        selected <- gene_df[ens_id %in% selected, 3],"")
    )

    seurat_obj <- seurat_obj[,IDtype() %in% input$cellIdentsRdg]
  
    g <- RidgePlot(seurat_obj, selected, combine = FALSE,
      group.by = input$selectGrpRdg, cols = cluster_clrs)

    for(k in 1:length(g)) {
      g[[k]] <- g[[k]] + theme(legend.position = "none")
    }

    # return(plot_grid(plotlist = g, ncol = 1))
    
    pg <- plot_grid(plotlist = g, ncol = 1) +
      labs(title = paste("Selected analysis:",
        as.character(input$Analysis)), subtitle = "", caption = "") +
        theme(plot.title = element_text(face = "bold", size = 15, hjust = 0))
    
    return(pg)
  })

  output$cellSelectRdg <- renderUI({ # New cell type select
    pickerInput("cellIdentsRdg", "Add or remove clusters:",
      choices = as.character(printIdents()), multiple = TRUE,
      selected = as.character(printIdents()), options = list(
       `actions-box` = TRUE), width = "85%")
  })

  mismatchRdg <- function() {
    selected <- unlist(strsplit(input$rdgGenes, " "))

    mismatch <- ifelse(!selected %in% c(com_name,ens_id),
      selected[!selected %in% c(com_name,ens_id)],"")
    return(mismatch)
  }

  output$notInRdg <- renderText({input$runRdgPlot
    isolate({mismatchRdg()})
  })

  output$SelectedDataRdg <- renderText({input$runRdgPlot
    isolate({input$Analysis})
  })

  output$myRdgPlotF <- renderPlot({input$runRdgPlot
    isolate({withProgress({p <- RdgPlotF(); print(p)},
      message = "Rendering plot..", min = 0, max = 10, value = 10)})
  })

  getHeightRdg <- function() {
    l <- getLenInput(input$rdgGenes)
    if (l == 1) {h <- "600px"
    } else {
      h <- as.character(ceiling(l) * 600)
      h <- paste0(h, "px")
    }
    return(h)
  }

  output$plot.uiRdgPlotF <- renderUI({input$runRdgPlot
    isolate({h <- getHeightRdg(); plotOutput("myRdgPlotF",
      width = "800px", height = h)})
  })

  output$downloadRdgPlot <- downloadHandler(
    filename = "Ridge_plot.pdf", content = function(file) {
      pdf(file, onefile = FALSE,
        width = 12,
        height = 10 * getLenInput(input$rdgGenes))
      print(RdgPlotF())
      dev.off()
    }
  )


  # ======== Dot Plot ======== #
  DotPlotF <- reactive({
    clustering <- input$dPlotClust
    if (clustering == TRUE) {
          seurat_obj <- SelectDataset()
      selected <- unlist(strsplit(input$dotGenes, " "))
        
      ifelse(selected %in% com_name,
        selected <- selected[selected %in% com_name],
      
        ifelse(selected %in% ens_id,
          selected <- gene_df[ens_id %in% selected, 3],"")
      )
      
      seurat_obj <- seurat_obj[,IDtype() %in% input$cellIdentsDot]

      seurat_obj_sub <- seurat_obj[rownames(seurat_obj) %in% selected,]
      dist_mat <- dist(seurat_obj_sub@assays$RNA@data)
      clust <- hclust(dist_mat)
      markers_clust <- clust$labels
      
      g <- DotPlot(seurat_obj, features = markers_clust,
        cols = "RdYlBu", dot.scale = input$dotScale,
        group.by = input$selectGrpDot)
        
      g <- g + labs(title = paste("Selected analysis:",
        as.character(input$Analysis)), subtitle = "", caption = "") +
        theme(plot.title = element_text(face = "plain", size = 14))

      g <- g + coord_flip() + theme(
        axis.text.x = element_text(angle = 90, hjust = 1))

    } else {
      seurat_obj <- SelectDataset()
      selected <- unlist(strsplit(input$dotGenes, " "))
      
      ifelse(selected %in% com_name,
        selected <- selected[selected %in% com_name],
      
        ifelse(selected %in% ens_id,
          selected <- gene_df[ens_id %in% selected, 3],"")
        )

      seurat_obj <- seurat_obj[,IDtype() %in% input$cellIdentsDot]
      print(input$cellIdentsDot)

      g <- DotPlot(seurat_obj, features = selected,
        cols = "RdYlBu", dot.scale = input$dotScale,
        group.by = input$selectGrpDot)

      g <- g + labs(title = paste("Selected analysis:",
        as.character(input$Analysis)), subtitle = "", caption = "") +
        theme(plot.title = element_text(face = "plain", size = 14))

      g <- g + coord_flip() + theme(
        axis.text.x = element_text(angle = 90, hjust = 1))
    }
    return(g)
  })

  output$cellSelectDot <- renderUI({ # New cell type select
    pickerInput("cellIdentsDot", "Add or remove clusters:",
      choices = as.character(printIdents()), multiple = TRUE,
      selected = as.character(printIdents()), options = list(
       `actions-box` = TRUE), width = "85%")
  })

  mismatchDot <- function() {
    selected <- unlist(strsplit(input$dotGenes, " "))

    mismatch <- ifelse(!selected %in% c(com_name,ens_id),
      selected[!selected %in% c(com_name,ens_id)],"")
    return(mismatch)
  }

  output$notInDot <- renderText({input$runDotPlot
    isolate({mismatchDot()})
  })

  output$SelectedDataDot <- renderText({input$runDotPlot
      isolate({input$Analysis})
  })

  output$myDotPlotF <- renderPlot({input$runDotPlot
    isolate({withProgress({p <- DotPlotF(); print(p)},
      message = "Rendering plot..", min = 0, max = 10, value = 10)
    })
  })

  getHeightDot <- function() {
    l <- getLenInput(input$dotGenes)
    h <- paste0(as.character(l * 35), "px")
    return(h)
  }
  
  # ! check/change for project
  # TODO create formula for n clusters/treats and dplot width
  dplotWidth <- function () {
    if(input$selectGrpDot == "data.set") {
      w <- "500px"
    } else {
      w <- "800px"
    }
  }

  dotHeight <- function() {
    l <- getLenInput(input$dotGenes)
    l <- as.numeric(l)
    return(l)
  }

  # output$plot.uiDotPlotF <- renderUI({input$runDotPlot
  #   isolate({h <- getHeightDot(); plotOutput("myDotPlotF",
  #     width = dplotWidth(), height = h)
  #   })
  # })

  output$plot.uiDotPlotF <- renderUI({input$runDotPlot
  isolate({h <- getHeightDot(); plotOutput("myDotPlotF",
    width = paste0(input$manAdjustDotW, "px"),
    height = paste0(input$manAdjustDotH, "px"))})
  })
  
  output$downloadDotPlot <- downloadHandler(
    filename = "dot_plot.pdf", content = function(file) {
      pdf(file, onefile = FALSE, width = 12, height = dotHeight() * 0.5)
      print(DotPlotF())
      dev.off()
    }
  )

  # ======== pHeatmap ======== #
  selectedCellsHmap <- reactive({
    multiGrep2(input$cellIdentsHmap, colnames(hmap_list[[1]]))
  })

  pHeatmapF <- reactive({
    selected <- unlist(strsplit(input$PhmapGenes, " "))

    ifelse(selected %in% com_name,
      selected <- selected[selected %in% com_name],
    
      ifelse(selected %in% ens_id,
        selected <- gene_df[ens_id %in% selected, 3],"")
    )

    avg_mtx <- hmap_list[[input$mtxSelectHmap]]
    goi_mat <- avg_mtx[rownames(avg_mtx) %in% selected, selectedCellsHmap()]
    goi_mat <- goi_mat[match(selected, rownames(goi_mat)),]

    n_trt <- length(unique(file_list[[1]]@meta.data$data.set))
    mtx_cols <- ncol(goi_mat) - n_trt

    hmapColors <- colorRampPalette(
      rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)

    pheatmap::pheatmap(goi_mat, cluster_rows = input$pHmapClust,
      cluster_cols = FALSE, color = hmapColors, annotation_col = NULL,
      legend = FALSE, annotation_colors = anno_cols,
      gaps_col = seq(n_trt, mtx_cols, by = n_trt),
      annotation_names_col = FALSE, annotation_legend = FALSE)
  })

  mismatchPhmap <- function() {
    selected <- unlist(strsplit(input$PhmapGenes, " "))

    mismatch <- ifelse(!selected %in% c(com_name, ens_id),
      selected[!selected %in% c(com_name, ens_id)],"")
    return(mismatch)
  }

  output$notInPhmap <- renderText({input$runPhmap
    isolate({mismatchPhmap()})
  })

  output$SelectedDataPhmap <- renderText({input$runPhmap
    isolate({input$Analysis})
  })

  avg_mtx_names <- unique(unlist(lapply(seq_along(colnames(hmap_list[[1]])),
    function(i){strsplit(colnames(hmap_list[[1]]), "_")[[i]][1]})))

  output$cellSelectHmap <- renderUI({ # New cell type selected
  pickerInput("cellIdentsHmap", "Add or remove clusters:",
    choices = avg_mtx_names, multiple = TRUE,
    selected = avg_mtx_names, options = list(
      `actions-box` = TRUE), width = "80%")
  })

  output$myPhmapF <- renderPlot({input$runPhmap
    isolate({withProgress({p <- pHeatmapF(); print(p)},
      message = "Rendering plot..", min = 0, max = 10, value = 10)
    })
  })

  getWidthPhmap <- reactive({
    if(input$pHmapClust == TRUE ) {
      w <- (length(selectedCellsHmap()) * 14) + 150
      return(w)
    } else {
      w <- (length(selectedCellsHmap()) * 14) + 90
      return(w)
    }
  })

  getHeightPhmap <- reactive({
    l <- getLenInput(input$PhmapGenes)
    h <- (l * 13) + 85
    return(h)
  })

  output$plot.uiPheatmapF <- renderUI({input$runPhmap
    isolate({
      w <- paste0(getWidthPhmap(), "px"); h <- paste0(getHeightPhmap(), "px")
        plotOutput("myPhmapF", width = w, height = h)
    })
  })

  output$downloadPhmap <- downloadHandler(
    filename = "heatmap.pdf", content = function(file) {
      pdf(file, width = (getWidthPhmap() / 90),
          height = (getHeightPhmap() / 90))
      print(pHeatmapF())
      dev.off()
    }
  )


  # ======== Differential Expression ======== #
  diffExp <- reactive({
    seurat_obj <- SelectDataset()
    seurat_obj <- seurat_obj[,IDtype() %in% input$cellIdentsDiff]
    meta <- seurat_obj@meta.data

    print(input$identText1)
    print(input$identText1)
    subset1 <- as.character(input$identText1)
    subset2 <- as.character(input$identText2)

    if ("data.set" %in% colnames(meta)) {
      group1 <- rownames(meta[meta$data.set %in% subset1,])
      group2 <- rownames(meta[meta$data.set %in% subset2,])
    } else {
      group1 <- rownames(meta[meta$cell.type.ident %in% subset1,])
      group2 <- rownames(meta[meta$cell.type.ident %in% subset2,])
    }

    diff_results <- FindMarkers(test.use = input$statSelectDiff,
      seurat_obj, ident.1 = group1, ident.2 = group2)

    diff_results$Gene.name.uniq <- ""
    diff_results$Gene.name.uniq <- rownames(diff_results)

    pval <- as.numeric(input$pValCutoff)
    diff_results <- diff_results[
      diff_results$p_val_adj < pval, c(6,1:5)]
    diff_results <<- diff_results[
      order(diff_results$avg_logFC, decreasing = TRUE),]
  })

  # Requires input$identText to execute before diffExp()
  diffReact <- eventReactive(c(input$identText1, input$identText2), diffExp())

  output$diffTable <- renderTable({input$runDiffExp
    isolate({withProgress(diffReact(), message = "Calculating..",
      min = 0, max = 10, value = 10)}
  )}, digits = -5)

  output$diffOut1 <- renderUI({
    pickerInput("identText1", tags$b("Group 1 - positive FC"),
      choices = as.character(printTreats()), multiple = TRUE,
      selected = as.character(printTreats())[1], options = list(
       `actions-box` = TRUE), width = "80%")
  })

  output$diffOut2 <- renderUI({
    pickerInput("identText2", tags$b("Group 2 - negative FC"),
      choices = as.character(printTreats()), multiple = TRUE,
      selected = as.character(printTreats())[2],options = list(
        `actions-box` = TRUE), width = "80%")
  })

  output$cellSelectDiff <- renderUI({ # New cell type select
    pickerInput("cellIdentsDiff", "Add or remove clusters:",
      choices = as.character(printIdents()), multiple = TRUE,
      selected = as.character(printIdents()), options = list(
       `actions-box` = TRUE), width = "80%")
  })

  output$SelectedDataDiff <- renderText({input$runDiffExp
    isolate({input$Analysis})
  })

  makeDiffTable <- function() {
    markerTable <<- inner_join(diff_results,
      gene_df, by = "Gene.name.uniq")
    return(markerTable)
  }

  # qmethod "double" rids the table of escape backslashes
  output$downloadDiffExp <- downloadHandler(
    filename = "diff_exp_results.tsv",
    content = function(file) {
      write.table(makeDiffTable(), file,
        row.names = FALSE, col.names = TRUE,
        qmethod = "double", sep = "\t")
    }
  )


  # ======== Download meta data ======== #
  output$downloadClusterMarkers <- downloadHandler(
    filename = "cluster_markers.xlsx",
    content = function(file) {
      file.copy("./data/cluster_markers.xlsx", file)
    }
  )
} # Server close






