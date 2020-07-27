library(shiny)
library(cowplot)
library(Seurat)
library(ggplot2)
library(shinythemes)
library(shinyWidgets)
library(dplyr)


# ============ Functions
multiGrep2 <- function(toMatch, toSearch, ...) {
  toMatch <- ifelse(grepl("*", toMatch),
    gsub("\\*","\\\\*", toMatch), toMatch <- toMatch)
  
  toMatch <- paste(toMatch, collapse = "|")
  inCommon <- grep(toMatch, toSearch, value = FALSE)
  return(inCommon)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

getLenInput <- function(input) {
  selected <- unlist(strsplit(input, " "))
  
  ifelse(selected %in% com_name,
    selected <- gene_df[com_name %in% selected, 3],
    
    ifelse(selected %in% ens_id,
      selected <- gene_df[ens_id %in% selected, 3],"")
  )
  len <- length(selected)
  return(len)
}

# ============ Reading in data
files <- list.files("./data", pattern = "TRIMMED", full.names = TRUE)

file_list <- list()
for (i in 1:length(files)) {
    file_list[[i]] <- readRDS(files[i])
    DefaultAssay(file_list[[i]]) <- "RNA"
}

# ! check/change for project
Idents(file_list[[1]]) <- file_list[[1]]@meta.data$cell.type.ident
seurat_obj <- file_list[[1]]

print(object.size(file_list), units = "MB")
names(file_list) <- as.character(c("All cells"))

gene_df <- read.table("./data/Danio_Features_unique_Ens91_v2.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE)

gene_df <- gene_df[gene_df$Gene.name.uniq %in% rownames(seurat_obj),]

ens_id <- gene_df$Gene.stable.ID
com_name <- gene_df$Gene.name.uniq

gene_df$in.dataset <- ""
gene_df$in.dataset <- (gene_df$Gene.name.uniq %in% rownames(seurat_obj))
gene_df <- gene_df[,c(1:3,6,4:5)]
gene_df$ZFIN.ID <- paste0("=HYPERLINK(", '"', gene_df$ZFIN.ID, '"',")")

# ! check/change for project
avg_mtx <- readRDS(paste0("./data/mtx_CLR_nrml",
  "_scld_tmpts_in_cell_type_mphg_regen_14k_cells_seurat3_v1.0_.RDS"))

# rearrange hmap to Tatjana's order
colnames(avg_mtx) <- sub(" ", "",colnames(avg_mtx))
to_match <- c("kng1", "actn3b", "gata3", "rbp4", "cldnh",
  "spock3", "eomesa", "mpx", "pcna", "fabp3", "f13a1b", "mcamb", "tspan10",
  "runx3", "irg1", "stat1b")

new_order <- unlist(lapply(1:16,
  function(i){grep(to_match[i], colnames(avg_mtx))}))
avg_mtx <- avg_mtx[,new_order]

clusterColors <- gg_color_hue(
  length(levels(seurat_obj@active.ident)))

# ! check/change for project
trt_colors <- c("green3", "darkorange",
  "deeppink", "mediumorchid1", "deepskyblue")


# ================================== Server ===================================


server <- function(input, output) {

  # ======== Dataset selection ======== #
  SelectDataset <- reactive({
    seurat_obj <- file_list[[input$DataSet]]
    print(names(file_list[input$DataSet]))
    return(seurat_obj)
  })

  whichDataset <- function() {
  seurat_obj <- SelectDataset()
  if ("data.set" %in% colnames(seurat_obj@meta.data)) {
    "data.set"} else {"cell.type.ident"}
  }

  printTreats <- reactive({
    seurat_obj <- SelectDataset()
    print(seurat_obj)
      if (whichDataset() == "data.set") {
        sort(unique(seurat_obj@meta.data$data.set))
      } else {
        sort(unique(seurat_obj@meta.data$cell.type.ident))
    }
  })

  printIdents <- reactive({
    seurat_obj <- SelectDataset()
    print(seurat_obj)
      if (whichDataset() == "cell.type.ident") {
        sort(unique(seurat_obj@meta.data$data.set))
      } else {
        sort(unique(seurat_obj@meta.data$cell.type.ident))
    }
  })

  sliceDataset <- function() {
    seurat_obj <- seurat_obj[,
      seurat_obj@meta.data$cell.type.ident %in% input$CellTypeSelect]
  }

  # ======== Gene Database ======== #
  GeneDB <- function() {
    seurat_obj <- SelectDataset()
    selected <- unlist(strsplit(input$dbGenes, " "))

    ifelse(selected %in% gene_df$Gene.name.uniq,
      ind <- multiGrep2(selected, gene_df$Gene.name.uniq),

      ifelse(selected %in% gene_df$Gene.stable.ID,
        ind <- multiGrep2(selected, gene_df$Gene.stable.ID),
        "gene not in database")
    )
    gene_df[ind,]
  }
  output$GeneDB <- renderTable({GeneDB()})


  # ======== UMAP Cluster plot ======== #
  ClusterPlotF <- function() {
    seurat_obj <- SelectDataset()

    umap_clusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.10,
      label = TRUE, label.size = 4)

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
    
    umap_clusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.10,
    label = TRUE, label.size = 0)

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
  output$plot.uiDatFeatPlotH1 <- renderUI(
    {plotOutput("myDatFeatPlotH1", width = "850px", height = "450px")})

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
  FeaturePlotF <- function() {
    seurat_obj <- SelectDataset()
    selected <- unlist(strsplit(input$featureGenes, " "))
    
    ifelse(selected %in% com_name,
      selected <- selected[selected %in% com_name],
    
      ifelse(selected %in% ens_id,
        selected <- gene_df[ens_id %in% selected, 3],"")
    )

    cells_to_plt <- rownames(seurat_obj@meta.data[
      seurat_obj@meta.data$data.set %in% input$cellIdentsFeat, ])
    
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
  return(plot_grid(plotlist = feat, ncol = 1))
  }
# 
  output$cellSelectFeat <- renderUI({
    pickerInput("cellIdentsFeat", "Add/remove treatments from plot:",
      choices = as.character(printTreats()), multiple = TRUE,
      selected = as.character(printTreats()), options = list(
       `actions-box` = TRUE), width = "85%")
  })

  mismatchFeat <- function() {
    selected <- unlist(strsplit(input$featureGenes, " "))

    mismatch <- ifelse(!selected %in% c(com_name,ens_id),
      selected[!selected %in% c(com_name,ens_id)],"")
    return(mismatch)
  }

  output$notInFeat <- renderText({input$runFeatPlot
    isolate({mismatchFeat()})
  })

  output$SelectedDataFeat <- renderText({input$runFeatPlot
    isolate({input$DataSet})
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
  VlnPlotF <- function() {
    seurat_obj <- SelectDataset()
    selected <- unlist(strsplit(input$vlnGenes, " "))
    
    ifelse(selected %in% com_name,
      selected <- selected[selected %in% com_name],
    
      ifelse(selected %in% ens_id,
        selected <- gene_df[ens_id %in% selected, 3],"")
    )

    seurat_obj <- seurat_obj[, # New
      seurat_obj@meta.data$cell.type.ident %in% input$cellIdentsVln]

    g <- VlnPlot(seurat_obj, selected,
      pt.size = input$ptSizeVln, combine = FALSE,
      group.by = input$selectGrpVln, cols = clusterColors)

    for(k in 1:length(g)) {
      g[[k]] <- g[[k]] + theme(legend.position = "none")
    }
    return(plot_grid(plotlist = g, ncol = 1))
  }

  output$cellSelectVln <- renderUI({ # New cell type select
    pickerInput("cellIdentsVln", "Add/remove cell types:",
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
      isolate({input$DataSet})
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
  RdgPlotF <- function() {
    seurat_obj <- SelectDataset()
    selected <- unlist(strsplit(input$rdgGenes, " "))
    
    ifelse(selected %in% com_name,
      selected <- selected[selected %in% com_name],
    
      ifelse(selected %in% ens_id,
        selected <- gene_df[ens_id %in% selected, 3],"")
    )

    seurat_obj <- seurat_obj[, # New
      seurat_obj@meta.data$cell.type.ident %in% input$cellIdentsRdg]
  
    g <- RidgePlot(seurat_obj, selected, combine = FALSE,
      group.by = input$selectGrpRdg, cols = clusterColors)

    for(k in 1:length(g)) {
      g[[k]] <- g[[k]] + theme(legend.position = "none")
    }

    return(plot_grid(plotlist = g, ncol = 1))
  }

  output$cellSelectRdg <- renderUI({ # New cell type select
    pickerInput("cellIdentsRdg", "Add/remove cell types:",
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
    isolate({input$DataSet})
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
  DotPlotF <- function() {
    clustering <- input$dPlotClust
    if (clustering == TRUE) {
          seurat_obj <- SelectDataset()
      selected <- unlist(strsplit(input$dotGenes, " "))
        
      ifelse(selected %in% com_name,
        selected <- selected[selected %in% com_name],
      
        ifelse(selected %in% ens_id,
          selected <- gene_df[ens_id %in% selected, 3],"")
      )
      
      seurat_obj <- seurat_obj[, # New
        seurat_obj@meta.data$cell.type.ident %in% input$cellIdentsDot]

      seurat_obj_sub <- seurat_obj[rownames(seurat_obj) %in% selected,]
      dist_mat <- dist(seurat_obj_sub@assays$RNA@data)
      clust <- hclust(dist_mat)
      markers_clust <- clust$labels
      
      g <- DotPlot(seurat_obj, features = markers_clust,
        cols = "RdYlBu", dot.scale = input$dotScale,
        group.by = input$selectGrpDot)
      g <- g + coord_flip() + theme(
        axis.text.x = element_text(angle = 90, hjust = 1))
      # g <- g + ggtitle(as.character(input$DataSet)) #

    } else { 
      seurat_obj <- SelectDataset()
      selected <- unlist(strsplit(input$dotGenes, " "))
      
      ifelse(selected %in% com_name,
        selected <- selected[selected %in% com_name],
      
        ifelse(selected %in% ens_id,
          selected <- gene_df[ens_id %in% selected, 3],"")
        )

      seurat_obj <- seurat_obj[, # New
        seurat_obj@meta.data$cell.type.ident %in% input$cellIdentsDot]

      g <- DotPlot(seurat_obj, features = selected,
        cols = "RdYlBu", dot.scale = input$dotScale,
        group.by = input$selectGrpDot)
      g <- g + coord_flip() + theme(
        axis.text.x = element_text(angle = 90, hjust = 1))
      # g <- g + ggtitle(as.character(input$DataSet))
    }
    return(g)
  }

  output$cellSelectDot <- renderUI({ # New cell type select
    pickerInput("cellIdentsDot", "Add/remove cell types:",
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
      isolate({input$DataSet})
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

  dplotWidth <- function () {
    if(input$selectGrpDot == "data.set") {
      w <- "450px"
    } else {
      w <- "800px"
    }
  }

  output$plot.uiDotPlotF <- renderUI({input$runDotPlot
    isolate({h <- getHeightDot(); plotOutput("myDotPlotF",
      width = dplotWidth(), height = h)
    })
  })

  dotHeight <- function() {
    l <- getLenInput(input$dotGenes)
    l <- as.numeric(l)
    return(l)
  }
  
  output$downloadDotPlot <- downloadHandler(
    filename = "dot_plot.pdf", content = function(file) {
      pdf(file, onefile = FALSE, width = 12, height = dotHeight() * 0.5)
      print(DotPlotF())
      dev.off()
    }
  )


  # ======== pHeatmap ======== #
  selectedCellsHmap <- reactive({
    multiGrep2(input$cellIdentsHmap, colnames(avg_mtx))
  })

  pHeatmapF <- function() {
    selected <- unlist(strsplit(input$PhmapGenes, " "))

    ifelse(selected %in% com_name,
      selected <- selected[selected %in% com_name],
    
      ifelse(selected %in% ens_id,
        selected <- gene_df[ens_id %in% selected, 3],"")
    )
    
    goi_mat <- avg_mtx[rownames(avg_mtx) %in% selected, selectedCellsHmap()]
    n_trt <- length(unique(seurat_obj@meta.data$data.set))
    mtx_cols <- ncol(avg_mtx) - n_trt

    n_trt <- length(unique(file_list[[1]]@meta.data$data.set))
    mtx_cols <- ncol(goi_mat) - n_trt

    pheatmap::pheatmap(goi_mat, cluster_rows = input$pHmapClust,
      cluster_cols = FALSE, annotation_col = NULL,
      legend = FALSE, annotation_colors = anno_cols,
      gaps_col = seq(n_trt, mtx_cols, by = n_trt),
      annotation_names_col = FALSE, annotation_legend = FALSE)
  }

  avg_mtx_names <- unique(unlist(lapply(seq_along(colnames(avg_mtx)),
    function(i){strsplit(colnames(avg_mtx), "_")[[i]][1]})))

  output$cellSelectHmap <- renderUI({ # New cell type selected
  pickerInput("cellIdentsHmap", "Add or remove clusters:",
    choices = avg_mtx_names, multiple = TRUE,
    selected = avg_mtx_names, options = list(
      `actions-box` = TRUE), width = "80%")
  })

  mismatchPhmap <- function() {
    selected <- unlist(strsplit(input$PhmapGenes, " "))

    mismatch <- ifelse(!selected %in% c(com_name,ens_id),
      selected[!selected %in% c(com_name,ens_id)],"")
    return(mismatch)
  }

  output$notInPhmap <- renderText({input$runPhmap
    isolate({mismatchPhmap()})
  })

  output$SelectedDataPhmap <- renderText({input$runPhmap
      isolate({input$DataSet})
  })

  output$myPhmapF <- renderPlot({input$runPhmap
    isolate({withProgress({p <- pHeatmapF(); print(p)},
      message = "Rendering plot..", min = 0, max = 10, value = 10)
    })
  })

  getHeightPhmap <- function() {
    l <- getLenInput(input$PhmapGenes)
    h <- paste0(as.character((l * 13) + 85), "px")
    return(h)
  }

  output$plot.uiPheatmapF <- renderUI({input$runPhmap
    isolate({h <- getHeightPhmap(); plotOutput("myPhmapF",
      width = "1125px", height = h)
    })
  })

  pHmapHeight <- function() {
    l <- getLenInput(input$PhmapGenes)
    l <- as.numeric(l)
    print(l)
    return(l)
  }
  
  output$downloadPhmap <- downloadHandler(
    filename = "heatmap.pdf", content = function(file) {
      pdf(file, width = 12, height = (pHmapHeight() * 0.20) + 1)
      print(pHeatmapF())
      dev.off()
    }
  )


  # ======== Differential Expression ======== #
  diffExp <- function() {
    seurat_obj <- SelectDataset()
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

    diff_results <- FindMarkers(
      seurat_obj, ident.1 = group1, ident.2 = group2)

    diff_results$Gene.name.uniq <- ""
    diff_results$Gene.name.uniq <- rownames(diff_results)

    pval <- as.numeric(input$pValCutoff)
    diff_results <- diff_results[
      diff_results$p_val_adj < pval, c(6,1:5)]
    diff_results <<- diff_results[
      order(diff_results$avg_logFC, decreasing = TRUE),]
  }

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

  output$SelectedDataDiff <- renderText({input$runDiffExp
    isolate({input$DataSet})
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

  output$downloadDatasetGenes <- downloadHandler(
    filename = "genes_in_data_set.xlsx",
    content = function(file) {
      file.copy("./data/genes_in_data_set.xlsx", file)
    }
  )

  output$downloadUMAPtreatment <- downloadHandler(
    filename = "UMAP_treatment.pdf",
    content = function(file) {
      file.copy("./data/UMAP_treatment.pdf", file)
    }
  )

  output$downloadUMAPclusters <- downloadHandler(
    filename = "UMAP_clusters.pdf",
    content = function(file) {
      file.copy("./data/UMAP_clusters.pdf", file)
    }
  )

  output$downloadAllAmbGenes <- downloadHandler(
    filename = "ambiguous_genes_Ens91.xlsx",
    content = function(file) {
      file.copy("./data/ambiguous_genes_Ens91.xlsx", file)
    }
  )

  output$downloadAmbGenesDataset <- downloadHandler(
    filename = "ambiguous_genes_in_dataset.xlsx",
    content = function(file) {
      file.copy("./data/ambiguous_genes_in_dataset.xlsx", file)
    }
  )

  # ! https://www.r-bloggers.com/a-little-trick-for-debugging-shiny/
  # observeEvent(input$browser, {
  #   browser()
  # })
} # Server close


# =============================== User Interface ==============================


ui <- fixedPage(theme = shinytheme("lumen"), # paper lumen cosmo spacelab
  tags$head(includeCSS(paste0("./www/styles.css"))),
  div(headerPanel("Macrophage scRNA-seq"), style = 'width:1560px;'),
  div(tabsetPanel(
    
    # ================ #
    tabPanel("Welcome!", fluid = TRUE,
      mainPanel(class = "welcome",

        fluidRow(tags$br()),
        fluidRow(tags$br()),

        column(12, tags$br()),

        column(12, align = "center",
          uiOutput("plot.uiDatFeatPlotH1"), tags$br()),

        fluidRow(tags$br()),
        fluidRow(tags$br()),
        column(12, align = "center",
          tags$b("Select Data Set"),
          pickerInput("DataSet", label = "",
            choices = list(Combined = names(file_list)),
            selected = "All cell types", width = "50%")
        ),
        fluidRow(tags$br()),
        fluidRow(tags$br()),

        column(12, tags$hr()),
        fluidRow(tags$br()),
        fluidRow(tags$br()),
        column(10, align = "center", offset = 1,
          column(12, align = "left", tags$b("Instructions")),
          column(12, align = "left",
              'All genes available for plotting can be downloaded in the
              Excel spreadsheet below labeled "genes in dataset", using either
              Ensembl IDs or common names from the',
              tags$b("Gene.name.unique"),'column as input. You cannot, however,
              mix common names with Ensembl IDs in the same query. Groups of
              genes can be directly copied/pasted from the spreadsheet into 
              the app input field and will have the necessary spacing 
              by default. Please note that this data set was produced with the
              Ensembl 91 gene annotation in zebrafish (genome version 10).
              We therefore recommend using Ensembl gene IDs as input, 
              as common gene names can change with annotation updates.',
              'Cluster markers and related figures can be downloaded in',
              tags$a(href = "http://bioinfo/n/projects/ddiaz/Analysis/Scripts/sb2191-regen/regen-summary/site/IntegratedData/",
                tags$b("this notebook")),
              '. All genes used for this dataset can be downloaded below:'),

          fluidRow(tags$br()),
          fluidRow(tags$br()),
          column(12, align = "center", offset = 0,
            downloadButton("downloadDatasetGenes", "Genes in Data Set",
              style = "padding:8px; font-size:80%")),
          fluidRow(tags$br()),
          fluidRow(tags$br()),
          
          fluidRow(tags$br()),
          fluidRow(tags$br()),
          column(12, align = "left", tags$b("An important note on ambiguous gene names")),
          column(12, align = "left",
              'Gene expression in this data set is quantfied by the number of
              deduplicated UMIs (unique molecular index) that map to Ensembl 
              gene IDs, which are identify unique coding sequences (CDS) in the 
              zebrafish genome. In some cases different Ensembl IDs
              will have the same common gene name. This may occur
              when there is no consensus over which CDS represents
              the common name, or because the product of a particular
              CDS has not been characterized in detail. These repeated
              common names are denoted by an asterisk followed by an
              integer value (e.g. sox4a*1). The asterisk is not a part of the
              common name by default; it was added to signify that the name
              is repeated in this data set and needs further clarification.
              The integer after the asterisk highlights which occurrence
              of the repeat you are viewing, but does not carry any
              functional significance. For example, sox4a has two
              different Ensembl IDs in version 91 of the annotation,
              ENSDARG00000096389 - sox4, and ENSDARG00000004588 - sox4a*1, but
              only ENSDARG00000004588 - sox4a*1 is expressed in this data set.',
              
              fluidRow(tags$br()),
              fluidRow(tags$br()),

              'The most important item to consider when referencing the
              nucleotide sequence of a gene is',
                tags$b("which Ensembl ID corresponds to 
                your expression pattern of interest."),
              'This ensures that you are targeting the same CDS that reads 
              are mapping to for this particular data set. You can easily 
              check if a gene name is ambiguous by copying any portion
              of the common name into the Gene Database section of the app
              (sox4a is shown as an example by default). Details on
              how Ensembl IDs are curated can be found at the folowing link:',
                tags$a(
                  href = "http://www.ensembl.org/info/genome/genebuild/genome_annotation.html",
                    "http://www.ensembl.org/info/genome/genebuild/genome_annotation.html"),
              '. Additionally, there are two spreadsheets below with all
              of the repeated common names in both this data set, and the
              Ensembl 91 zebrafish annotation.')),
        fluidRow(tags$br()),
        fluidRow(tags$br())
        
        # ! https://www.r-bloggers.com/a-little-trick-for-debugging-shiny/
        # type $('#browser').show(); in browser java console in web browser
        # actionButton("browser", "browser"),
        # tags$script("$('#browser').hide();")
      )
    ),


    # ================ #
    # tabPanel("Gene Database", fluid = TRUE,
    #   sidebarLayout(
    #     sidebarPanel(
    #       textInput("dbGenes", "Insert gene name or ensembl ID:",
    #         value = "gadd45gb.1 slc1a3a znf185 si:ch73-261i21.5"),
    #       fluidRow(tags$br())
    #     ),
    #       mainPanel(fluidRow(
    #           uiOutput("GeneDB")
    #       )
    #     )
    #   )
    # ),


    # ================ #
    tabPanel("Feature Plots", fluid = FALSE,
      sidebarLayout(fluid = TRUE,
        
        sidebarPanel(width = 4,
              
          column(12, align = "left",
            textInput("featureGenes", "Insert gene name or ensembl ID:",
              value = "mpeg1.1 mfap4")),
          
          column(12, align = "center",
            actionButton("runFeatPlot", "Generate Plots",
            style = 'padding:5px; font-size:80%')),
          
          column(12, tags$hr(width = "50%"), align = "center"),
          column(12, align = "center", uiOutput("cellSelectFeat")),

          column(12, tags$br()),
          column(12,
            column(6, textInput("CellBackCol",
              "Cell background color:", value = "azure3")),
            column(6, textInput("CellForeCol",
              "Cell foreground color:", value = "blue3"))
            ),

          column(12, tags$br()),
          column(12, align = "center",
            column(6, align = "left",
              numericInput("featDPI", "Download quality (DPI):",
                value = 200, min = 50, step = 25, max = 400, width = "100%")),
            column(6, align = "left",
              numericInput("ptSizeFeature", "Input cell size:", value = 1.00,
                min = 0.25, step = 0.25, max = 2.00, width = "100%"))
          ),
          
          column(12, tags$br()),
          column(12, align = "center",
            downloadButton("downloadFeaturePlotF", "Feature plot.png",
            style = 'padding:5px; font-size:80%')),

          fluidRow(tags$br()),
          fluidRow(tags$br()),
          fluidRow(tags$br()),
          column(12, uiOutput("plot.uiDatFeatPlotV1"), align = "center"),
          fluidRow(tags$br()),
          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(8, tags$br()),
            column(8, tags$b("Mismatches or genes not present"),
                "(if applicable)", tags$b(":")),
            column(8,uiOutput("notInFeat")),
            column(8, tags$hr()),
            
            fluidRow(tags$br()),
            column(12, uiOutput("plot.uiFeaturePlotF")
              )
          )
        )
      )
    ),
    
    
    # ================ #
    tabPanel("Violin Plots", #fluid = FALSE,
      sidebarLayout(fluid = TRUE,
        
        sidebarPanel(fluid = FALSE, width = 4,
          column(12, textInput("vlnGenes", width = "100%",
            "Insert gene name or ensembl ID:",
            value = "mpeg1.1 mfap4")),
          
          column(12, align = "center",
            actionButton("runVlnPlot", "Generate Plots",
            style = 'padding:5px; font-size:80%')),

          column(12, tags$hr(width = "50%"), align = "center"),
          column(12, align = "center", uiOutput("cellSelectVln")), # New
          
          column(12, tags$br()),
          column(12, align = "center",
            column(6,
              radioGroupButtons("selectGrpVln",
                "Group cells by:", choices = list(Dataset = "data.set",
                  Cluster = "cell.type.ident"), width = "100%")),
            column(6,
              numericInput("ptSizeVln", "Input cell size:", value = 0.25,
              min = 0.00, step = 0.75, max = 2.00, width = "80%"))
          ),

          column(12, tags$br()),
          column(12, align = "center", "Plot download (pdf):"),
          column(12, tags$br()),
          column(12, align = "center", downloadButton(
            "downloadVlnPlot", "Violin plot.pdf",
            style = 'padding:5px; font-size:80%')),

          fluidRow(tags$br()),
          fluidRow(tags$br()),
          fluidRow(tags$br()),
          column(12, uiOutput("plot.uiDatFeatPlotV2"), align = "center"),
          fluidRow(tags$br()),
          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(8, tags$br()),
            column(8, tags$b("Gene mismatches"), "(if present)", tags$b(":")),
            column(8,uiOutput("notInVln")),
            column(8, tags$hr()),
            # column(8, tags$b(uiOutput("SelectedDataVln"))), 
            column(12, uiOutput("plot.uiVlnPlotF")
            )
          )
        )
      )
    ),


    # ================ #
    tabPanel("Ridge Plots", #fluid = FALSE,
      sidebarLayout(fluid = TRUE,
        
        sidebarPanel(fluid = FALSE, width = 4,
          column(12, textInput("rdgGenes", width = "100%",
            "Insert gene name or ensembl ID:",
            value = "mpeg1.1 mfap4")),
          
          column(12, align = "center",
            actionButton("runRdgPlot", "Generate Plots",
            style = 'padding:5px; font-size:80%')),

          column(12, tags$hr(width = "50%"), align = "center"),
          column(12, align = "center", uiOutput("cellSelectRdg")), # New
          
          column(12, tags$br()),
           column(12, align = "center",
            column(6,
              radioGroupButtons("selectGrpRdg",
                "Group cells by:", choices = list(Dataset = "data.set",
                  Cluster = "cell.type.ident"), width = "100%")),
            column(6,
              numericInput("ptSizeRdg", "Input cell size:", value = 0.25,
              min = 0.00, step = 0.75, max = 2.00, width = "80%"))
          ),

          column(12, tags$br()),
          column(12, align = "center", "Plot download (pdf):"),
          column(12, tags$br()),
          column(12, align = "center", downloadButton(
            "downloadRdgPlot", "Ridge plot.pdf",
            style = 'padding:5px; font-size:80%')),
          
          fluidRow(tags$br()),
          fluidRow(tags$br()),
          fluidRow(tags$br()),
          column(12, uiOutput("plot.uiDatFeatPlotV3"), align = "center"),
          fluidRow(tags$br()),
          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(8, tags$br()),
            column(8, tags$b("Gene mismatches"), "(if present)", tags$b(":")),
            column(8,uiOutput("notInRdg")),
            column(8, tags$hr()),
            # column(8, tags$b(uiOutput("SelectedDataRdg"))),
            column(12, uiOutput("plot.uiRdgPlotF")
            )
          )
        )
      )
    ),


    # ================ #
    tabPanel("Dot Plot", #fluid = FALSE,
      sidebarLayout(fluid = TRUE,
        
        sidebarPanel(fluid = FALSE, width = 4,
          column(12, align = "left  ",
          textInput("dotGenes",
            "Insert gene name or ensembl ID:",
            value = "mpeg1.1 mfap4 lcp1 f13a1b tnfa lyz ctss1 txn ccl34b.1"),
          
          checkboxInput("dPlotClust",
            label = "Check box to enable row clustering.", value = FALSE)),

          column(12, align = "center",
            actionButton("runDotPlot", "Generate Plots",
            style = 'padding:5px; font-size:80%')),
          column(12, tags$hr(width = "50%"), align = "center"),
          column(12, align = "center", uiOutput("cellSelectDot")), # New
          
          column(12, tags$br()),
          column(12, align = "center",
            column(6,
              radioGroupButtons("selectGrpDot",
                "Group cells by:", choices = list(Dataset = "data.set",
                  Cluster = "cell.type.ident"), width = "100%")),
            column(6,
              numericInput("dotScale", "Dot diameter:", value = 10, min = 4,
                step = 1, max = 20, width = "80%"), align = "center")
          ),

          column(12, tags$br()),
          column(12, align = "center", "Plot download (pdf):"),
          column(12, tags$br()),
          column(12, align = "center", downloadButton(
            "downloadDotPlot", "dot plot.pdf",
            style = 'padding:5px; font-size:80%')),

          fluidRow(tags$br()),
          fluidRow(tags$br()),
          fluidRow(tags$br()),
          column(12, uiOutput("plot.uiDatFeatPlotV4"), align = "center"),
          fluidRow(tags$br()),
          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(8, tags$br()),
            column(8, tags$b("Mismatches or genes not present"),
                "(if applicable)", tags$b(":")),
            column(8, uiOutput("notInDot")),
            column(8, tags$hr()),

            fluidRow(tags$br()),
            column(12, uiOutput("plot.uiDotPlotF"))
          )
        )
      )
    ),


    # ================ #
    tabPanel("Heatmap", #fluid = FALSE,
      sidebarLayout(fluid = TRUE,
        
        sidebarPanel(fluid = FALSE, width = 4,
          column(12, align = "left  ",
            textInput("PhmapGenes",
              "Insert gene name or ensembl ID:",
              value = "mpeg1.1 mfap4 lcp1 f13a1b tnfa lyz ctss1 txn ccl34b.1"),
            checkboxInput("pHmapClust",
              label = "Check box to enable row clustering.", value = FALSE)),

          column(12, align = "center",
            column(12, align = "center", uiOutput("cellSelectHmap")),
              column(12, tags$br()),
            actionButton("runPhmap", "Generate Plots",
              style = 'padding:5px; font-size:80%')),
          column(12, tags$hr(width = "50%"), align = "center"),

          column(12, align = "center", downloadButton(
            "downloadPhmap", "Download Heatmap",
            style = 'padding:5px; font-size:80%')),

          # column(12, tags$br()),
          # column(12, align = "center", "Plot download (png):"),
          # column(12, tags$br()),
          # column(12, align = "center", downloadButton(
          #   "downloadPhmap", "heatmap.png",
          #   style = "padding:5px; font-size:80%")),

          fluidRow(tags$br()),
          fluidRow(tags$br()),
          column(12, align = "left", tags$b('Note:'),'Highly expressed genes have a
          tendency to "wash out" the color values of other genes on this 
          heatmap. It might be useful to remove the higher expressed genes
          to get a better visualization of genes with less extreme values.'),
          fluidRow(tags$br()),
          fluidRow(tags$br()),
          column(12, uiOutput("plot.uiDatFeatPlotV6"), align = "center"),
          fluidRow(tags$br()),
          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(8, tags$br()),
            column(8, tags$b("Mismatches or genes not present"),
                "(if applicable)", tags$b(":")),
            column(8, uiOutput("notInPhmap")),
            column(8, tags$hr()),

            fluidRow(tags$br()),
            column(12, uiOutput("plot.uiPheatmapF"))
          )
        )
      )
    ),


    # ================ #
    tabPanel("Differential Expression", fluid = TRUE,
      sidebarLayout(
        
        sidebarPanel(
          uiOutput("idents"),

          column(12, align = "center",
            uiOutput("diffOut1"),
            fluidRow(tags$br()),
            uiOutput("diffOut2")),
            
          column(12, tags$hr(width = "50%"), align = "center"),

          fluidRow(tags$br()),
          column(12, align = "center", numericInput("pValCutoff",
            "Input adjusted p-value cutoff:", value = 0.05, min = 0.00,
            step = 0.001, max = 1.00, width = "50%")),

          fluidRow(tags$br()),
          column(12, align = "center",
            actionButton("runDiffExp", "Run Differential Expression",
            style = "padding:5px; font-size:80%")),

          column(12, tags$hr(width = "50%"), align = "center"),

          column(12, align = "center",
            downloadButton("downloadDiffExp",
            "Download Results",
            style = 'padding:5px; font-size:80%')),

          fluidRow(tags$br()),
          fluidRow(tags$br()),
          fluidRow(tags$br()),
          column(12, uiOutput("plot.uiDatFeatPlotV5"), align = "center"),
          fluidRow(tags$br()),
          fluidRow(tags$br())
        ),    
        mainPanel(
          fluidRow(
            column(8, tags$br()),
            # column(8, tags$b(uiOutput("SelectedDataDiff "))), 
            column(12, align = "left", class = "diffExpMain",
              uiOutput("diffTable"),
                tags$b('Click "Run Differential Expression"')
            )
          )
        )
      )
    )
  ), style = 'width:1580px;')#,
# shinyDebuggingPanel::withDebuggingPanel()
)


# =============================================================================
shinyApp(ui = ui, server = server)


# ==== Command line tools
# Upload app to shinyaqpps.io (must be directory that contains app.R)
# start R session

# Deploy to shinyapps.io
# rsconnect::deployApp('/Volumes/projects/ddiaz/Analysis/Scripts/rsconnect/shinyapps.io/all_macrophage_datasets_Nicolas', account = 'piotrowskilab')

# Execute app locally
# options(shiny.reactlog=TRUE, shiny.fullstacktrace = TRUE); shiny::runApp('/Volumes/projects/ddiaz/Analysis/Scripts/rsconnect/shinyapps.io/all_macrophage_datasets_Nicolas/app.R')
# options(shiny.reactlog=TRUE, shiny.fullstacktrace = TRUE); shiny::runApp('/Users/ddiaz/Desktop/macrophage-old')

# Logs
# rsconnect::showLogs(account = 'piotrowskilab', appName = 'all_macrophage_datasets_Nicolas')
