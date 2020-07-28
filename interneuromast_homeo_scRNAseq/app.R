library(shiny)
library(cowplot)
library(Seurat)
library(ggplot2)
library(shinythemes)
library(shinyWidgets)
library(dplyr)
library(rsconnect)
library(profvis)
library(hrbrthemes)
library(reshape2)

"%||%" <- devtools:::`%||%`


# ===================== installing devel version of Seurat to help rsconnect to shin.io 

#install.packages('devtools')
library(devtools)
dev_mode(on=T)
#devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)


# when finished do:

#dev_mode(on=F)  #and you are back to having stable Seurat

# ===================== 

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
  selected <- unlist(strsplit(input, " ")) #splits the input genes into a list, sep by  " "
  
  ifelse(selected %in% com_name,
         selected <- gene_df[com_name %in% selected, 3],
         
         ifelse(selected %in% ens_id,
                selected <- gene_df[ens_id %in% selected, 3],"")
  )
  len <- length(selected)
  return(len)
}


files <- list.files("./data", pattern = "TRIMMED", full.names = TRUE)
file_list <- list()

print("Loading Seurat objects...")
for (i in 1:length(files)) {
  file_list[[i]] <- readRDS(files[i])
  DefaultAssay(file_list[[i]]) <- "RNA"
  #change early-HC to young-HCs
  file_list[[i]]@meta.data$cell.type.ident <- plyr::revalue(
    file_list[[i]]@meta.data$cell.type.ident, c("early-HCs" = "young-HCs"))
  print("checking cells per data.set")
  print(addmargins(table(Idents(file_list[[i]]),file_list[[i]]$data.set)))
}
print("done.")

#hmap_files <- list.files("./data", pattern = "mtx", full.names = TRUE)
#hmap_list <- list()

#print("Loading heatmap matrices...")
#for (i in 1:length(hmap_files)) {
#  hmap_list[[i]] <- readRDS(hmap_files[i])
#}
#print("done.")


# ! =========== items to check/change for project {START}
#file_list <- file_list[c(6,5,1:4)]
#hmap_list <- hmap_list[c(2,1,3)]

names(file_list) <- as.character(c(
  "interneuromast"))
#names(hmap_list) <- as.character(c("LOG", "CLR", "RC"))

trt_colors <- c("green3", "gold", "darkorange",
                "deeppink", "mediumorchid1", "deepskyblue", "blue")

smpl_genes_sm <- paste0("atoh1a her4.1")
smpl_genes_lg <- paste0("atoh1a her4.1 hes2.2 dld sox4a*1 myclb gadd45gb.1",
                        " insm1a wnt2 sost sfrp1a pcna mki67 isl1 slc1a3a glula lfng cbln20 ebf3a",
                        " znf185 si:ch211-229d2.5 si:ch73-261i21.5 spaca4l foxp4 crip1")
smpl_gene_single <- paste0("atoh1a")

app_title <- "Interneuromast Homeo scRNA-seq Analysis"

gene_df <- read.table("./data/Danio_Features_unique_Ens91_v2.tsv",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)

branch <- "master" # CHECK BEFORE DEPLOYMENT!
app_name <- "interneuromast_homeo_scRNAseq"
# ! =========== {END}


ens_id <- gene_df$Gene.stable.ID
com_name <- gene_df$Gene.name.uniq

# =========== Server
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
    if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      "seurat_clusters"
    } else if ("cell.type.ident" %in% colnames(seurat_obj@meta.data)) {
      "cell.type.ident"
    } else if ("data.set" %in% colnames(seurat_obj@meta.data)){
      "data.set"
    #}
    }else{
      "tree.ident"}
  }
  
  whichTreat <- function(){
    seurat_obj <- SelectDataset()
   if ("data.set" %in% colnames(seurat_obj@meta.data)){
      "data.set"
   }
  }
  
  printTreats <- reactive({
    seurat_obj <- SelectDataset()
    print(seurat_obj)
    if (whichTreat() == "data.set") {
      sort(unique(seurat_obj@meta.data$data.set))
    } else {
      NULL # single data set
    }
  })
  
  printSubClusters <- reactive({
    seurat_obj <- SelectDataset()
    print(seurat_obj)
    print("hi")
    if (whichDataset() == "seurat_clusters") {
      sort(unique(seurat_obj@meta.data$seurat_clusters))
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
  
  printDownSampleOptions <- reactive({
    percentage <- as.numeric(c(1.00,.75,.50,.25))
    percentage
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
                               label = TRUE, label.size = 0,
                               cols = cluster_clrs, group.by = "cell.type.ident")
    } else {
      umap_clusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.10,
                               label = TRUE, label.size = 0, group.by = "cell.type.ident")
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
                               label = TRUE, label.size = 0,
                               cols = cluster_clrs, group.by = "cell.type.ident")
    } else {
      umap_clusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.10,
                               label = TRUE, label.size = 0, group.by = "cell.type.ident")
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
    
    umap_subclusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.10,
                                label = TRUE, label.size = 0)
    
    umap_subclusters <- umap_subclusters + labs(x = "UMAP 1", y = "UMAP 2") + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.line.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), axis.line.y = element_blank(),
            axis.title = element_text(size = 12), legend.position = "bottom")
    
    datfeat_list <- list(umap_subclusters, umap_clusters, umap_dataset)
    plot_h <- plot_grid(plotlist = datfeat_list, ncol = 3)
    
    plot_v <- plot_grid(plotlist = datfeat_list, ncol = 1)
    plot_v <- plot_grid(plot_v) + theme(
      plot.background = element_rect(size = 2, color = "#DCDCDC"))
    
    datfeat_list <- list(plot_h, plot_v)
    datfeat_list
  }
  
  output$myDatFeatPlotH1 <- renderPlot({DatFeatPlotF()[[1]]})
  output$plot.uiDatFeatPlotH1 <- renderUI({
    plotOutput("myDatFeatPlotH1", width = "950px", height = "450px")
  })
  
  n_panels <- 1:8
  
  #generates umaps on sidebar of each tab
  lapply(n_panels, function(i) {
    output[[paste0("myDatFeatPlotV", i)]] <- 
      renderPlot({DatFeatPlotF()[[2]]})
  })
  
  #generates umaps on sidebar of each tab
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
  
  # # ======== Stacked Violin Plot ======== #
  # 
  StkdVlnPlotF <- reactive({
    seurat_obj <- SelectDataset()
    selected <- unlist(strsplit(input$vlnStkdGenes, " "))
    
    ifelse(selected %in% com_name,
           selected <- selected[selected %in% com_name],
           
           ifelse(selected %in% ens_id,
                  selected <- gene_df[ens_id %in% selected, 3],"")
    )
    
    #seurat_obj <- seurat_obj[,IDtype() %in% input$cellIdentsStkdVln]
    
    
    
    ids <- as.list(levels(seurat_obj$cell.type.ident))
    
    
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
      obj_trt_list[[i]] <- seurat_obj[,seurat_obj[["cell.type.ident"]] == ids[[i]]]
    }
    
    stacked_violin_plot <- function(goi, obj_trt_list){
      trt_plot_list <- list()[1:length(ids)]
      names(trt_plot_list) <- ids
      for (i in 1:length(ids)) {
        vln_obj <- VlnPlot(
          obj_trt_list[[i]], features = goi, pt.size = input$ptSizeStkdVln) +
          xlab("") + ylab(ids[i]) + ggtitle("") +
          theme(legend.position = "none", axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_text(size = rel(1), angle = 0),
                axis.text.y = element_text(size = rel(1)),
                plot.margin = unit(c(-1.0, 0.5, -1.0, 0.5), "cm"))
        trt_plot_list[[i]] <- vln_obj
      }
      
      trt_plot_list[[length(trt_plot_list)]]<- trt_plot_list[[length(trt_plot_list)]] +
        theme(axis.text.x=element_text(), axis.ticks.x = element_line())
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
      options(repr.plot.width = 12, repr.plot.height = 3)
      
      grid_obj[[i]] <- stacked_violin_plot(goi = selected[[i]], obj_trt_list = obj_trt_list)
      
    }
    names(grid_obj) <- selected
    
    final_grid <- cowplot::plot_grid(plotlist = grid_obj, nrow = length(grid_obj), axis = "l", align = "hv", scale = 0.9) +
      theme(plot.margin = margin(.2, .2, .2, .2, unit = "in"))
    
    return(final_grid)
    
  })
  
  output$cellSelectStkdVln <- renderUI({ # New cell type select
    pickerInput("cellIdentsVln", "Add or remove clusters:",
                choices = as.character(printIdents()), multiple = TRUE,
                selected = as.character(printIdents()), options = list(
                  `actions-box` = TRUE), width = "85%")
  })
  
  mismatchStkdVln <- function() {
    selected <- unlist(strsplit(input$vlnStkdGenes, " "))
    
    mismatch <- ifelse(!selected %in% c(com_name,ens_id),
                       selected[!selected %in% c(com_name,ens_id)],"")
    return(mismatch)
  }
  
  output$notInStkdVln <- renderText({input$runStkdVlnPlot
    isolate({mismatchStkdVln()})
  })
  
  output$SelectedDataStkdVln <- renderText({input$runStkdVlnPlot
    isolate({input$Analysis})
  })
  
  output$myStkdVlnPlotF <- renderPlot({input$runStkdVlnPlot
    isolate({withProgress({p <- StkdVlnPlotF(); print(p)},
                          message = "Rendering plot..",
                          min = 0, max = 10, value = 10)
    })
  })
  
  getHeightStkdVln <- function() {
    l <- getLenInput(input$vlnStkdGenes)
    if (l == 1) {h <- "800"
    } else {
      h <- as.numeric(ceiling(l) * 400)
      #h <- paste0(h, "px")
    }
    return(h)
  }
  
  output$plot.uiStkdVlnPlotF <- renderUI({input$runStkdVlnPlot
    isolate({h <- getHeightStkdVln(); plotOutput("myStkdVlnPlotF",
                                                 width = "800px", height = paste0(h, "px"))})
  })
  
  output$downloadStkdVlnPlot <- downloadHandler(
    filename = "StkdViolin_plot.pdf", content = function(file) {
      png(file,
          width = 800,
          height = getHeightStkdVln() ,units = "px" )#* getLenInput(input$vlnStkdGenes))
      print(StkdVlnPlotF())
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
      
      seurat_obj_sub <- seurat_obj[rownames(seurat_obj) %in% selected,] #filters the obj obly by the selected input
      dist_mat <- dist(seurat_obj_sub@assays$RNA@data)  #turns the selected genes into a matrix
      clust <- hclust(dist_mat)
      markers_clust <- clust$labels
      
      # if (input$selectGrpDot == "data.set") {
      #     caption_txt <- paste(
      #       "selected cells:", paste(input$cellIdentsDot, collapse = ", "))
      #     stringr::str_wrap(caption_txt, width = 10)
      #   } else {
      #     ""
      #   }
      
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
      
      # if (input$selectGrpDot == "data.set") {
      #     caption_txt <- paste(
      #       "selected cells:", paste(input$cellIdentsDot, collapse = ", "))
      #     stringr::str_wrap(caption_txt, width = 10)
      #   } else {
      #     ""
      #   }
      
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
  
  output$plot.uiDotPlotF <- renderUI({input$runDotPlot
    isolate({h <- getHeightDot(); plotOutput("myDotPlotF",
                                             width = paste0(input$manAdjustDotW, "px"),
                                             height = paste0(input$manAdjustDotH, "px"))})
  })
  
  dotHeight <- function() {
    l <- getLenInput(input$dotGenes)
    l <- as.numeric(l)
    return(l)
  }
  
  output$downloadDotPlot <- downloadHandler(
    filename = "dot_plot.pdf", content = function(file) {
      png(file, height = as.numeric(input$manAdjustDotH),
          width = as.numeric(input$manAdjustDotW), units = "px")
      print(DotPlotF())
      dev.off()
    }
  )
  
  # # ======== ggplot Heatmap ======== #
  pHeatmapF <- reactive({
    clustering <- input$pHmapClust  #enable row clustering
    if (clustering == TRUE){
      seurat_obj <- SelectDataset()
      selected <- unlist(strsplit(input$PhmapGenes, " "))
      
      ifelse(selected %in% com_name,
             selected <- selected[selected %in% com_name],
             
             ifelse(selected %in% ens_id,
                    selected <- gene_df[ens_id %in% selected, 3],"")
      )
      
      seurat_obj <- seurat_obj[,IDtype() %in% input$cellIdentsHmap]
      
      seurat_obj_sub <- seurat_obj[rownames(seurat_obj) %in% selected,]
      dist_mat <- dist(seurat_obj_sub@assays$RNA@data)
      clust <- hclust(dist_mat)   #reorder genes
      markers_clust <- clust$labels
      
      dotplot <- DotPlot(seurat_obj, features = markers_clust,
                         group.by = input$selectGrpHmap)
      
      g <- ggplot(dotplot$data, aes(id, features.plot,fill= avg.exp.scaled, width = 1, height = 1)) +
        geom_tile(color = "gray", size = 1) +
        scale_fill_distiller(
          palette = "RdYlBu") +
        theme_ipsum()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5,size = 13),
              axis.title.y.right = element_text(size=13))
      
      
      g <- g + labs(title = paste("Selected analysis:",
                                  as.character(input$Analysis)), subtitle = "", caption = "") +
        theme(plot.title = element_text(face = "plain", size = 14))
      
    } else {
      seurat_obj <- SelectDataset()
      selected <- unlist(strsplit(input$PhmapGenes, " "))
      
      ifelse(selected %in% com_name,
             selected <- selected[selected %in% com_name],
             
             ifelse(selected %in% ens_id,
                    selected <- gene_df[ens_id %in% selected, 3],"")
      )
      
      seurat_obj <- seurat_obj[,IDtype() %in% input$cellIdentsHmap]
      print(input$cellIdentsHmap)
      
      
      dotplot <- DotPlot(seurat_obj, features = selected,
                         group.by = input$selectGrpHmap)
      
      g <- ggplot(dotplot$data, aes(id, features.plot,fill= avg.exp.scaled, width = 1, height = 1)) +
        geom_tile(color = "gray", size = 1) +
        scale_fill_distiller(
          palette = "RdYlBu") +
        theme_ipsum()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5,size = 13),
              axis.title.y.right = element_text(size=13))
      
      g <- g + labs(title = paste("Selected analysis:",
                                  as.character(input$Analysis)), subtitle = "", caption = "") +
        theme(plot.title = element_text(face = "plain", size = 14))
      
    }
    
    return(g)
    
  })
  
  #renders the drop-down box w/ Ident choices
  output$cellSelectHmap <- renderUI({ # New cell type select
    pickerInput("cellIdentsHmap", "Add or remove clusters:",
                choices = as.character(printIdents()), multiple = TRUE,
                selected = as.character(printIdents()), options = list(
                  `actions-box` = TRUE), width = "85%")
  })
  
  mismatchPhmap <- function() {
    selected <- unlist(strsplit(input$PhmapGenes, " "))
    
    mismatch <- ifelse(!selected %in% c(com_name, ens_id),
                       selected[!selected %in% c(com_name, ens_id)],"")
    return(mismatch)
  }
  
  #prints the mismatches or genes not present (for ui.R)
  output$notInPhmap <- renderText({input$runPhmap
    isolate({mismatchPhmap()})
  })
  
  output$SelectedDataPhmap <- renderText({input$runPhmap
    isolate({input$Analysis})
  })
  
  #renders plot w/ progress bar
  output$myPhmapF <- renderPlot({input$runPhmap
    isolate({withProgress({p <- pHeatmapF(); print(p)},
                          message = "Rendering plot..", min = 0, max = 10, value = 10)
    })
  })
  
  getHeightPhmap <- reactive({
    l <- getLenInput(input$PhmapGenes)
    h <- as.numeric(l * 35)
    return(h)
  })
  
  getWidthPhmap <- function() {
    if(input$selectGrpHmap == "cell.type.ident") {
      w <- "1200"
    } else {
      w <- "800"
    }
    return(w)
  }
  
  output$plot.uiPheatmapF <- renderUI({input$runPhmap
    isolate({h <- getHeightPhmap(); plotOutput("myPhmapF",
                                               width = paste0(input$manAdjustHmapW, "px"),
                                               height = paste0(input$manAdjustHmapH, "px"))})
  })
  
  #download
  output$downloadhmap <- downloadHandler(
    filename = "heatmap.png", content = function(file) {
      png(file, height = as.numeric(input$manAdjustHmapH),
          width = as.numeric(input$manAdjustHmapW), units = "px")
      print(pHeatmapF())
      dev.off()
    }
  )
  
  # # # ======== Individual Cell ggplot Heatmap ======== #
  IndvpHeatmapF <- reactive({
    clustering <- input$IndvpHmapClust  #enable row clustering
    if (clustering == TRUE){
      seurat_obj <- SelectDataset()
      selected <- unlist(strsplit(input$IndvPhmapGenes, " "))
      
      ifelse(selected %in% com_name,
             selected <- selected[selected %in% com_name],
             
             ifelse(selected %in% ens_id,
                    selected <- gene_df[ens_id %in% selected, 3],"")
      )
      
      seurat_obj <- seurat_obj[,IDtype() %in% input$cellIdentsIndvHmap]
      
      seurat_obj <- seurat_obj[, sample(Cells(seurat_obj), size = round(as.numeric(input$cellDownSampleIndvHmap)*length(colnames(seurat_obj))), replace=F)]
      
      print(addmargins(table(seurat_obj$cell.type.ident))) #check
      
      seurat_obj_sub <- seurat_obj[rownames(seurat_obj) %in% selected,]
      dist_mat <- dist(seurat_obj_sub@assays$RNA@data)
      clust <- hclust(dist_mat)   #reorder genes
      markers_clust <- clust$labels
      
      group.by <- input$selectGrpIndvHmap #choose group.by parameter
      cells <- NULL
      col.min = -2.5
      col.max = 2.5
      
      cells <- cells %||% colnames(x = seurat_obj)
      
      data <- as.data.frame(x = t(x = as.matrix(x = GetAssayData(
        object = seurat_obj, slot = "data")[markers_clust, cells, drop = FALSE])))
      
      
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
      
      #preserve identity order
       if (group.by == "seurat_clusters"){
        data$id <- factor(data$id, levels = levels(seurat_obj$seurat_clusters))
      }else{
        data$id <- factor(data$id, levels = levels(seurat_obj$cell.type.ident))
      }
      
      g <- ggplot(data, aes(Cell, Feature,fill= value)) +
        geom_tile(height = .95, width = 2) +
        scale_fill_distiller(
          palette = "RdYlBu") +
        theme_ipsum()+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y.right = element_text(size=13),panel.spacing = unit(.25, "lines"),
              strip.text.x  = element_text(angle = 90, vjust = 0.5, hjust=.5,size = 12)) + 
        facet_grid( ~ id, space = 'free', scales = 'free')
      
    } else {
      seurat_obj <- SelectDataset()
      selected <- unlist(strsplit(input$IndvPhmapGenes, " "))
      
      ifelse(selected %in% com_name,
             selected <- selected[selected %in% com_name],
             
             ifelse(selected %in% ens_id,
                    selected <- gene_df[ens_id %in% selected, 3],"")
      )
      
      seurat_obj <- seurat_obj[,IDtype() %in% input$cellIdentsIndvHmap]
      print(input$cellIdentsIndvHmap)
      
      seurat_obj <- seurat_obj[, sample(Cells(seurat_obj), size = round(as.numeric(input$cellDownSampleIndvHmap)*length(colnames(seurat_obj))), replace=F)]
      
      print(input$cellDownSampleIndvHmap)
      print(addmargins(table(seurat_obj$cell.type.ident))) #check
      
      
      group.by <- input$selectGrpIndvHmap #choose group.by parameter
      cells <- NULL
      col.min = -2.5
      col.max = 2.5
      
      cells <- cells %||% colnames(x = seurat_obj)
      
      data <- as.data.frame(x = t(x = as.matrix(x = GetAssayData(
        object = seurat_obj, slot = "data")[selected, cells, drop = FALSE])))
      
      
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
      
      #preserve identity order
      if (group.by == "seurat_clusters"){
        data$id <- factor(data$id, levels = levels(seurat_obj$seurat_clusters))
      }else{
        data$id <- factor(data$id, levels = levels(seurat_obj$cell.type.ident))
      }
      
      g <- ggplot(data, aes(Cell, Feature,fill= value)) +
        geom_tile(height = .95, width = 2) +
        scale_fill_distiller(
          palette = "RdYlBu") +
        theme_ipsum()+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y.right = element_text(size=13),panel.spacing = unit(.25, "lines"),
              strip.text.x  = element_text(angle = 90, vjust = 0.5, hjust=.5,size = 12)) + 
        facet_grid( ~ id, space = 'free', scales = 'free')
      
      
    }
    return(g)
    
  })
  
  #renders the drop-down box w/ Ident choices
  output$cellSelectIndvHmap <- renderUI({ # New cell type select
    pickerInput("cellIdentsIndvHmap", "Add or remove clusters:",
                choices = as.character(printIdents()), multiple = TRUE,
                selected = as.character(printIdents()), options = list(
                  `actions-box` = TRUE), width = "85%")
  })
  
  #renders the drop-down box w/ downsample choices
  output$SelectDownSamplePropIndvHmap <- renderUI({ # New cell type select
    pickerInput("cellDownSampleIndvHmap", "Choose downsample proportion:",
                choices = as.character(printDownSampleOptions()), multiple = FALSE,
                selected = as.character(printDownSampleOptions()[1]), options = list(
                  `actions-box` = TRUE), width = "85%")
  })
  
  
  mismatchIndvPhmap <- function() {
    selected <- unlist(strsplit(input$IndvPhmapGenes, " "))
    
    mismatch <- ifelse(!selected %in% c(com_name, ens_id),
                       selected[!selected %in% c(com_name, ens_id)],"")
    return(mismatch)
  }
  
  #prints the mismatches or genes not present (for ui.R)
  output$notInIndvPhmap <- renderText({input$runIndvPhmap
    isolate({mismatchIndvPhmap()})
  })
  
  output$SelectedDataIndvPhmap <- renderText({input$runIndvPhmap
    isolate({input$Analysis})
  })
  
  #renders plot w/ progress bar
  output$myIndvPhmapF <- renderPlot({input$runIndvPhmap
    isolate({withProgress({p <- IndvpHeatmapF(); print(p)},
                          message = "Rendering plot..", min = 0, max = 10, value = 10)
    })
  })
  
  getHeightIndvPhmap <- reactive({
    l <- getLenInput(input$IndvPhmapGenes)
    h <- as.numeric(l * 35)
    return(h)
  })
  
  getWidthIndvPhmap <- function() {
    if(input$selectGrpIndvHmap == "cell.type.ident.by.data.set") {
      w <- "1600"
    } else {
      w <- "800"
    }
    return(w)
  }
  
  # output$plot.uiIndvpHeatmapF <- renderUI({input$runIndvPhmap
  #   isolate({
  #     w <- paste0(getWidthIndvPhmap()); h <- paste0(getHeightIndvPhmap())
  #     plotOutput("myIndvPhmapF", width = paste0(w, "px"), height = paste0(h, "px"))
  #   })
  # })
  
  output$plot.uiIndvpHeatmapF <- renderUI({input$runIndvPhmap
    isolate({h <- getHeightIndvPhmap(); plotOutput("myIndvPhmapF",
                                                   width = paste0(input$manAdjustIndvHmapW, "px"),
                                                   height = paste0(input$manAdjustIndvHmapH, "px"))})
  })
  
  #download
  output$downloadIndvhmap <- downloadHandler(
    filename = "IndvHeatmap.png", content = function(file) {
      png(file, height = as.numeric(input$manAdjustIndvHmapH),
          width = as.numeric(input$manAdjustIndvHmapW), units = "px")
      print(IndvpHeatmapF())
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
      print("hello")
      group1 <- rownames(meta[meta$seurat_clusters %in% subset1,])
      group2 <- rownames(meta[meta$seurat_clusters %in% subset2,])
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
                choices = as.character(printSubClusters()), multiple = TRUE,
                selected = as.character(printSubClusters())[1], options = list(
                  `actions-box` = TRUE), width = "80%")
  })
  
  output$diffOut2 <- renderUI({
    pickerInput("identText2", tags$b("Group 2 - negative FC"),
                choices = as.character(printSubClusters()), multiple = TRUE,
                selected = as.character(printSubClusters())[2],options = list(
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

# =========== UI
ui <- fixedPage(theme = shinythemes::shinytheme("lumen"), # paper lumen cosmo
                tags$head(includeCSS(paste0("./www/styles.css"))),
                div(headerPanel(app_title), style = 'width:1560px;'),
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
                                            tags$b("Select Analysis"),
                                            
                                            column(12, tags$br()),
                                            pickerInput("Analysis", label = "",
                                                        choices = list(Combined = names(file_list)),
                                                        selected = "all she-pos. cells", width = "50%")
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
                                                   )
                           ),
                  
                  
                  # ================ #
                  tabPanel("Gene Database", fluid = TRUE,
                           sidebarLayout(
                             sidebarPanel(
                               textInput("dbGenes", "Insert gene name or ensembl ID:",
                                         value = "gadd45gb.1 slc1a3a znf185 si:ch73-261i21.5"),
                               fluidRow(tags$br()),
                               column(12, uiOutput("plot.uiDatFeatPlotV6"), align = "center"),
                               fluidRow(tags$br())
                             ),
                             mainPanel(fluidRow(
                               column(11, tags$br()),
                               uiOutput("GeneDB")
                             )
                             )
                           )
                  ),
                  
                  
                  # ================ #
                  tabPanel("Feature Plots", fluid = FALSE,
                           sidebarLayout(fluid = TRUE,
                                         
                                         sidebarPanel(width = 4,
                                                      
                                                      column(12, align = "left",
                                                             textInput("featureGenes", "Insert gene name or ensembl ID:",
                                                                       value = smpl_genes_sm)),
                                                      
                                                      column(12, align = "center",
                                                             actionButton("runFeatPlot", "Generate Plots",
                                                                          style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$hr(width = "50%"), align = "center"),
                                                      column(12, align = "center",
                                                             downloadButton("downloadFeaturePlotF", "Download png",
                                                                            style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$br()),
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
                                                                    numericInput("ptSizeFeature", "Input cell size:", value = 0.50,
                                                                                 min = 0.25, step = 0.25, max = 2.00, width = "100%"))
                                                      ),
                                                      
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
                                                                           value = smpl_genes_sm)),
                                                      
                                                      column(12, align = "center",
                                                             actionButton("runVlnPlot", "Generate Plots",
                                                                          style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$hr(width = "50%"), align = "center"),
                                                      column(12, align = "center", downloadButton(
                                                        "downloadVlnPlot", "Download pdf",
                                                        style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$br()),
                                                      column(12, align = "center", uiOutput("cellSelectVln")), # New
                                                      
                                                      column(12, tags$br()),
                                                      column(12, align = "center",
                                                             column(6,
                                                                    radioGroupButtons("selectGrpVln",
                                                                                      "Group cells by:", choices = list(Cluster = "cell.type.ident",
                                                                                                                        Subcluster = "seurat_clusters"
                                                                                                                        ), width = "100%")),
                                                             column(6,
                                                                    numericInput("ptSizeVln", "Input cell size:", value = 0.25,
                                                                                 min = 0.00, step = 0.75, max = 2.00, width = "80%"))
                                                      ),
                                                      
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
                  
                  
                  #    # ================ # Stacked Vln Plts
                  #    tabPanel("Stacked Violin Plots", #fluid = FALSE,
                  #       sidebarLayout(fluid = TRUE,
                  # 
                  #       sidebarPanel(fluid = FALSE, width = 4,
                  #         column(12, textInput("vlnStkdGenes", width = "100%",
                  #                         "Insert gene name or ensembl ID:",
                  #                         value = smpl_genes_lg)),
                  #    
                  #         column(12, align = "center",
                  #           actionButton("runStkdVlnPlot", "Generate Plots",
                  #                        style = 'padding:5px; font-size:80%')),
                  #    
                  #         column(12, tags$hr(width = "50%"), align = "center"),
                  #         column(12, align = "center", downloadButton(
                  #      "downloadStkdVlnPlot", "Download pdf",
                  #      style = 'padding:5px; font-size:80%')),
                  #    
                  #         column(12, tags$br()),
                  #         column(12, align = "center", uiOutput("cellSelectStkdVln")), # New
                  #    
                  #         column(12, tags$br()),
                  #           column(12, align = "center",
                  #           column(6,
                  #                  radioGroupButtons("selectGrpStkdVln",
                  #                                    "Group cells by:", choices = list(Time = "data.set",
                  #                                                                      Cluster = "cell.type.ident"), width = "100%")),
                  #           column(6,
                  #                  numericInput("ptSizeStkdVln", "Input cell size:", value = 0.25,
                  #                               min = 0.00, step = 0.75, max = 2.00, width = "80%"))
                  #    ),
                  #    
                  #    fluidRow(tags$br()),
                  #    fluidRow(tags$br()),
                  #    column(12, uiOutput("plot.uiDatFeatPlotV8"), align = "center"),
                  #    fluidRow(tags$br()),
                  #    fluidRow(tags$br())
                  #   ),
                  # 
                  #   mainPanel(
                  #     fluidRow(
                  #     column(8, tags$br()),
                  #     column(8, tags$b("Gene mismatches"), "(if present)", tags$b(":")),
                  #     column(8,uiOutput("notInStkdVln")),
                  #     column(8, tags$hr()),
                  # # column(8, tags$b(uiOutput("SelectedDataVln"))),
                  #     column(12, uiOutput("plot.uiStkdVlnPlotF")
                  #     )
                  #     )
                  #     )
                  #     )
                  #     ),
                  
                  # ================ # Stacked Vln Plts
                  tabPanel("Stacked Violin Plots", #fluid = FALSE,
                           sidebarLayout(fluid = TRUE,
                                         
                                         sidebarPanel(fluid = FALSE, width = 4,
                                                      column(12, textInput("vlnStkdGenes", width = "100%",
                                                                           "Insert gene name or ensembl ID:",
                                                                           value = smpl_genes_sm)),
                                                      #column(12, tags$br()),
                                                      column(12, em("Please select no more than 3 genes at a time for computational efficiency")),
                                                      column(12, tags$br()),
                                                      
                                                      column(12, align = "center",
                                                             actionButton("runStkdVlnPlot", "Generate Plots",
                                                                          style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$hr(width = "50%"), align = "center"),
                                                      column(12, align = "center", downloadButton(
                                                        "downloadStkdVlnPlot", "Download pdf",
                                                        style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$br()),
                                                      column(12, align = "center", uiOutput("cellSelectStkdVln")), # New
                                                      
                                                      column(12, tags$br()),
                                                      column(12, align = "center",
                                                             column(6,
                                                                    radioGroupButtons("selectGrpStkdVln",
                                                                                      "Group cells by:", choices = list(Cluster = "cell.type.ident"), 
                                                                                      width = "100%")),
                                                             column(6,
                                                                    numericInput("ptSizeStkdVln", "Input cell size:", 
                                                                                 value = 0.00, min = 0.00, step = 0.75, 
                                                                                 max = 2.00, width = "80%"))
                                                      ),
                                                      
                                                      fluidRow(tags$br()),
                                                      fluidRow(tags$br()),
                                                      column(12, uiOutput("plot.uiDatFeatPlotV8"), align = "center"),
                                                      fluidRow(tags$br()),
                                                      fluidRow(tags$br())
                                         ),
                                         
                                         mainPanel(
                                           fluidRow(
                                             column(8, tags$br()),
                                             column(8, tags$b("Gene mismatches"), "(if present)", tags$b(":")),
                                             column(8,uiOutput("notInStkdVln")),
                                             column(8, tags$hr()),
                                             # column(8, tags$b(uiOutput("SelectedDataVln"))),
                                             column(12, uiOutput("plot.uiStkdVlnPlotF")
                                             )
                                           )
                                         )
                           )
                  ),
                  
                  # # ================ #
                  # tabPanel("Ridge Plots", #fluid = FALSE,
                  #   sidebarLayout(fluid = TRUE,
                  #     
                  #     sidebarPanel(fluid = FALSE, width = 4,
                  #       column(12, textInput("rdgGenes", width = "100%",
                  #         "Insert gene name or ensembl ID:",
                  #         value = smpl_genes_sm)),
                  #       
                  #       column(12, align = "center",
                  #         actionButton("runRdgPlot", "Generate Plots",
                  #         style = 'padding:5px; font-size:80%')),
                  # 
                  #       column(12, tags$hr(width = "50%"), align = "center"),
                  #       column(12, align = "center", downloadButton(
                  #         "downloadRdgPlot", "Download pdf",
                  #         style = 'padding:5px; font-size:80%')),
                  # 
                  #       column(12, tags$br()),
                  #       column(12, align = "center", uiOutput("cellSelectRdg")), # New
                  #       
                  #       column(12, tags$br()),
                  #        column(12, align = "center",
                  #         column(6,
                  #           radioGroupButtons("selectGrpRdg",
                  #             "Group cells by:", choices = list(Time = "data.set",
                  #               Cluster = "cell.type.ident"), width = "100%")),
                  #         column(6,
                  #           numericInput("ptSizeRdg", "Input cell size:", value = 0.25,
                  #           min = 0.00, step = 0.75, max = 2.00, width = "80%"))
                  #       ),          
                  #       
                  #       fluidRow(tags$br()),
                  #       fluidRow(tags$br()),
                  #       column(12, uiOutput("plot.uiDatFeatPlotV3"), align = "center"),
                  #       fluidRow(tags$br()),
                  #       fluidRow(tags$br())
                  #     ),
                  #     
                  #     mainPanel(
                  #       fluidRow(
                  #         column(8, tags$br()),
                  #         column(8, tags$b("Gene mismatches"), "(if present)", tags$b(":")),
                  #         column(8,uiOutput("notInRdg")),
                  #         column(8, tags$hr()),
                  #         # column(8, tags$b(uiOutput("SelectedDataRdg"))),
                  #         column(12, uiOutput("plot.uiRdgPlotF")
                  #         )
                  #       )
                  #     )
                  #   )
                  # ),
                  # 
                  
                  # ================ #
                  tabPanel("Dot Plot", #fluid = FALSE,
                           sidebarLayout(fluid = TRUE,
                                         
                                         sidebarPanel(fluid = FALSE, width = 4,
                                                      column(12, align = "left  ",
                                                             textInput("dotGenes",
                                                                       "Insert gene name or ensembl ID:",
                                                                       value = smpl_genes_lg),
                                                             checkboxInput("dPlotClust",
                                                                           label = "Check box to enable row clustering.", value = FALSE)),
                                                      
                                                      column(12, align = "center",
                                                             actionButton("runDotPlot", "Generate Plots",
                                                                          style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$hr(width = "50%"), align = "center"),
                                                      column(12, align = "center", downloadButton(
                                                        "downloadDotPlot", "Download pdf",
                                                        style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$br()),
                                                      column(12, align = "center", uiOutput("cellSelectDot")), # New
                                                      
                                                      column(12, tags$br()),
                                                      column(12, align = "center",
                                                             column(6,
                                                                    radioGroupButtons("selectGrpDot",
                                                                                      "Group cells by:", choices = list(Cluster = "cell.type.ident",
                                                                                                                        Subcluster = "seurat_clusters"), width = "100%")),
                                                             column(6,
                                                                    numericInput("dotScale", "Dot diameter:", value = 10, min = 4,
                                                                                 step = 1, max = 20, width = "80%"), align = "center")
                                                      ),
                                                      
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
                                             
                                             column(8, align = "left",
                                                    # column(4,  align = "center", "Manual figure adjustment:",
                                                    #   column(11, style = "padding-top: 8px;",
                                                    #     switchInput("manAdjustDot", value = FALSE))),
                                                    column(3, align = "left", numericInput(
                                                      "manAdjustDotW", label = "Width (pixels):", value = 600, step = 50,
                                                      width = "100%")),
                                                    column(3,  align = "left", numericInput(
                                                      "manAdjustDotH", label = "Height (pixels):", value = 900, step = 50,
                                                      width = "100%"))
                                             ),
                                             fluidRow(tags$br()),
                                             column(12, uiOutput("plot.uiDotPlotF"))
                                           )
                                         )
                           )
                  ),
                  
                  # ================ # DoHeatmap
                  tabPanel("Heat Map", #fluid = FALSE,
                           sidebarLayout(fluid = TRUE,
                                         
                                         sidebarPanel(fluid = FALSE, width = 4,
                                                      column(12, align = "left  ",
                                                             textInput("PhmapGenes",
                                                                       "Insert gene name or ensembl ID:",
                                                                       value = smpl_genes_lg),
                                                             checkboxInput("pHmapClust",
                                                                           label = "Check box to enable row clustering.", value = FALSE)),
                                                      
                                                      column(12, align = "center",
                                                             actionButton("runPhmap", "Generate Plots",
                                                                          style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$hr(width = "50%"), align = "center"),
                                                      column(12, align = "center", downloadButton(
                                                        "downloadhmap", "Download pdf",
                                                        style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$br()),
                                                      column(12, align = "center", uiOutput("cellSelectHmap")), # New
                                                      
                                                      column(12, tags$br()),
                                                      column(12, align = "center",
                                                             column(12,
                                                                    radioGroupButtons("selectGrpHmap",
                                                                                      "Group cells by:", 
                                                                                      choices = list(Cluster = "cell.type.ident",
                                                                                                     Subcluster = "seurat_clusters"), 
                                                                                      width = "100%"))
                                                             
                                                      ),
                                                      
                                                      fluidRow(tags$br()),
                                                      fluidRow(tags$br()),
                                                      column(12, uiOutput("plot.uiDatFeatPlotV7"), align = "center"),
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
                                             column(8, align = "left",
                                                    column(3, align = "left", numericInput(
                                                      "manAdjustHmapW", label = "Width (pixels):", value = 600, step = 50,
                                                      width = "100%")),
                                                    column(3,  align = "left", numericInput(
                                                      "manAdjustHmapH", label = "Height (pixels):", value = 900, step = 50,
                                                      width = "100%"))
                                             ),
                                             fluidRow(tags$br()),
                                             column(12, uiOutput("plot.uiPheatmapF"))
                                           )
                                         )
                           )
                  ),
                  
                  #================ # ggplot Indv. Cell heatmap
                  tabPanel("Indv. Cell Heatmap", #fluid = FALSE,
                           sidebarLayout(fluid = TRUE,
                                         
                                         sidebarPanel(fluid = FALSE, width = 4,
                                                      column(12, align = "left  ",
                                                             textInput("IndvPhmapGenes",
                                                                       "Insert gene name or ensembl ID:",
                                                                       value = smpl_genes_lg),
                                                             checkboxInput("IndvpHmapClust",
                                                                           label = "Check box to enable row clustering.", value = FALSE)),
                                                      
                                                      column(12, align = "center",
                                                             actionButton("runIndvPhmap", "Generate Plots",
                                                                          style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$hr(width = "50%"), align = "center"),
                                                      column(12, align = "center", downloadButton(
                                                        "downloadIndvhmap", "Download pdf",
                                                        style = 'padding:5px; font-size:80%')),
                                                      
                                                      column(12, tags$br()),
                                                      column(12, align = "center", uiOutput("cellSelectIndvHmap")), # New
                                                      
                                                      column(12, tags$br()),
                                                      column(12, align = "center", uiOutput("SelectDownSamplePropIndvHmap")), #downsample drop down
                                                      
                                                      column(12, tags$br()),
                                                      column(12, align = "center",
                                                             column(12,
                                                                    radioGroupButtons("selectGrpIndvHmap",
                                                                                      "Group cells by:",
                                                                                      choices = list(
                                                                                      Cluster = "cell.type.ident", 
                                                                                      Subcluster = "seurat_clusters"),
                                                                                      width = "100%"))
                                                             
                                                      ),
                                                      
                                                      fluidRow(tags$br()),
                                                      fluidRow(tags$br()),
                                                      column(12, uiOutput("plot.uiDatFeatPlotV3"), align = "center"),
                                                      fluidRow(tags$br()),
                                                      fluidRow(tags$br())
                                         ),
                                         
                                         mainPanel(
                                           fluidRow(
                                             column(8, tags$br()),
                                             column(8, tags$b("Mismatches or genes not present"),
                                                    "(if applicable)", tags$b(":")),
                                             column(8, uiOutput("notInIndvPhmap")),
                                             column(8, tags$hr()),
                                             column(8, align = "left",
                                                    column(3, align = "left", numericInput(
                                                      "manAdjustIndvHmapW", label = "Width (pixels):", value = 600, step = 50,
                                                      width = "100%")),
                                                    column(3,  align = "left", numericInput(
                                                      "manAdjustIndvHmapH", label = "Height (pixels):", value = 900, step = 50,
                                                      width = "100%"))
                                             ),
                                             fluidRow(tags$br()),
                                             column(12, uiOutput("plot.uiIndvpHeatmapF"),style = "overflow-y: scroll;overflow-x: scroll;")
                                           )
                                         )
                           )
                  ),
                  
                  
                  
                  # # ================ # Phmap - omit
                  # tabPanel("Heatmap", #fluid = FALSE,
                  #   fixedRow(
                  #     column(12, tags$br()),
                  # 
                  #     column(5, align = "left",
                  #       column(12, align = "left",
                  #         column(12,
                  #           textInput("PhmapGenes", width = "100%",
                  #           "Insert gene name or ensembl ID:",
                  #             value = smpl_genes_lg),
                  #           checkboxInput("pHmapClust",
                  #             label = "Check box to enable row clustering.",
                  #             value = FALSE),
                  #           column(12, align = "center", uiOutput("cellSelectHmap")),
                  #           column(12, tags$br())
                  #           ),
                  #         
                  #         column(12, align = "center",
                  #           actionButton("runPhmap", "Generate Plots",
                  #             style = 'padding:5px; font-size:80%'),
                  #         downloadButton("downloadPhmap", "Download pdf",
                  #             style = 'padding:5px; font-size:80%')),
                  #         column(12, tags$br())
                  #       )
                  #     ),
                  # 
                  #     column(7, align = "left",
                  #       column(12, tags$b("Mismatches or genes not present"),
                  #         "(if applicable)", tags$b(":")),
                  #       column(12, uiOutput("notInPhmap")),
                  # 
                  #       column(12, tags$hr()),
                  #       column(9, align = "left", tags$b('Note:'),
                  #       'Highly expressed genes have a tendency to "wash out" the color 
                  #       values of genes with lower expression on this heatmap. It might 
                  #       be useful to remove the higher expressed genes to get a better 
                  #       visualization of genes with less extreme values. You can also 
                  #       change the expression normalization method to decrease/increase
                  #       the effect highly expressed genes. You can find details on each 
                  #       method in the ',
                  #         tags$a(href = "https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/NormalizeData",
                  #           tags$b("Seurat documentation")), "."),
                  #       
                  #       column(12, tags$br()),
                  #       column(12, align = "left",
                  #         radioGroupButtons("mtxSelectHmap", "Normalization method:",
                  #           choices = list(Log = "LOG", CLR = "CLR", RC = "RC"),
                  #           width = "100%")
                  #       )
                  #     ),
                  #     
                  #     column(12, align = "center", tags$hr(width = "100%")),
                  #     column(12, tags$b("Selected analysis: all she-pos. cells")),
                  #     column(12, tags$br()),
                  #     column(12, class = "hmapID", uiOutput("plot.uiPheatmapF"))
                  #   )
                  # ),
                  
                  
                  # ================ #
                  tabPanel("Differential Expression", fluid = TRUE,
                           sidebarLayout(
                             
                             sidebarPanel(
                               uiOutput("idents"),
                               
                               column(12, align = "center",
                                      uiOutput("diffOut1"),
                                      fluidRow(tags$br()),
                                      uiOutput("diffOut2")),
                               
                               column(12, tags$br()),
                               column(12, align = "center", uiOutput("cellSelectDiff")),
                               column(12, tags$hr(width = "50%"), align = "center"),
                               
                               fluidRow(tags$br()),
                               column(12, align = "center",
                                      pickerInput("statSelectDiff", label = "Select statistical test:",
                                                  multiple = FALSE, selected = "wilcox", width = "210px",
                                                  choices = list(wilcox = "wilcox", bimodal = "bimod", ROC = "roc",
                                                                 t = "t", negbinom = "negbinom", poisson = "poisson", LR = "LR",
                                                                 MAST = "MAST", DESeq2 = "DESeq2"))),
                               
                               fluidRow(tags$br()),
                               column(12, align = "center", numericInput("pValCutoff",
                                                                         "Input adjusted p-value cutoff:", value = 0.05, min = 0.00,
                                                                         step = 0.001, max = 1.00, width = "210px")),
                               
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
                                        column(12, align = "center",
                                               fluidRow(tags$br()),
                                               tags$b('Click "Run Differential Expression"'))
                                 )
                               )
                             )
                           )
                  )
), style = 'width:1500px;')#,
# shinyDebuggingPanel::withDebuggingPanel()
                )

print("Size of all Seurat objects:")
print(object.size(file_list), units = "MB")


# ======================================================== Deploy/execute tools 


if (FALSE) { # Not run
  # # Deploy from local
  # if(branch == "master") {
  # rsconnect::deployApp(paste0("/Volumes/easystore/SIMR_2019/shiny-apps-main/", app_name),
  #   account = "piotrowskilab")
  # }
  
  # Deploy from server
  # if(branch == "master") {
  #   rsconnect::deployApp(paste0("/home/ntran2/bgmp/shiny-apps-main/", app_name),
  #     account = "piotrowskilab")
  # }
  
  #Execute app locally - please define app name first
  # options(shiny.reactlog = TRUE, shiny.fullstacktrace = TRUE)
  # shiny::runApp(paste0("home/ntran2/bgmp/shiny-apps-main/", app_name, "/app.R"))
  
  # Logs
  rsconnect::showLogs(account = 'piotrowskilab',
                      appName = app_name)
}


# =========== # Execute app
shinyApp(ui = ui, server = server) # MUST be at last line of app.R

# profvis({
#   runApp("../../../shiny-apps-main/smrtseq_vs_10X_scRNAseq/")
# },prof_output = "./data")







