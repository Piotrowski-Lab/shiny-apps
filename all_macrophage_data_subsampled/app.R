library(shiny)
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
# change for dataset
files <- list.files("./data", pattern = "TRIMMED", full.names = TRUE)
files <- files[c(1,5,2:4)]
seurat_obj <- readRDS(files[1])

# Load times
start <- proc.time()
file_list <- list()[1:5] # change for dataset
for (i in 1:length(files)) {
    file_list[[i]] <- readRDS(files[i])
}
elapsed <- proc.time() - start
print(elapsed)

# change for dataset
names(file_list) <- c("Homeostasis", "15 minutes post neo", 
  "1 hour post neo", "3 hours post neo", "5 hours post neo")

# Cell and gene stats
if (FALSE) {
  trt <- levels(file_list[[13]]@meta.data$data.set)
  for (i in 1:length(file_list)) {
    print(paste(names(file_list[i]), "total cells:", sep = " "))
    print(ncol(file_list[[i]]))
    m <- file_list[[i]]@meta.data
    print("Distribution:")
    print(summary(m[, "nFeature_RNA"]))
    cat("\n")
    if ("data.set" %in% colnames(m)) {
      for (j in 1:length(trt)) {
        print(names(file_list[i]))
        print(trt[j])
        print(paste(length(m$data.set[m$data.set == trt[j]]), "cells"))
        print("median gene count:")
        print(median(m[m$data.set == trt[j], "nFeature_RNA"]))
        cat("\n")
      }
    }
  }
}

gene_df <- read.table("./data/Danio_Features_unique_Ens91_v2.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gene_df <- gene_df[gene_df$Gene.name.uniq %in% rownames(seurat_obj),]

ens_id <- gene_df$Gene.stable.ID
com_name <- gene_df$Gene.name.uniq

gene_df$in.dataset <- ""
gene_df$in.dataset <- (gene_df$Gene.name.uniq %in% rownames(seurat_obj))
gene_df <- gene_df[,c(1:3,6,4:5)]
gene_df$ZFIN.ID <- paste0("=HYPERLINK(", '"', gene_df$ZFIN.ID, '"',")")

clusterColors <- gg_color_hue(
  length(levels(seurat_obj@active.ident)))


# ================================== Server ===================================


server <- function(input, output) {

  # ======== Dataset selection ======== #
  SelectDataset <- reactive({
    # browser()
    seurat_obj <- file_list[[input$DataSet]]
    print(names(file_list[input$DataSet]))
    return(seurat_obj)
  })

  whichDataset <- function() {
  seurat_obj <- SelectDataset()
  if ("data.set" %in% colnames(seurat_obj@meta.data)) {
    "data.set"} else {"tree.ident"}
  }

  printIdents <- reactive({
    seurat_obj <- SelectDataset()
    print(seurat_obj)
      if (whichDataset() == "data.set") {
        sort(unique(seurat_obj@meta.data$data.set))
      } else {
        sort(unique(seurat_obj@meta.data$tree.ident))
    }
  })

  # # ======== Gene Database ======== #
  # GeneDB <- function() {
  #   seurat_obj <- SelectDataset()
  #   selected <- unlist(strsplit(input$dbGenes, " "))

  #   ifelse(selected %in% gene_df$Gene.name.uniq,
  #     ind <- multiGrep2(selected, gene_df$Gene.name.uniq),

  #     ifelse(selected %in% gene_df$Gene.stable.ID,
  #       ind <- multiGrep2(selected, gene_df$Gene.stable.ID),
  #       "gene not in database")
  #   )
  #   gene_df[ind,]
  # }
  # output$GeneDB <- renderTable({GeneDB()})


  # ======== UMAP Cluster plot ======== #
  ClusterPlotF <- function() {
    seurat_obj <- SelectDataset()

    umap_clusters <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.25,
    label = TRUE, label.size = 3)

    umap_clusters <- umap_clusters + labs(x = "UMAP 1", y = "UMAP 2") + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.title = element_text(size = 10))
    umap_clusters
  }
  
  output$myClusterPlotF <- renderPlot({ClusterPlotF()})
  output$plot.uiClusterPlotF <- renderUI({plotOutput("myClusterPlotF",
    width = "600px", height = "500px")})


  # ======== Feature Plot ======== #
  FeaturePlotF <- function() {
    seurat_obj <- SelectDataset()
    selected <- unlist(strsplit(input$featureGenes, " "))
    
    ifelse(selected %in% com_name,
      selected <- selected[selected %in% com_name],
    
      ifelse(selected %in% ens_id,
        selected <- gene_df[ens_id %in% selected, 3],"")
    )
  
    feat <- FeaturePlot(seurat_obj, selected, reduction = "umap",
      cols = c(input$CellBackCol, input$CellForeCol), combine = FALSE,
      pt.size = input$ptSizeFeature)

    for(k in 1:length(feat)) {
      feat[[k]] <- feat[[k]] + labs(x = "UMAP 1", y = "UMAP 2") + 
      theme(axis.text.x = element_blank(), legend.position="none",
      axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.y = element_blank(), axis.title = element_text(size = 14),
      panel.border = element_rect(colour = "black", fill = NA, size = 1))
    }
  return(cowplot::plot_grid(plotlist = feat, ncol = 2))
  }

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
    if (l == 1) {h <- "750px"
    } else {
    h <- as.character(ceiling(l) * 400)
    h <- paste0(h, "px")
  }
    return(h)
  }

  output$plot.uiFeaturePlotF <- renderUI({input$runFeatPlot
    isolate({h <- getHeightFeat()
      plotOutput("myFeaturePlotF", width = "1500px", height = h)
    })
  })

  getFeatSizeRatio <- function() {
    l <- getLenInput(input$featureGenes)
    r <- as.numeric(1000 / (ceiling(l) * 750))
    return(r)
  }

  output$downloadFeaturePlotF <- downloadHandler(
    filename = "Feature_plot.pdf", content = function(file) {
      pdf(file, onefile = FALSE,
        width = getFeatSizeRatio() * as.numeric(input$featScaleFactor),
        height = as.numeric(input$featScaleFactor))
      FeaturePlotF()
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
  
    g <- VlnPlot(seurat_obj, selected,
      pt.size = input$ptSizeVln, combine = FALSE, group.by = whichDataset(),
      cols = clusterColors)

    for(k in 1:length(g)) {
      g[[k]] <- g[[k]] + theme(legend.position = "none")
    }
    return(cowplot::plot_grid(plotlist = g, ncol = 2))
  }

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
    h <- as.character(ceiling(l/2) * 500)
    h <- paste0(h, "px")
    return(h)
  }

  output$plot.uiVlnPlotF <- renderUI({input$runVlnPlot
    isolate({h <- getHeightVln(); plotOutput("myVlnPlotF",
      width = "1500px", height = h)})
  })

  getVlnSizeRatio <- function() {
    l <- getLenInput(input$vlnGenes)
    r <- as.numeric(1250 / (ceiling(l/2) * 500))
    return(r)
  }
  
  output$downloadVlnPlot <- downloadHandler(
    filename = "Violin_plot.pdf", content = function(file) {
      pdf(file, onefile = FALSE,
        width = getVlnSizeRatio() * as.numeric(input$VlnScaleFactor),
        height = as.numeric(input$VlnScaleFactor))
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
  
    g <- RidgePlot(seurat_obj, selected,
      combine = FALSE,
      group.by = whichDataset(),
      cols = clusterColors)

    for(k in 1:length(g)) {
      g[[k]] <- g[[k]] + theme(legend.position = "none")
    }

    return(cowplot::plot_grid(plotlist = g, ncol = 2))
  }

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
    h <- as.character(ceiling(l/2) * 500)
    h <- paste0(h, "px")
    return(h)
  }

  output$plot.uiRdgPlotF <- renderUI({input$runRdgPlot
    isolate({h <- getHeightRdg(); plotOutput("myRdgPlotF",
      width = "1500px", height = h)})
  })

  getRdgSizeRatio <- function() {
    l <- getLenInput(input$rdgGenes)
    r <- as.numeric(1250 / (ceiling(l/2) * 500))
    return(r)
  }
  
  output$downloadRdgPlot <- downloadHandler(
    filename = "Ridge_plot.pdf", content = function(file) {
      pdf(file, onefile = FALSE,
        width = getRdgSizeRatio() * as.numeric(input$RdgScaleFactor),
        height = as.numeric(input$RdgScaleFactor))
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

      seurat_obj_sub <- seurat_obj[rownames(seurat_obj) %in% selected,]
      dist_mat <- dist(seurat_obj_sub@assays$RNA@data)
      clust <- hclust(dist_mat)
      markers_clust <- clust$labels
      
      g <- DotPlot(seurat_obj, features = markers_clust,
        cols = "RdYlBu", dot.scale = input$dotScale,
        group.by = whichDataset())
      g <- g + coord_flip()
      g <- g + ggtitle(as.character(input$DataSet))

    } else { 
      seurat_obj <- SelectDataset()
      selected <- unlist(strsplit(input$dotGenes, " "))
      
      ifelse(selected %in% com_name,
        selected <- selected[selected %in% com_name],
      
        ifelse(selected %in% ens_id,
          selected <- gene_df[ens_id %in% selected, 3],"")
        )
      g <- DotPlot(seurat_obj, features = selected,
        cols = "RdYlBu", dot.scale = input$dotScale,
        group.by = whichDataset(),)
      g <- g + coord_flip()
      g <- g + ggtitle(as.character(input$DataSet))
    }
    return(g)
  }

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

  output$plot.uiDotPlotF <- renderUI({input$runDotPlot
    isolate({h <- getHeightDot(); plotOutput("myDotPlotF",
      width = "800px", height = h)
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
      group1 <- rownames(meta[meta$tree.ident %in% subset1,])
      group2 <- rownames(meta[meta$tree.ident %in% subset2,])
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
# test ==========================

  output$diffOut1 <- renderUI({
    pickerInput("identText1", tags$b("Group 1 - positive FC"),
      choices = as.character(printIdents()), multiple = TRUE,
      selected = as.character(printIdents())[1], options = list(
       `actions-box` = TRUE))
  })

  output$diffOut2 <- renderUI({
    pickerInput("identText2", tags$b("Group 2 - negative FC"),
      choices = as.character(printIdents()), multiple = TRUE,
      selected = as.character(printIdents())[2],options = list(
        `actions-box` = TRUE))
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
    filename = "genes_in_dataset.xlsx",
    content = function(file) {
      file.copy("./data/genes_in_dataset.xlsx", file)
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
  # Debug
  observeEvent(input$browser,{
    browser()
  })
} # Server close


# =============================== User Interface ==============================


ui <- fixedPage(theme = shinytheme("lumen"),
  tags$head(includeCSS(paste0("./www/styles.css"))),
  headerPanel("Macrophage regeneration scRNA-seq"), # change for Dataset
  tabsetPanel(
    
    # ================ #
    tabPanel("Welcome!", fluid = TRUE,
      mainPanel(

        fluidRow(tags$br()),
        fluidRow(tags$br()),

        # change for Dataset
        column(12, align = "center",
          tags$b("Select Data Set"),
          pickerInput("DataSet", label = "",
            choices = list(
              updated = names(file_list)[1:5],
            selected = "Homeostasis", width = "60%")
          )
        ),

        column(12, tags$br()),
        column(12, align = "center",
          uiOutput("plot.uiClusterPlotF")),

        column(12, tags$hr()),
        column(12, tags$b("Instructions")),
            column(12,
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
              as common gene names can change with annotation updates.
              All genes used for this dataset, as well genes enriched in
              each cluster, can be downloaded below:'),

        fluidRow(tags$br()),
        fluidRow(tags$br()),
        column(12, tags$b("An important note on ambiguous gene names")),
        column(12,
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
              Ensembl 91 zebrafish annotation.'),
        fluidRow(tags$br()),
        fluidRow(tags$br()),

        # Debug - type $('#browser').show(); in browser java console
        actionButton("browser", "browser"),
        tags$script("$('#browser').hide();")
      )
    ),


    # ================ #
    # tabPanel("Gene Database", fluid = TRUE,
    #   sidebarLayout(
    #     sidebarPanel(
    #       textInput("dbGenes", "Insert gene name or ensembl ID:",
    #         value = "mpeg1.1 mfap4"), # change for data set
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
                     
          textInput("featureGenes",
            "Insert gene name or ensembl ID:",
            value = "mpeg1.1 mfap4"), # change for data set
          
          column(12, align = "center",
            actionButton("runFeatPlot", "Generate Plots",
            style = 'padding:5px; font-size:80%')),
          
          column(12, tags$hr(width = "50%"), align = "center"),

          column(12, numericInput("ptSizeFeature",
            "Input cell size:", value = 1.00, min = 0.25,
            step = 0.25, max = 2.00, width = "50%"),
          align = "center"),

          column(12, tags$br()),
            "Choose feature plot cell color.",
          textInput("CellBackCol", "Background:", value = "azure3"),
          textInput("CellForeCol", "Foreground:", value = "blue3"),
          
          column(12, tags$br()),
          textInput("featScaleFactor",
            "Set scale factor for plot download:", value = "40"),

          column(12, align = "center",
            downloadButton("downloadFeaturePlotF", "Feature plot.pdf",
            style = 'padding:5px; font-size:80%')),

          column(12, tags$br(), align = "left",
            "For UMAP plot pdf downloads, set the scale factor",
            "to the number of genes queried multiplied by 10.",
            "This will produce a reasonably sized pdf",
            "proportional to what's displayed on the app.",
            "Increasing or decreasing the scale factor will",
            "blow-up or compress the image within the pdf.",
            "This changes how 'crowded' the cells will appear",
            "when the pdf is rendered."),

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
          column(8, tags$b(uiOutput("SelectedDataFeat"))),
          uiOutput("plot.uiFeaturePlotF")
          )
        )
      )
    ),
    
    
    # ================ #
    tabPanel("Violin Plots", #fluid = FALSE,
      sidebarLayout(fluid = TRUE,
        
        sidebarPanel(fluid = FALSE, width = 4,
          textInput("vlnGenes",
            "Insert gene name or ensembl ID:",
            value = "mpeg1.1 mfap4"), # change for data set
          
          column(12, align = "center",
            actionButton("runVlnPlot", "Generate Plots",
            style = 'padding:5px; font-size:80%')),

          column(12, tags$hr(width = "50%"), align = "center"),
          
          column(12, numericInput("ptSizeVln",
              "Input cell size:", value = 0.25, min = 0.00,
              step = 0.75, max = 2.00, width = "50%"),
          align = "center"),

          column(12, tags$br()),
          textInput("VlnScaleFactor",
            "Set scale factor for plot download:", value = "8"),
          
          column(12, align = "center", downloadButton(
            "downloadVlnPlot", "Violin plot.pdf",
            style = 'padding:5px; font-size:80%')),
            fluidRow(tags$br()),

          column(12, tags$br(), align = "left",
            "For violin plot pdf downloads, set the scale factor",
            "to the number of genes queried multiplied by 4.",
            "This will produce a reasonably sized pdf",
            "proportional to what's displayed on the app.",
            "Increasing or decreasing the scale factor will",
            "blow-up or compress the image within the pdf."),

          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(8, tags$br()),
            column(8, tags$b("Gene mismatches"), "(if present)", tags$b(":")),
            column(8,uiOutput("notInVln")),
            column(8, tags$hr()),
            column(8, tags$b(uiOutput("SelectedDataVln"))),
            uiOutput("plot.uiVlnPlotF")
          )
        )
      )
    ),


    # ================ #
    tabPanel("Ridge Plots", #fluid = FALSE,
      sidebarLayout(fluid = TRUE,
        
        sidebarPanel(fluid = FALSE, width = 4,
          textInput("rdgGenes",
            "Insert gene name or ensembl ID:",
            value = "mpeg1.1 mfap4"), # change for data set
          
          column(12, align = "center",
            actionButton("runRdgPlot", "Generate Plots",
            style = 'padding:5px; font-size:80%')),

          column(12, tags$hr(width = "50%"), align = "center"),
          
          column(12, align = "center",
            numericInput("ptSizeRdg",
              "Input cell size:", value = 0.25, min = 0.00,
              step = 0.75, max = 2.00, width = "50%")
          ),

          column(12, tags$br()),
          textInput("RdgScaleFactor",
            "Set scale factor for plot download:", value = "8"),
          
          column(12, align = "center", downloadButton(
            "downloadRdgPlot", "Ridge plot.pdf",
            style = 'padding:5px; font-size:80%')),
            fluidRow(tags$br()),

          column(12, tags$br(), align = "left",
            "For Ridge plot pdf downloads, set the scale factor",
            "to the number of genes queried multiplied by 8.",
            "This will produce a reasonably sized pdf",
            "proportional to what's displayed on the app.",
            "Increasing or decreasing the scale factor will",
            "blow-up or compress the image within the pdf."),

          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(8, tags$br()),
            column(8, tags$b("Gene mismatches"), "(if present)", tags$b(":")),
            column(8,uiOutput("notInRdg")),
            column(8, tags$hr()),
            column(8, tags$b(uiOutput("SelectedDataRdg"))),
            uiOutput("plot.uiRdgPlotF")
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
            "Insert gene name or ensembl ID:", # change for dataset
            value = "mpeg1.1 mfap4 lcp1 f13a1b tnfa lyz ctss1 txn ccl34b.1"),
          
          checkboxInput("dPlotClust",
            label = "Check box to enable row clustering.", value = FALSE)),

          column(12, align = "center",
            actionButton("runDotPlot", "Generate Plots",
            style = 'padding:5px; font-size:80%')),
          column(12, tags$hr(width = "50%"), align = "center"),

          fluidRow(tags$br()),
          column(12, numericInput("dotScale",
              "Dot diameter scaling:", value = 8, min = 4,
              step = 1, max = 20, width = "50%"),
          align = "center"),

          fluidRow(tags$br()),
          fluidRow(tags$br()),
          column(12, align = "center", downloadButton(
            "downloadDotPlot", "dot plot.pdf",
            style = 'padding:5px; font-size:80%')),
          fluidRow(tags$br()),
          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(8, tags$br()),
            column(8, tags$b("Gene mismatches"), "(if present)", tags$b(":")),
            column(8, uiOutput("notInDot")),
            column(8, tags$hr()),
            column(8, tags$b(uiOutput("SelectedDataDot"))),
            uiOutput("plot.uiDotPlotF")
          )
        )
      )
    ),


    # ================ #
    tabPanel("Differential Expression", fluid = TRUE,
      sidebarLayout(
        
        sidebarPanel(
          uiOutput("idents"),

          uiOutput("diffOut1"),
          fluidRow(tags$br()),
          uiOutput("diffOut2"),

          column(12, tags$hr(width = "50%"), align = "center"),

          fluidRow(tags$br()),
          column(12, align = "left", numericInput("pValCutoff",
            "Input adjusted p-value cutoff:", value = 0.05, min = 0.00,
            step = 0.001, max = 1.00, width = "100%")),

          fluidRow(tags$br()),
          column(12, align = "center",
            actionButton("runDiffExp", "Run Differential Expression",
            style = "padding:5px; font-size:80%")),

          column(12, tags$hr(width = "50%"), align = "center"),

          column(12, align = "center",
            downloadButton("downloadDiffExp",
            "Download Results",
            style = 'padding:5px; font-size:80%')),

          fluidRow(tags$br())
        ),
        mainPanel(
          fluidRow(
            column(8, tags$br()),
            column(8, tags$b(uiOutput("SelectedDataDiff "))),
            uiOutput("diffTable"),
            column(8, tags$br()),
            column(8, align = "center",
              tags$b('Click "Run Differential Expression"'))
          )
        )
      )
    )
  )#,
# shinyDebuggingPanel::withDebuggingPanel()
)


# =============================================================================
shinyApp(ui = ui, server = server) #

# change for dataset (app folder name)

# ==== Command line tools
# Upload app to shinyaqpps.io (must be directory that contains app.R)
# start R session

# Deploy to shinyapps.io
# rsconnect::deployApp('/Volumes/projects/ddiaz/Analysis/Scripts/rsconnect/shinyapps.io/all_macrophage_data_subsampled/', account = 'piotrowskilab')

# Execute app locally
# options(shiny.reactlog=TRUE); shiny::runApp('/Volumes/projects/ddiaz/Analysis/Scripts/rsconnect/shinyapps.io/all_macrophage_data_subsampled/app.R')

# Logs
# rsconnect::showLogs(account = 'piotrowskilab', appName = 'all_macrophage_data_subsampled')
