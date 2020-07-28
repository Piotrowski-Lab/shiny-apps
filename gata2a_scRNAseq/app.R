library(shiny)
library(ggplot2)
library(Seurat)
library(shinythemes)
library(shinycssloaders)
library(shinyBS)

# Grep multiple items
multiGrep2 <- function(toMatch, toSearch, ...) {
  toMatch <- ifelse(grepl("*", toMatch),
    gsub("\\*","\\\\*", toMatch), toMatch <- toMatch)
  
  toMatch <- paste(toMatch, collapse = "|")
  inCommon <- grep(toMatch, toSearch, value = FALSE)
  return(inCommon)
}

# ======== Reading in data
# A complete Seurat object with cluster info
seurat_obj <- readRDS(
  "./data/SeurObj_cell_update_gata2a_pos_Seurat3_v1.0_.RDS")

# For converting ensembl IDs to comman genes names
gene_df <- read.table("./data/Danio_Features_unique_Ens98_v1.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE)

gene_df <- gene_df[gene_df$Gene.name.uniq %in% rownames(seurat_obj),]
ens_id <- gene_df$Gene.stable.ID
com_name <- gene_df$Gene.name.uniq

smpl_genes_sm <- paste0("atoh1a her4.1 dld sox4a*1 foxp4 crip1")


# ================================== Server ===================================
server <- function(input, output) {


  # ======== Gene Database ======== #
  GeneDB <- function() {
    selected <- unlist(strsplit(input$MyText, " "))
    
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


  # ======== Feature Plot ======== #
  FeaturePlotF <- function(){
    seurat_obj <- seurat_obj
    genes_select <- unlist(strsplit(input$MyText, " "))
    
    ifelse(genes_select %in% gene_df$gene_name,
      genes_select <- as.character(genes_select),
        
      ifelse(genes_select %in% gene_df$gene_id,
        genes_select <- as.character(
          gene_df[gene_df$gene_id %in% genes_select, 2]),
            "gene not in database")
    )

    feat <- FeaturePlot(seurat_obj, genes_select,
      reduction = "umap", cols = c("azure3", "blue3"),
      combine = FALSE, pt.size = input$CellSize)

    for(k in 1:length(feat)) {
      feat[[k]] <- feat[[k]] + labs(x = "UMAP 1", y = "UMAP 2") + 
      theme(axis.text.x = element_blank(), legend.position="none",
      axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.y = element_blank(), axis.title = element_text(size = 18),
      panel.border = element_rect(colour = "#FFFFFF", fill = NA, size = 1))
    }
  return(cowplot::plot_grid(plotlist = feat, ncol = 1))
  }

  output$myFeaturePlotF <- renderPlot({
    input$runPlots
    isolate({
      p <- FeaturePlotF()
      print(p)
    })
  })

  getHeightFeat <- function(){
    l <- length(unlist(strsplit(input$MyText, " ")))
    h <- as.character(ceiling(l) * 1100)
    h <- paste0(h, "px")
    return(h)
  }

  output$plot.uiFeaturePlotF <- renderUI({
    input$runPlots
    isolate({
      h <- getHeightFeat()
      plotOutput("myFeaturePlotF",
        width = "1100px",
        height = h)
    })
  })

  output$downloadFeaturePlotF <- downloadHandler(
    filename = "UMAP_plot.pdf",
    content = function(file){
      pdf(file,
        height = as.numeric(input$MyHeight),
        width = as.numeric(input$MyWidth),
        onefile = FALSE)
      FeaturePlotF()
      dev.off()
    }
  )


  # ======== Heatmap ======== #
  HmapF <- function() {
    seurat_obj <- seurat_obj 
    genes_select <- unlist(strsplit(input$MyText, " "))
    
    ifelse(genes_select %in% gene_df$gene_name,
      genes_select <- as.character(genes_select),
        
      ifelse(genes_select %in% gene_df$gene_id,
        genes_select <- as.character(
          gene_df[gene_df$gene_id %in% genes_select, 2]),
            "gene not in database")
    )

    DoHeatmap(seurat_obj, genes_select)
  }

  output$myHmapF <- renderPlot({
    input$runPlots
    isolate({p <- HmapF()
    print(p)})
  })

  getHeightHmap <- function(){
    l <- length(unlist(strsplit(input$MyText, " ")))
    h <- paste0(as.character(l * 35), "px")
    return(h)
  }

  output$plot.uiHmapF <- renderUI({
    input$runPlots
    isolate({
      h <- getHeightHmap()
      plotOutput("myHmapF",
        width = paste0(input$SlideWidthHmap, "px"),
        height = h)
    })
  })


  # ======== Violin Plot ======== #
  VlnPlotF <- function(){
    seurat_obj <- seurat_obj
    genes_select <- unlist(strsplit(input$MyText, " "))
    
    ifelse(genes_select %in% gene_df$gene_name,
      genes_select <- as.character(genes_select),
        
      ifelse(genes_select %in% gene_df$gene_id,
        genes_select <- as.character(
          gene_df[gene_df$gene_id %in% genes_select, 2]),
            "gene not in database")
    )

    # VlnPlot(
    #   seurat_obj,
    #   genes_select,
    #   point.size.use = input$CellSize,
    #   nCol = 2) #, use.raw=T

    g <- VlnPlot(seurat_obj, genes_select,
    pt.size = input$CellSize, combine = FALSE)

    for(k in 1:length(g)) {
      g[[k]] <- g[[k]] + theme(legend.position = "none")
    }
    return(cowplot::plot_grid(plotlist = g, ncol = 1))
  }

  getHeightVln <- function(){
    l <- length(unlist(strsplit(input$MyText, " ")))
    h <- as.character(ceiling(l) * 500)
    h <- paste0(h, "px")
    return(h)
  }

  output$plot.uiVlnPlotF <- renderUI({
    input$runPlots
    isolate({
      h <- getHeightVln()
      plotOutput("myVlnPlotF",
        width = "1200px",
        height = h)
    })
  })

  output$myVlnPlotF <- renderPlot({
        input$runPlots
    isolate({
      p <- VlnPlotF()
      print(p)})
  })

  output$downloadVlnPlot <- downloadHandler(
    filename = "violin_plot.pdf",
    content = function(file){
      pdf(file,
        height = as.numeric(input$MyHeight),
        width = as.numeric(input$MyWidth),
        onefile = FALSE) 
      p <- VlnPlotF()
      print(p)
      dev.off()
      }
    )

  getHmapSizeRatio <- function(){
    l <- length(unlist(strsplit(input$MyText, " ")))
    r <- as.numeric(input$SlideWidthHmap / (l * 35))
    return(r)
  }

  output$downloadHmapF <- downloadHandler(
    filename = "heatmap.pdf",
    content = function(file){
      pdf(file,
        height = as.numeric(input$MyHeight),
        width = as.numeric(input$MyWidth),
        onefile = FALSE)
      p <- HmapF()
      print(p)
      dev.off()
    }
  )

  # ======== Download meta data ======== #
  output$downloadClusterMarkers <- downloadHandler(
    filename = "cluster_markers.xlsx",
    content = function(file){
      file.copy("./data/cluster_markers.xlsx", file)
    }
  )

  output$downloadDatasetGenes <- downloadHandler(
    filename = "genes_in_dataset.xlsx",
    content = function(file){
      file.copy("./data/genes_in_dataset.xlsx", file)
    }
  )

  output$downloadUMAPclusters <- downloadHandler(
    filename = "UMAP_clusters.pdf",
    content = function(file){
      file.copy("./data/UMAP_clusters.pdf", file)
    }
  )
#shinyDebuggingPanel::makeDebuggingPanelOutput()
}

# ================================== UI ======================================


ui <- fixedPage(theme = shinytheme("paper"),
  tags$head(includeCSS("./www/styles.css")),
  headerPanel("gata2a scRNA-seq"),
  
  sidebarLayout(
    sidebarPanel(
      
      # box to paste ensembl ids
      textInput("MyText", 
        tags$b("Insert genes for all plots:"),
        value = smpl_genes_sm),

      sliderInput("CellSize",
        "Select point size on UMAP and violin plots:",
        min = 0.25, max = 2, 
        value = 0.5, step = 0.25,
        ticks = FALSE, width = "95%"),

      column(12, align = "center",
        actionButton("runPlots", "Generate Plots",
          style = 'padding:8px; font-size:80%')),

      column(12, tags$br()),
      column(12, 'Type at least 1 gene name or Ensembl ID
        separated by a space. You can enter a list of common names,
        or a list of Ensembl IDs, but not both. Please refer to the',
        tags$b('instructions'),'on the welcome page for details.'),
      
      column(12, tags$hr()),
      tags$b("PDF Download Dimensions"),

      textInput("MyHeight",
        "Select plot height (inch) for plot download:", value = "10"),
      
      textInput("MyWidth",
        "Select plot width (inch) for plot download:", value = "10"),

      column(12,"Please note that the download 
        dimensions are independent of the figure dimensions
        displayed on the app itself.",
        align = "left"),

      column(12, tags$br()),
      fluidRow(
        column(12, align = "center",
          downloadButton(
            "downloadFeaturePlotF", "UMAP.pdf",
            style = "padding:8px; font-size:80%"),
          downloadButton(
            "downloadHmapF", "Heatmap.pdf",
            style = "padding:8px; font-size:80%"),
          downloadButton(
            "downloadVlnPlot", "Violin.pdf",
            style = "padding:8px; font-size:80%")
        )
      ),
      column(12, tags$br()),
      fluidRow(
        column(12, id = "umapSidebar", align = "center",
          tags$br(),
          tags$img(src='UMAP_clusters.png',
            width = "100%", height = "100%")
          ),
        column(12, tags$br()),
        column(12, align = "left",
          "Right click and select",
            '"open image in new tab"', "for a larger view")
      )
    ),
  
    mainPanel(
      tabsetPanel(
        tabPanel("Welcome!",
        column(12, id = "welcomePanel",
          fixedRow(
            column(12, tags$br()),
            column(10, align = "center", 
              downloadButton(
                "downloadUMAPclusters", "UMAP clusters.pdf",
                style = "padding:8px; font-size:80%")),
            
            column(12, tags$br()),
            column(10, align = "center",
              tags$img(
                src='UMAP_clusters.png',
                width = "75%",
                height = "75%")
            ),
            column(10, tags$hr()),

            column(12, tags$b("Instructions")),
            column(10,
              'All genes available for plotting can be downloaded in the
              Excel spreadsheet below labeled "genes in dataset", using either
              Ensembl IDs or common names from the',
              tags$b("Gene.name.unique"),'column as input. You cannot, however,
              mix common names with Ensembl IDs in the same query. Groups of
              genes can be directly copied/pasted from the spreadsheet into 
              the app input field and will have the necessary spacing 
              by default. Please note that this data set was produced with the
              Ensembl 98 gene annotation in zebrafish (genome version 10).
              We therefore recommend using Ensembl gene IDs as input, 
              as common gene names can change with annotation updates.
              All genes used for this dataset, as well genes enriched in
              each cluster, can be downloaded below:'),

            column(12, tags$br()),
            column(10, align = "center",
              column(5, align = "center", offset = 1,
                downloadButton("downloadClusterMarkers",
                  "Cluster Markers", style = "padding:8px; font-size:80%")
              ),
              column(5, align = "center",
                downloadButton("downloadDatasetGenes",
                  "Genes in Data Set", style = "padding:8px; font-size:80%")
              )
            ),

            column(12, tags$br()),
            column(10,
              'Additionally, there may be individual common names that have
              multiple Ensemble gene IDs. In most cases, these repeated 
              common names will have a period appended with a number to
              distinguish each repeat expressed in the data set
              (e.g. sox4a, sox4a.1). You can verify if a gene name is repeated
              by comparing the "Associated.Gene.Name" and "Gene.name.unique"
              columns in both the "cluster marker" and "genes in data set"
              spreadsheets.',
              offset = 0),

            column(10, tags$br()),
            column(12, tags$b("Gene Name DB")),
            column(10,
              'The gene name DB (database) tab is a way to quickly determine 
              whether or not a list of genes with common names can be plotted.
              This database contains all 32,105 zebrafish coding and non-coding
              features in the Ensembl 98 annotation. When you enter at least
              two genes/features with one correct name (atoh1a works great as 
              a default value) the table will populate with
              any common name that matches a portion of the text,
              and whether or not those features can be plotted. If the
              table returns "TRUE" in the third column that gene
              can be plotted (FALSE indicates otherwise).
              For example, if you enter "atoh1a wnt notch" in the search
              field, the table will dislpay atoh1a and all features with wnt and
              notch in their common name, as well as their corresponding
              Ensembl IDs. Only genes that are expressed in at least 3
              cells can be plotted (i.e. those that return TRUE in the 3rd column).',
              offset = 0),
            column(10, tags$hr())
          )
        )
        ),
        tabPanel("Gene Name DB",
          fluidRow(
            column(12, id = "geneDB", tags$br(),
              uiOutput("GeneDB"))
          )
        ),
        tabPanel("Feature Plot",
          fixedRow(
            column(12, id = "featPlot", tags$br(),
              uiOutput("plot.uiFeaturePlotF"))
          )
        ),
        tabPanel("Heatmap",
          fluidRow(
            column(12, tags$br()),
            column(4, sliderInput("SlideWidthHmap",
              "Change plot width (pixels):", ticks = FALSE,
              min = 0, max = 7000, value = 2000)),
            column(12, id = "hmapPlot",
              uiOutput("plot.uiHmapF"))
          )
        ),
        tabPanel("Violin Plot",
          fluidRow(
            column(12, id = "vlnPlot", tags$br(),
              uiOutput("plot.uiVlnPlotF"))
            )
        )
      )
    )
  )#,
#shinyDebuggingPanel::withDebuggingPanel()
)


#__________
#run the shinyApp
shinyApp(ui = ui, server = server)

# bash command to run locally
# options(shiny.reactlog=TRUE, shiny.fullstacktrace = TRUE); shiny::runApp("/Volumes/projects/ddiaz/Analysis/Scripts/rsconnect/shinyapps.io/gata2a_scRNAseq/app.R")

# Deploy to shinyapps.io
# rsconnect::deployApp('/Volumes/projects/ddiaz/Analysis/Scripts/rsconnect/shinyapps.io/gata2a_scRNAseq', account = 'piotrowskilab')
