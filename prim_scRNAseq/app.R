library(shiny)
library(Seurat)
library(shinythemes)
library(shinycssloaders)
library(shinyBS)

# Grep multiple items
multiGrep <- function(toMatch, toSearch, ...){
  toMatch <- paste(toMatch, collapse = "|")
  inCommon <- grep(toMatch, toSearch, ...)
  return(inCommon)
}

# ======== Reading in data
# A complete Seurat object with cluster info
seurat_unt <- readRDS("./data/seurat_obj_MOLNG2352_v1.3_.RDS")

# For converting ensembl IDs to comman genes names
gene_names_df <- read.delim("./data/Danio_Features_unique_Ens87.tsv",
  header = TRUE, stringsAsFactors = FALSE, sep = "\t")

gene_names_df$in_dataset <- ""
gene_names_df$in_dataset <- rownames(seurat_unt@data) %in% gene_names_df$gene_name
gene_names_df <- gene_names_df[,c(1,2,5,3,4)]



# ================================== Server ===================================
server <- function(input, output) {

  # ======== Gene Name DB ======== #
  GeneDB <- function(){
    dataset_seurat <- seurat_unt 
    genes_to_select <- unlist(strsplit(input$MyText, " "))
    
    ifelse(genes_to_select %in% gene_names_df$gene_name,
      ind <- multiGrep(genes_to_select, gene_names_df$gene_name),
        ifelse(genes_to_select %in% gene_names_df$gene_id,
          ind <- multiGrep(genes_to_select, gene_names_df$gene_id),
            "gene not in dataset")
    )
    gene_names_df[ind,]
  }
  
  output$GeneDB <- renderTable({
    input$runPlots
    isolate({
      GeneDB()
    })
  })


  # ======== Violin Plot ======== #
  VlnPlotF <- function(){
    dataset_seurat <- seurat_unt
    genes_to_select <- unlist(strsplit(input$MyText, " "))
    
    ifelse(genes_to_select %in% gene_names_df$gene_name,
      genes_select <- as.character(genes_to_select),
        ifelse(genes_to_select %in% gene_names_df$gene_id,
          genes_select <- as.character(
            gene_names_df[gene_names_df$gene_id %in% genes_to_select, 2]),
              "gene not in database")
    )

    VlnPlot(
      dataset_seurat,
      genes_select,
      point.size.use = input$CellSize,
      nCol = 2) #, use.raw=T 
  }

  getHeightVln <- function(){
    l <- length(unlist(strsplit(input$MyText, " ")))
    h <- as.character(ceiling(l/2) * 500)
    h <- paste0(h, "px")
    return(h)
  }

  output$plot.uiVlnPlotF <- renderUI({
    input$runPlots
    isolate({
      h <- getHeightVln()
      plotOutput("myVlnPlotF",
        width = "1000px",
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


  # ======== Feature Plot ======== #
  FeaturePlotF <- function(){
    dataset_seurat <- seurat_unt
    genes_to_select <- unlist(strsplit(input$MyText, " "))
    
    ifelse(genes_to_select %in% gene_names_df$gene_name,
      genes_select <- as.character(genes_to_select),
        ifelse(genes_to_select %in% gene_names_df$gene_id,
          genes_select <- as.character(
            gene_names_df[gene_names_df$gene_id %in% genes_to_select, 2]),
              "gene not in database")
    )

    FeaturePlot(dataset_seurat, genes_select,
      nCol = 2, 
      cols.use = c("azure3", "blue3"),
      pt.size = input$CellSize,
      no.axes = TRUE)
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
    h <- as.character(ceiling(l/2) * 500)
    h <- paste0(h, "px")
    return(h)
  }

  output$plot.uiFeaturePlotF <- renderUI({
    input$runPlots
    isolate({
      h <- getHeightFeat()
      plotOutput("myFeaturePlotF",
        width = "1000px",
        height = h)
    })
  })

  output$downloadFeaturePlotF <- downloadHandler(
    filename = "tSNE_plot.pdf",
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
  HmapF <- function(){
    dataset_seurat <- seurat_unt 
    genes_to_select <- unlist(strsplit(input$MyText, " "))
    
    ifelse(genes_to_select %in% gene_names_df$gene_name,
      genes_select <- as.character(genes_to_select),
        ifelse(genes_to_select %in% gene_names_df$gene_id,
          genes_select <- as.character(
            gene_names_df[gene_names_df$gene_id %in% genes_to_select, 2]),
              "gene not in database")
    )

    DoHeatmap(
      object = dataset_seurat,
      genes.use = genes_select,
      slim.col.label = TRUE,
      group.label.rot = TRUE,
      group.label.loc = "top",
      col.low = "darkblue",
      col.mid = "lightblue",
      col.high = "red",
      cex.row = 12,
      cex.col = 12,
      remove.key = FALSE)
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

  output$downloadtSNEclusters <- downloadHandler(
    filename = "tsne_clusters.pdf",
    content = function(file){
      file.copy("./data/tsne_clusters.pdf", file)
    }
  )
#shinyDebuggingPanel::makeDebuggingPanelOutput()
}

# ================================== UI ======================================


ui <- fixedPage(theme = shinytheme("paper"),
  tags$head(includeCSS("./www/styles.css")),
  headerPanel("Zebrafish Prim scRNA-seq"),
  
  sidebarLayout(
    sidebarPanel(
      
      # box to paste ensembl ids
      textInput("MyText", 
        tags$b("Insert genes for all plots:"),
        value = "atoh1a pcna dld tekt3 slc1a3a otofb"),

      sliderInput("CellSize",
        "Select point size on t-SNE and violin plots:",
        min = 0.25, max = 4, 
        value = 1.5, step = 0.25,
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
            "downloadFeaturePlotF", "t-SNE.pdf",
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
        column(12, align = "center",
          tags$img(src='tsne_clusters.png',
            width = "100%", height = "100%"))
      )
    ),
  
    mainPanel(
      tabsetPanel(
        tabPanel("Welcome!",
          fluidRow(

            column(12, tags$br()),
            column(10, align = "center", 
              downloadButton(
                "downloadtSNEclusters", "tSNE clusters.pdf",
                style = "padding:8px; font-size:80%")),
            
            column(12, tags$br()),
            column(10, align = "center",
              tags$img(
                src='tsne_clusters.png',
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
              Ensembl 84 gene annotation in zebrafish (genome version 10).
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
              features in the Ensembl 84 annotation. When you enter at least
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
        ),
        tabPanel("Gene Name DB",
          fluidRow(
            uiOutput("GeneDB")
          )
        ),
        tabPanel("t-SNE Plot",
          fluidRow(
            column(12, tags$br()),
            uiOutput("plot.uiFeaturePlotF")
          )
        ),
        tabPanel("Heatmap",
          fluidRow(
            column(12, tags$br()),
            column(4, sliderInput("SlideWidthHmap",
              "Change plot width (pixels):", ticks = FALSE,
              min = 0, max = 7000, value = 1000)),
            column(12, tags$br()),
            uiOutput("plot.uiHmapF")
          )
        ),
        tabPanel("Violin Plot",
          fluidRow(
            column(12, tags$br()),
            uiOutput("plot.uiVlnPlotF")
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
# R -e "shiny::runApp('/Volumes/projects/ddiaz/Analysis/Scripts/rsconnect/shinyapps.io/prim_scRNAseq')"

# start R session
# rsconnect::deployApp('/Volumes/projects/ddiaz/Analysis/Scripts/rsconnect/shinyapps.io/prim_scRNAseq', account = 'piotrowskilab')
