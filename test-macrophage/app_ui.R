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
                "Group cells by:", choices = list(Time = "data.set",
                  Cluster = "cell.type.ident"), width = "100%")),
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


    # ================ #
    tabPanel("Ridge Plots", #fluid = FALSE,
      sidebarLayout(fluid = TRUE,
        
        sidebarPanel(fluid = FALSE, width = 4,
          column(12, textInput("rdgGenes", width = "100%",
            "Insert gene name or ensembl ID:",
            value = smpl_genes_sm)),
          
          column(12, align = "center",
            actionButton("runRdgPlot", "Generate Plots",
            style = 'padding:5px; font-size:80%')),

          column(12, tags$hr(width = "50%"), align = "center"),
          column(12, align = "center", downloadButton(
            "downloadRdgPlot", "Download pdf",
            style = 'padding:5px; font-size:80%')),

          column(12, tags$br()),
          column(12, align = "center", uiOutput("cellSelectRdg")), # New
          
          column(12, tags$br()),
           column(12, align = "center",
            column(6,
              radioGroupButtons("selectGrpRdg",
                "Group cells by:", choices = list(Time = "data.set",
                  Cluster = "cell.type.ident"), width = "100%")),
            column(6,
              numericInput("ptSizeRdg", "Input cell size:", value = 0.25,
              min = 0.00, step = 0.75, max = 2.00, width = "80%"))
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
                "Group cells by:", choices = list(Time = "data.set",
                  Cluster = "cell.type.ident"), width = "100%")),
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

            fluidRow(tags$br()),
            column(12, uiOutput("plot.uiDotPlotF"))
          )
        )
      )
    ),


    # ================ #
    tabPanel("pBulk Heatmap", #fluid = FALSE,
      fixedRow(
        column(12, tags$br()),

        column(5, align = "left",
          column(12, align = "left",
            column(12,
              textInput("PhmapGenes", width = "100%",
              "Insert gene name or ensembl ID:",
                value = smpl_genes_lg),
              checkboxInput("pHmapClust",
                label = "Check box to enable row clustering.",
                value = FALSE),
              column(12, align = "center", uiOutput("cellSelectHmap")),
              column(12, tags$br())
              ),
            
            column(12, align = "center",
              actionButton("runPhmap", "Generate Plots",
                style = 'padding:5px; font-size:80%'),
            downloadButton("downloadPhmap", "Download pdf",
                style = 'padding:5px; font-size:80%')),
            column(12, tags$br())
          )
        ),

        column(7, align = "left",
          column(12, tags$b("Mismatches or genes not present"),
            "(if applicable)", tags$b(":")),
          column(12, uiOutput("notInPhmap")),

          column(12, tags$hr()),
          column(9, align = "left", tags$b('Note:'),
          'Highly expressed genes have a tendency to "wash out" the color 
          values of genes with lower expression on this heatmap. It might 
          be useful to remove the higher expressed genes to get a better 
          visualization of genes with less extreme values. You can also 
          change the expression normalization method to decrease/increase
          the effect highly expressed genes. You can find details on each 
          method in the ',
            tags$a(href = "https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/NormalizeData",
              tags$b("Seurat documentation")), "."),
          
          column(12, tags$br()),
          column(12, align = "left",
            radioGroupButtons("mtxSelectHmap", "Normalization method:",
              choices = list(Log = "LOG", CLR = "CLR", RC = "RC"),
              width = "100%")
          )
        ),
        
        column(12, align = "center", tags$hr(width = "100%")),
        column(12, tags$b("Selected analysis: all she-pos. cells")),
        column(12, tags$br()),
        column(12, class = "hmapID", uiOutput("plot.uiPheatmapF"))
      )
    ),


    # ================ #
    tabPanel("Cell Heatmap", #fluid = FALSE,
      sidebarLayout(fluid = TRUE,
        
        sidebarPanel(fluid = FALSE, width = 4,
          column(12, align = "left  ",
          textInput("dotGenes",
            "Insert gene name or ensembl ID:",
            value = smpl_genes_lg),
          checkboxInput("cellHmapClust",
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
                "Group cells by:", choices = list(Time = "data.set",
                  Cluster = "cell.type.ident"), width = "100%")),
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

            fluidRow(tags$br()),
            column(12, uiOutput("plot.uiDotPlotF"))
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