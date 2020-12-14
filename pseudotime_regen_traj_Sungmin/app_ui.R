# ================================================================================ ui 
ui <- fixedPage(theme = shinythemes::shinytheme("lumen"), # paper lumen cosmo
		tags$head(includeCSS(paste0("./www/styles.css"))),
		div(headerPanel(app_title), style = 'width:1660px;'),
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
# pickerInput("Analysis", label = "",
#             choices = list(Combined = names(file_list)),
#             selected = "all she-pos. cells", width = "50%")
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
									tags$a(href = "https://webfs/n/projects/ddiaz/Analysis/Scripts/sb2191-regen/regen-summary/site/",
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
					column(12, uiOutput("plot.uiDatFeatPlotV3"), align = "right"),
					fluidRow(tags$br())
					),
				mainPanel(fluidRow(
						column(6, tags$br()),
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
          
					#download svg button
					column(12, tags$hr(width = "50%"), align = "center"),
					column(12, align = "center",
						downloadButton("downloadSVGFeaturePlotF", "Download SVG",
							style = 'padding:5px; font-size:80%')),

					#download pdf button
					column(12, tags$hr(width = "50%"), align = "center"),
					column(12, align = "center",
					       downloadButton("downloadPDFFeaturePlotF", "Download PDF",
					                      style = 'padding:5px; font-size:80%')),
					
					#download png button
					column(12, tags$hr(width = "50%"), align = "center"),
					column(12, align = "center",
					       downloadButton("downloadPNGFeaturePlotF", "Download PNG",
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
								numericInput("nodeSize", "Input Node Size:",
									value = 3, min = 1, step = 0.25, max = 5, width = "100%")),
							column(6, align = "left",
								numericInput("ptSizeFeature", "Input cell size:", value = 1.0,
									min = 0.25, step = 0.25, max = 4.00, width = "100%"))
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
	tabPanel("Pseudotime Line Plots", fluid = FALSE,
			sidebarLayout(fluid = TRUE,

				sidebarPanel(width = 4,

					column(12, align = "left",
						textInput("ptimeLinePlotGenes", "Insert gene name or ensembl ID:",
							value = smpl_genes_sm)),

					column(12, align = "center",
						actionButton("runPtimeLinePlot", "Generate Plots",
							style = 'padding:5px; font-size:80%')),

					column(12, tags$hr(width = "50%"), align = "center"),
					column(12, align = "center",
						downloadButton("downloadSVGPtimeLinePlotF", "Download SVG",
							style = 'padding:5px; font-size:80%')),
					
					column(12, tags$hr(width = "50%"), align = "center"),
					column(12, align = "center",
					       downloadButton("downloadPDFPtimeLinePlotF", "Download PDF",
					                      style = 'padding:5px; font-size:80%')),
					
					column(12, tags$hr(width = "50%"), align = "center"),
					column(12, align = "center",
					       downloadButton("downloadPNGPtimeLinePlotF", "Download PNG",
					                      style = 'padding:5px; font-size:80%')),

# column(12, tags$br()),
# column(12, align = "center", uiOutput("cellSelectFeat")),

					column(12, tags$br()),
					column(12,  align = "center", radioGroupButtons("selectGrpPtimeLinePlot",
							"Graph Choices:", choices = list(WithLegend = "WithLegend", NoLegend = "NoLegend"), width = "100%")),

					fluidRow(tags$br()),
					fluidRow(tags$br()),
					column(12, uiOutput("plot.uiDatFeatPlotV2"), align = "center"),
					fluidRow(tags$br()),
					fluidRow(tags$br())
	),

	mainPanel(
			fluidRow(
				column(8, tags$br()),
				column(8, tags$b("Mismatches or genes not present"),
					"(if applicable)", tags$b(":")),
				column(8,uiOutput("notInPtimeLinePlot")),
				column(8, tags$hr()),

				fluidRow(tags$br()),
				column(12, uiOutput("plot.uiPtimeLinePlotF")
				      )
				)
		 )
	)
 	), #end PtimeLinePlot tab
# 
# ================ #Multi-Genes Pseudotime Line Plots tab
tabPanel("Multi-Genes Pseudotime Line Plots", fluid = FALSE,
         sidebarLayout(fluid = TRUE,

                       sidebarPanel(width = 4,

                                    column(12, align = "left",
                                           textInput("multptimeLinePlotGenes", "Insert gene name or ensembl ID:",
                                                     value = smpl_genes_lg)),

                                    column(12, align = "center",
                                           actionButton("runMultiPtimeLinePlot", "Generate Plots",
                                                        style = 'padding:5px; font-size:80%')),

                                    column(12, tags$hr(width = "50%"), align = "center"),
                                    column(12, align = "center",
                                           downloadButton("downloadSVGMultiGPtimeLinePlotF", "Download SVG",
                                                          style = 'padding:5px; font-size:80%')),

                                    column(12, tags$hr(width = "50%"), align = "center"),
                                    column(12, align = "center",
                                           downloadButton("downloadPDFMultiGPtimeLinePlotF", "Download PDF",
                                                          style = 'padding:5px; font-size:80%')),

                                    column(12, tags$hr(width = "50%"), align = "center"),
                                    column(12, align = "center",
                                           downloadButton("downloadPNGMultiGPtimeLinePlotF", "Download PNG",
                                                          style = 'padding:5px; font-size:80%')),

                                    # column(12, tags$br()),
                                    # column(12, align = "center", uiOutput("cellSelectFeat")),

                                    # column(12, tags$br()),
                                    # column(12,  align = "center", radioGroupButtons("selectGrpMultiPtimeLinePlot",
                                    #                                                 "Graph Choices:", choices = list(WithLegend = "WithLegend", NoLegend = "NoLegend"), width = "100%")),

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
                           column(8,uiOutput("notInMultiPtimeLinePlot")),
                           column(8, tags$hr()),

                           fluidRow(tags$br()),
                           column(12, uiOutput("plot.uiMultiGPtimeLinePlotF")
                           )
                         )
                       )
         )
), #end multiple gene tab

# ================ # Heatmap
tabPanel("Heatmap", fluid = FALSE,
         sidebarLayout(fluid = TRUE,
                       
                       sidebarPanel(width = 4,
                                    
                                    column(12, align = "left",
                                           textInput("hmapGenes", "Insert gene name or ensembl ID:",
                                                     value = smpl_genes_lg)),
                                    
                                    column(12, align = "center",
                                           actionButton("runHmap", "Generate Plots",
                                                        style = 'padding:5px; font-size:80%')),
                                    
                                    column(12, tags$hr(width = "50%"), align = "center"),
                                    column(12, align = "center",
                                           downloadButton("downloadSVGHeatmapF", "Download SVG",
                                                          style = 'padding:5px; font-size:80%')),
                                    
                                    column(12, tags$hr(width = "50%"), align = "center"),
                                    column(12, align = "center",
                                           downloadButton("downloadPDFHeatmapF", "Download PDF",
                                                          style = 'padding:5px; font-size:80%')),
                                    
                                    column(12, tags$hr(width = "50%"), align = "center"),
                                    column(12, align = "center",
                                           downloadButton("downloadPNGHeatmapF", "Download PNG",
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
                           column(8, tags$b("Mismatches or genes not present"),
                                  "(if applicable)", tags$b(":")),
                           column(8,uiOutput("notInHmap")),
                           column(8, tags$hr()),
                           
                           column(8, align = "left",
                                  column(3, align = "left", numericInput(
                                    "manAdjustHmapW", label = "Width (inches):", value = 12, step = 1,
                                    width = "100%")),
                                  column(3,  align = "left", numericInput(
                                    "manAdjustHmapH", label = "Height (inches):", value = 9, step = 1,
                                    width = "100%"))
                           ),
                           fluidRow(tags$br()),
                           column(12, uiOutput("plot.uiHeatmapF")
                           )
                         )
                       )
         )
) #end hmap tab
	) #tabsetPanel
#end div
	,style = 'width:1550px;')) #end

