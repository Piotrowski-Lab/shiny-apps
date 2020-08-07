library(shiny)
library(cowplot)
library(devtools)
dev_mode(on=T)
#devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
library(ggplot2)
library(shinythemes)
library(shinyWidgets)
library(dplyr)
library(rsconnect)
library(hrbrthemes)
library(reshape2)

	"%||%" <- devtools:::`%||%`
	
	`%notin%` <- Negate(`%in%`)

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

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
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

files <- list.files("./data", pattern = "TRIMMED", full.names = TRUE)
file_list <- list()

	print("Loading Seurat objects...")
	for (i in 1:length(files)) {
		file_list[[i]] <- readRDS(files[i])
			DefaultAssay(file_list[[i]]) <- "RNA"
	}
print("done.")


# ! =========== items to check/change for project {START}
file_list <- file_list[c(1,7,2:6,8:10)]

names(file_list) <- as.character(c("neuromast-cells","neuromast-and-skin", "AP-cells", "central-cells",
                                   "DV-cells", "HC-prog", "mantle-cells", "HC-lineage-homeo-and-regen",
                                   "HC-lineage-regen", "HC-lineage-homeo"))

multiple_idents_seurObj <- as.character(c("neuromast-and-skin","neuromast-cells","HC-lineage-homeo-and-regen","HC-lineage-regen", 
                                          "HC-lineage-homeo"))

ptime_analysis <- as.character(c("HC-lineage-homeo-and-regen","HC-lineage-regen", 
                                 "HC-lineage-homeo"))

trt_colors <- c("green3", "gold", "darkorange", "red", "magenta",
		"mediumpurple1", "lightseagreen", "deepskyblue", "blue")


smpl_genes_sm <- paste0("atoh1a her4.1")
smpl_genes_lg <- paste0("atoh1a her4.1 hes2.2 dld sox4a*1 myclb gadd45gb.1",
		" insm1a wnt2 sost sfrp1a pcna mki67 isl1 slc1a3a glula lfng cbln20 ebf3a",
		" znf185 si:ch211-229d2.5 si:ch73-261i21.5 spaca4l foxp4 crip1")

app_title <- "Neuromast Regeneration scRNA-seq"

gene_df <- read.table("./data/Danio_Features_unique_Ens91_v2.tsv",
		sep = "\t", header = TRUE, stringsAsFactors = FALSE)

branch <- "master" # CHECK BEFORE DEPLOYMENT!
app_name <- "regeneration_additional_timepoints"
# ! =========== {END}


ens_id <- gene_df$Gene.stable.ID
com_name <- gene_df$Gene.name.uniq


# =========== Server
source(paste0("https://raw.githubusercontent.com/ntran95/shiny-apps/",
			branch, "/", app_name, "/app_server.R"), local = TRUE)

#source(paste0("./app_server.R"), local = TRUE)

# =========== UI
source(paste0("https://raw.githubusercontent.com/ntran95/shiny-apps/",
			branch, "/", app_name, "/app_ui.R"), local = TRUE)

#source(paste0("./app_ui.R"), local = TRUE)

print("Size of all Seurat objects:")
print(object.size(file_list), units = "MB")


# ======================================================== Deploy/execute tools 


if (FALSE) { # Not run
# Deploy from local
	if(branch == "master") {
		rsconnect::deployApp(paste0("/Volumes/projects/ddiaz/Analysis/",
					"Scripts/rsconnect/shinyapps.io/", app_name),
				account = "piotrowskilab")
	}

# Deploy from server
	if(branch == "master") {
		rsconnect::deployApp(paste0("/n/projects/ddiaz/Analysis/",
					"Scripts/rsconnect/shinyapps.io/", app_name),
				account = "piotrowskilab")
	}

#Execute app locally - please define app name first
	options(shiny.reactlog = TRUE, shiny.fullstacktrace = TRUE)
		shiny::runApp(paste0("/Volumes/projects/ddiaz/Analysis/",
					"Scripts/rsconnect/shinyapps.io/", app_name, "/app.R"))

#Execute app from desktop
		options(shiny.reactlog = TRUE, shiny.fullstacktrace = TRUE)
		shiny::runApp(paste0("~/Desktop/", app_name, "/app.R"))

# Logs
		rsconnect::showLogs(account = 'piotrowskilab',
				appName = app_name)
}


# =========== # Execute app
shinyApp(ui = ui, server = server) # MUST be at last line of app.R
