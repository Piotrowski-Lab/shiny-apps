library(BiocManager)
options(repos = BiocManager::repositories())
#options(rsconnect.max.bundle.size=3145728000)
	library(shiny)
	library(cowplot)
	library(Seurat)
	library(monocle3)
	library(ggplot2)
	library(shinythemes)
	library(shinyWidgets)
	library(dplyr)
	library(rsconnect)
	library(reshape2)
	library(gridExtra)
	library(gtable)
	library(grid)
	library(devtools)
library(sctransform)

# devtools::load_all("/n/projects/nt2473/Analysis/Scripts/SeuratExtensions")
# devtools::load_all("/n/projects/nt2473/Analysis/Scripts/CellTrajectoryExtensions")
# local <- TRUE
# if(local){
# devtools::load_all("/Volumes/projects/nt2473/Analysis/Scripts/SeuratExtensions")
# devtools::load_all("/Volumes/projects/nt2473/Analysis/Scripts/CellTrajectoryExtensions")
# }
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

# =================== import SeuratExtension functios
cleanUMAP <- function(plot_obj, dark = FALSE, axis_title_size = 16,
		plot_title = "") {
	if (dark == TRUE) {
		plot_obj <- plot_obj + labs(title = element_text(plot_title), x = "UMAP 1", y = "UMAP 2") +
			theme(plot.background = element_rect(fill = "black",
						colour = "black", size = 0.5, linetype = "solid"),
					axis.text.x = element_blank(), axis.ticks.x = element_blank(),
					axis.line.x = element_blank(), axis.text.y = element_blank(),
					axis.ticks.y = element_blank(), axis.line.y = element_blank(),
					axis.title = element_text(size = axis_title_size, color = "white"),
					legend.text = element_text(color = c("white"), size = 10)) +
			ggtitle(plot_title)
	} else {
		plot_obj <- plot_obj + labs(title = element_text(plot_title), x = "UMAP 1", y = "UMAP 2") +
			theme(plot.title = element_text(plot_title), axis.text.x = element_blank(),
					axis.ticks.x = element_blank(), axis.line.x = element_blank(),
					axis.text.y = element_blank(), axis.ticks.y = element_blank(),
					axis.line.y = element_blank(),
					axis.title = element_text(size = axis_title_size))
	}
}

# ========================== import CellTrajectoryExtension functions
branch_nodes <- function(cds,reduction_method="UMAP"){
	g = principal_graph(cds)[[reduction_method]]
		branch_points <- which(igraph::degree(g) > 2)
		branch_points = branch_points[branch_points %in% root_nodes(cds, reduction_method) == FALSE]
		return(branch_points)
}

leaf_nodes <- function(cds,reduction_method="UMAP"){
	g = principal_graph(cds)[[reduction_method]]
		leaves <- which(igraph::degree(g) == 1)
		leaves = leaves[leaves %in% root_nodes(cds, reduction_method) == FALSE]
		return(leaves)
}

root_nodes <- function(cds, reduction_method="UMAP"){
	g = principal_graph(cds)[[reduction_method]]
		root_pr_nodes <- which(names(igraph::V(g)) %in%
				cds@principal_graph_aux[[reduction_method]]$root_pr_nodes)
		names(root_pr_nodes) <-
		cds@principal_graph_aux[[reduction_method]]$root_pr_nodes
		return(root_pr_nodes)
}


	files <- list.files("./data", pattern = ".RDS", full.names = TRUE)
file_list <- list()

	print("Loading Seurat objects...")
	for (i in 1:length(files)) {
		file_list[[i]] <- readRDS(files[i])
#DefaultAssay(file_list[[i]]) <- "RNA"
	}
print("done.")


# ! =========== items to check/change for project {START}
#name each element in obj list
names(file_list) <- as.character(c("ptime_main_traj", "ptime_central_traj", 
			"cds", "seurat_obj"))

#split elements in obj list into independent variables (by given names)
list2env(file_list,envir=.GlobalEnv)

#preserve dataset colors 
	cell_type_trt <- levels(cds$cell.type.and.trt)
type_trt_cols <- gg_color_hue(length(cell_type_trt))


	smpl_genes_single <- paste0("atoh1a")
	smpl_genes_sm <- paste0("atoh1a her4.1")
	smpl_genes_lg <- paste0("atoh1a her4.1 hes2.2 dld sox4a*1 myclb gadd45gb.1",
			" insm1a wnt2 sost sfrp1a pcna mki67 isl1 slc1a3a glula lfng cbln20 ebf3a",
			" znf185 si:ch211-229d2.5 si:ch73-261i21.5 spaca4l foxp4 crip1")

	app_title <- "Neuromast Regeneration Trajectory Inference"

	gene_df <- read.table("./data/Danio_Features_unique_Ens91_v2.tsv",
			sep = "\t", header = TRUE, stringsAsFactors = FALSE)

	branch <- "master" # CHECK BEFORE DEPLOYMENT!
	app_name <- "regen_pseudotime_traj_Sungmin"
# ! =========== {END}


	ens_id <- gene_df$Gene.stable.ID
	com_name <- gene_df$Gene.name.uniq


# =========== Server
	source(paste0("./app_server.R"), local = TRUE)

# # =========== UI

	source(paste0("./app_ui.R"), local = TRUE)

	print("Size of all Seurat objects:")
	print(object.size(file_list), units = "MB")


# =========== # Execute app
	shinyApp(ui = ui, server = server) # MUST be at last line of app.R
