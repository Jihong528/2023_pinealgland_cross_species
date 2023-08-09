#!/home/zhengjh/miniconda3/envs/r403/bin/Rscript
# parameter

library(optparse)
library(getopt)

option_list <- list(
	make_option(c("-o", "--output"), type = "character", default = FALSE,
							action = "store", help = "This is the output directory."
	),
	make_option(c("--seurat_object_path"), type = "character", default = FALSE,
							action = "store", help = "This is absoulte path of seurat_object_path"
	),
	make_option(c("--pseudodata_path"), type = "character", default = FALSE,
							action = "store", help = "This is absoulte path of pseudodata_path"
	),
	make_option(c("--gpcr_list_path"), type = "character", default = FALSE,
							action = "store", help = "This is absoulte path of gpcr_list_path"
	),
	make_option(c("--org_name"), type = "character", default = "Gene.Symbol.(Zebrafish)",
							action = "store", help = "This is org_name used to select the GPCR you want to used"
	),
	make_option(c("--lrdb_path"), type = "character", default = "FALSE",
							action = "store", help = "This is path of lrdb."
	),
	make_option(c("--filelabel"), type = "character", default = FALSE,
							action = "store", help = "This is label used to add to the name of saved files."
	),
	make_option(c( "--min_pct"), type = "double", default = 0.1,
							action = "store", help = "This is min_pct used in FindAllmarker function[default 0.1]."
	),
	make_option(c("--logfc_threshold"), type = "double", default = 0.25,
							action = "store", help = "This is logfc_threshold used in FindAllmarker function[default 0.25]."
	),
	make_option(c("--p_val_adj"), type = "double", default = 0.05,
							action = "p_val_adj", help = "This is p_val_adj used in FindAllmarker function[default 0.05]."
	)
)


# -help
opt = parse_args(OptionParser(option_list = option_list, usage = "select_cell_gpcr [-o] [--seurat_object_path] [--pseudodata_path] [--gpcr_list_path] [--org_name]......"))
print(opt)

timestart <- Sys.time() 

setwd(opt$output)
cat("Load used R package\n")
library(openxlsx)
library(corrplot)
library(Seurat)
library(SeuratObject)
library(pheatmap)
library(RColorBrewer)
library(reshape)
library(dplyr)
library(ggplotify)
library(scLink)
library(cowplot)
library(ggplot2)
library(igraph)
### used function
IdentifyCellTypeTopGPCR <- function(seurat_object, gpcr_list, min_pct, logfc_threshold, p_val_adj){
	
	cat("Num of Used GPCR list:", unique(length(gpcr_list$genename)), "\n")
	print(head(gpcr_list))
	###  use top marker to narrow down gpcr list
	gpcr_markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = min_pct, 
																 logfc.threshold = logfc_threshold, features = gpcr_list$genename)
	gpcr_markers <- gpcr_markers[gpcr_markers$p_val_adj < p_val_adj, ]
	gpcr_markers$Family <- gpcr_list$Family[match(gpcr_markers$gene, gpcr_list$genename)]
	gpcr_markers$Class <- gpcr_list$Class[match(gpcr_markers$gene, gpcr_list$genename)]
	gpcr_markers$Subclass <- gpcr_list$Subclass[match(gpcr_markers$gene, gpcr_list$genename)]
	cat("Num of Used GPCR list after filter min_pct_cutoff:",min_pct, "logfc.threshold:", logfc_threshold, "p_val_adj:", p_val_adj, "\n")
	cat(unique(length(gpcr_markers$gene)), "\n")
	return(gpcr_markers)
}

GenerateTopGPCRScaleData <- function(seurat_object, gpcr_markers){
	
	#p1 <- DoHeatmap(seurat_object, features = gpcr_markers$gene)
	#gene_order <- rev(levels(p1$data$Feature))
	
	unique_gpcr_markers <- dplyr::select(gpcr_markers, c("gene", "Family"))
	unique_gpcr_markers <- unique_gpcr_markers[!duplicated(unique_gpcr_markers$gene), ]
	#unique_gpcr_markers <- unique_gpcr_markers[gene_order, ]
	
	## gcpr scale data
	gpcr_data <- as.matrix(seurat_object@assays$RNA@scale.data[unique_gpcr_markers$gene, ])
	
  ## annotation_col information
	cat("Generate the annotation_col used in heatmap and sort the cellnames based on celltype_level\n")
	celltype_levels <- levels(seurat_object@active.ident)
	annotation_col = data.frame( celltype = seurat_object@active.ident)
	rownames(annotation_col) <- colnames(seurat_object)
	annotation_col <- dplyr::arrange(annotation_col, celltype)
	
	cat("Reorder gpcr_data based on annotation_col", "\n")
	gpcr_data <- gpcr_data[, rownames(annotation_col)]
	
	cat("Check the colnames of gpcr_data is the same as rownames of annotation_col:", "\n")
	print(all.equal(colnames(gpcr_data), rownames(annotation_col)))
	
	## annotation_row information
	cat("Generate the annotation_row used in heatmap", "\n")
	annotation_row = data.frame(gpcr_family = unique_gpcr_markers$Family)
	rownames(annotation_row) <- unique_gpcr_markers$gene
	cat("Check the rownames of gpcr_data is the same as rownames of annotation_row:",  "\n")
	print(all.equal(rownames(gpcr_data), rownames(annotation_row)))
	
	### anno colors information
	cat("Generate the ann_colors used in heatmap", "\n")
	getPalette <- colorRampPalette(RColorBrewer::brewer.pal(10, "Paired"))
	colorCount <- length(celltype_levels)
	col_celltype <- getPalette(colorCount)
	names(col_celltype) <- celltype_levels
	
	getPalette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))
	colorCount<- length(sort(unique(unique_gpcr_markers$Family)))
	col_gpcr_family<- getPalette(colorCount)
	names(col_gpcr_family) <- sort(unique(unique_gpcr_markers$Family))
	# Set3; Set2; Set1; Pastel2; Pastel1; Paired; Dark2; Accent
	
	ann_colors = list(celltype = col_celltype, gpcr_family = col_gpcr_family)
	
	heatmap_input <- list(plot_data = gpcr_data, annotation_col = annotation_col, 
												annotation_row = annotation_row, ann_colors = ann_colors)
	
}


##### folder name 
cat("Set the output path\n")

folder_name <- "CellGPCR"

if (file.exists(folder_name)) {
	print("CellGPCR existed.")
}else{
	print("Create folder CellGPCR")
	dir.create(folder_name)
} 


######  1.read in data ######
cat("Step1: read the seurat object, pseudo object and gpcr list\n")

seurat_object <- readRDS(opt$seurat_object_path)
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object)
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
print(table(seurat_object@active.ident))

pseducell_object <- readRDS(opt$pseudodata_path)
DefaultAssay(pseducell_object) <- "RNA"
pseducell_object <- ScaleData(pseducell_object, features = rownames(pseducell_object))
levels(pseducell_object) <- levels(seurat_object)
print(table(pseducell_object@active.ident))

gpcr_list <- readRDS(opt$gpcr_list_path)

gpcr_list_species <- gpcr_list[gpcr_list$group == opt$org_name, ]
intersect_genename <- intersect(gpcr_list_species$genename, rownames(seurat_object))
gpcr_list_species <- gpcr_list_species[gpcr_list_species$genename %in% intersect_genename, ]		
gpcr_list_species <- gpcr_list_species[!duplicated(gpcr_list_species$genename), ]	

cat("The number of GPCR first load\n")
print(dim(gpcr_list_species))

###### 2. identify celltype specific gpcr ####
cat("Step2: identify celltype specific gpcr using original seurat object\n")
celltype_gpcr_markers <- IdentifyCellTypeTopGPCR(seurat_object = seurat_object, gpcr_list = gpcr_list_species, 
																								 min_pct = opt$min_pct, logfc_threshold = opt$logfc_threshold, p_val_adj = opt$p_val_adj)
celltype_gpcr_markers <- arrange(celltype_gpcr_markers, cluster, Family)

filename1 <- paste0("celltype_gpcr_markers_", opt$filelabel, ".csv")
write.csv(celltype_gpcr_markers, filename1)
###### 3. visualize celltype specific gpcr using original seurat object and  ####
cat("Step3: visualize celltype specific gpcr using original seurat object and pseudo object\n")
heatmap_input <- GenerateTopGPCRScaleData(seurat_object, celltype_gpcr_markers)
p_realcell <- pheatmap::pheatmap(heatmap_input$plot_data, #fontsize_row=3,
																 color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
																 breaks=seq(-1, 1, length.out = 50),
																 treeheight_row=10, treeheight_col=2,border_color='grey',
																 cluster_cols = F,cluster_rows = F,fontsize_row = 18,
																 annotation_col = heatmap_input$annotation_col,
																 annotation_row = heatmap_input$annotation_row,
																 annotation_colors = heatmap_input$ann_colors,
																 show_colnames = F, scale='none',
																 border=TRUE,
																 main = "singlecell")

heatmap_input <- GenerateTopGPCRScaleData(pseducell_object, celltype_gpcr_markers)
p_pse <- pheatmap::pheatmap(heatmap_input$plot_data, #fontsize_row=3,
														color=colorRampPalette(c("navy", "white", "firebrick3"))(50),
														breaks=seq(-1, 1, length.out = 50),
														treeheight_row=10, treeheight_col=2,border_color='grey',
														cluster_cols = F,cluster_rows = F,fontsize_row = 18,
														annotation_col = heatmap_input$annotation_col,
														annotation_row = heatmap_input$annotation_row,
														annotation_colors = heatmap_input$ann_colors,
														show_colnames = F,scale='none',
														border=TRUE,
														main ="pseudocellsize_10")
dev.off()
plots <- plot_grid(as.grob(p_realcell), as.grob(p_pse), ncol = 2)
filename2 <- paste0("select_celltype_specific_gpcr_", opt$filelabel, ".pdf")
ggsave(filename = filename2, plot = plots, width = 40, height = 20)

###### 4. visualize celltype specific gpcr using original seurat object and  ####
cat("Step3: visualize celltype specific gpcr using original seurat object and pseudo object\n")
lrdbobject  <- readRDS(opt$lrdb_path)

plots_corr <- list()
plots_net <- list()
net_res <- list()
net_adj_res <- list()

gpcr_plot <-  unique(celltype_gpcr_markers$gene)
	
## find their ligand
index <- lrdbobject$interaction$receptor %in% gpcr_plot
ligands <- intersect(unique(lrdbobject$interaction$ligand[index]), rownames(seurat_object))
	
geneinfo <- data.frame(genename = c(ligands, gpcr_plot),
													 type = c(rep("ligand", length(ligands)), rep("GPCR", length(gpcr_plot))),
													 color = c(rep("#66C2A5", length(ligands)), rep("#FC8D62", length(gpcr_plot))))
count = t(as.matrix(seurat_object@assays$RNA@counts)) # 行是细胞，列是基因
count.norm = sclink_norm(count, scale.factor = 1e6, filter.genes = FALSE, gene.names = geneinfo$genename)
row_sum <- apply(count.norm, 1, sum)
col_sum <- apply(count.norm, 2, sum)
count.norm <- count.norm[row_sum > 1, col_sum > 1]

#### calculate the correlations of genes and Infer gene co-expression networks
lda_list <- seq(0.5, 0.1, -0.05)
networks = sclink_net(expr = count.norm, ncores = 6, lda = lda_list)
corr <- networks$cor
filename3 <- paste0("GPCR_ligands_corr_heatmap_scLink_", opt$filelabel, ".pdf")
pdf(filename3, width = 12, height = 12)
print(pheatmap(corr, main = "Correaltion of GPCRs with ligand of"))
dev.off()

net_adj_res <- list()

filename4 <- paste0("GPCR_ligands_scLink_network_", opt$filelabel, ".pdf")
pdf(filename4, width = 10, height = 10)
for (i in 1:length(networks$summary)){
	print(i)
	mat <- networks$summary[[i]]$adj
	mat[!upper.tri(mat, diag = TRUE)] <- 0
	net_adj <- reshape::melt(mat)
	net_adj <- net_adj[net_adj$value==1, ]
	net_adj <- net_adj[net_adj$X1!=net_adj$X2, ]
	net_adj_res[[i]] <- net_adj
	### show ligands and gpcrs correlation
	plots_net <- graph_from_data_frame(net_adj, directed=FALSE)
	V(plots_net)$color <- geneinfo$color[geneinfo$genename %in% names(V(plots_net))]
	### show gpcrs correlation
	net_adj_gpcr <- net_adj[as.character(net_adj$X1) %in% geneinfo$genename[geneinfo$type=="GPCR"],]
	net_adj_gpcr <- net_adj_gpcr[as.character(net_adj_gpcr$X2) %in% geneinfo$genename[geneinfo$type=="GPCR"],]
	plots_net_gpcr <- graph_from_data_frame(net_adj_gpcr, directed=FALSE)
	V(plots_net_gpcr)$color <- "#FC8D62"
	
	par(mfrow = c(1,2))
	print(plot(plots_net, layout=layout.circle, node.size = 20, vertex.size=20,
			 vertex.label.cex = 1, main = paste0("lda in sclink_net:", lda_list[i], "; G:ligands; R:gpcrs")))
	print(plot(plots_net_gpcr, layout=layout.circle, node.size = 20, vertex.size=20,
			 vertex.label.cex = 1, main = paste0("lda in sclink_net:", lda_list[i],"; R:gpcrs")))
	
}
dev.off()
opar <- par(no.readonly = TRUE) #保存默认设置
par(opar) #设置为默认参数

filename5 <-  paste0("net_adj_res", opt$filelabel, ".xlsx")
names(net_adj_res) <- paste0("ida_", lda_list)
write.xlsx(net_adj_res, filename5, colNames = T, rowNames = T)

files <- c(filename1, filename2, filename3,filename4, filename5)
file.copy(files, folder_name ,overwrite = T)
file.remove(files)


print("Happy~Finished!!")
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)


