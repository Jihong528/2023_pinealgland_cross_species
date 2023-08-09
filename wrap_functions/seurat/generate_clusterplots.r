#!/home/zhengjh/miniconda3/envs/r403/bin/Rscript
# parameter
library(optparse)
library(getopt)

option_list <- list(
	make_option(c("-o", "--output"), type = "character", default = FALSE,
							action = "store", help = "This is the output directory."
	),
	make_option(c( "--seurat_object_path"), type = "character", default = FALSE,
							action = "store", help = "This is absoulte path of species1_path."
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
	),
	make_option(c("--mt_rb_pattern"), type = "character", default = FALSE,
							action = "fileprefix", help = "This is pattern used to filter out the mt and rb genes in cell type enriched markers."
	),
	make_option(c("--top_num"), type = "integer", default = 20,
							action = "fileprefix", help = "This is top_num used to do the hierarchical clustering[default 20]."
	),
	make_option(c("--cluster_method"), type = "character", default = "ward.D2",
							action = "fileprefix", help = "This is the method used in the hclust function[default ward.D2]."
	)
)

# -help 
opt = parse_args(OptionParser(option_list = option_list, 
															usage = "GenerateClusterplots.R [-o] [--seurat_object_path] [--filelabel]....."))
print(opt)

###### step0.load package, set the output path and read seurat object ####
cat("Load used R packages and functions\n")
source("/home/zhengjh/scripts/seurat/SeuratWrapperFunction.R")

timestart<-Sys.time() 

cat("Set the output path\n")
setwd(opt$output)

seurat_object <- readRDS(opt$seurat_object_path)
###### step1.generate clustering plot ####
folder_name <- "GeneralFigs"

if (file.exists(folder_name)) {
	print("GeneralFigs existed.")
	}else{
	dir.create(folder_name)
	} 

cat("Generate the clutsering plot using DimPlot function\n")

if (nrow(table(seurat_object$orig.ident)) > 1) {
	p_cluster <- DimPlot(seurat_object,  reduction = "tsne",label = T, pt.size = 0.5) + NoLegend()
	combined_cluster_plotname <-  paste0("Clustering_", opt$filelabel, ".pdf")
	ggsave(filename = combined_cluster_plotname, plot = p_cluster, width = 7, height = 7)
	p_cluster_split <- DimPlot(seurat_object,  reduction = "tsne",label = T, pt.size = 0.5, split.by = "orig.ident") + NoLegend()
	split_clustering_plot_name <- paste0("Clustering_split_group", opt$filelabel, ".pdf")
	ggsave(filename = split_clustering_plot_name, plot = p_cluster_split, width = 14, height = 7)
	
	# Copy files to 2.Cluster
	file.copy(combined_cluster_plotname, folder_name ,overwrite = T)
	file.remove(combined_cluster_plotname) 
	file.copy(split_clustering_plot_name, folder_name ,overwrite = T)
	file.remove(split_clustering_plot_name) 
	
} else{
	p_cluster <- DimPlot(seurat_object,  reduction = "tsne",label = T, pt.size = 0.5) + NoLegend()
	combined_cluster_plotname <- paste0("Clustering_", opt$filelabel, ".pdf")
	ggsave(filename = combined_cluster_plotname, plot = p_cluster, width = 7, height = 7)
	
	file.copy(combined_cluster_plotname, folder_name ,overwrite = T)
	file.remove(combined_cluster_plotname) 
}

######  step2.cell percentage plot ####

cat("Generate the cell percentage plot\n")
x <- ggplot_build(p_cluster)
info = data.frame(colour = x$data[[1]]$colour, group = x$data[[1]]$group)
info <- unique((arrange(info, group)))
cols <- as.character(info$colour)

cellnum <- as.data.frame(table(seurat_object@active.ident, 
															 seurat_object$orig.ident))
colnames(cellnum) <- c("celltype", "group", "Freq")
percentageplot <- GeomcolWrap(plot_data = cellnum, 
															x_group = cellnum$group, 
															y_value = cellnum$Freq, fill_group = cellnum$celltype, 
															fill_cols = cols, 
															title_name = "Percent of Cells from each sample",
															x_axis_name = "sample", y_axis_name = "Percent")
ggsave( paste0("Cellpercentage_", opt$filelabel, ".pdf"), percentageplot, width = 8, height = 8)
write.csv(cellnum, paste0("Cellpercentage_num_", opt$filelabel, ".csv"))


######  step3. generate the qc plot ####
cat("Generate the qc plot\n")
seurat_object$celltype_assign <- seurat_object@active.ident
p_qc <- VlnPlot(object = seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), pt.size = 0,ncol = 4, cols = cols)
ggsave( paste0("QC_gene_UMI_percent_mt_rb_", opt$filelabel, ".pdf"), p_qc, width = 12, height = 6)

######  step4. dotplot of markers and cluster based on top20 markersot ####
cat("Generate dotplot of markers and cluster based on top20 markersot ")
hclustBasedTopNumMarkers(seurat_object = seurat_object, file_prefix = opt$filelabel,
												 min_pct = opt$min_pct, logfc_threshold = opt$logfc_threshold, p_val_adj = opt$p_val_adj, 
												 mt_rb_pattern = opt$mt_rb_pattern, top_num = opt$top_num, width_num = 12, height_num = 15,
												 cluster_method = opt$cluster_method)

files <- c(paste0("Cellpercentage_", opt$filelabel, ".pdf"), 
					 paste0("Cellpercentage_num_", opt$filelabel, "csv"),
					 paste0("QC_gene_UMI_percent_mt_rb_", opt$filelabel, ".pdf"),
					 paste0("Cluster_dendrogram_", opt$cluster_method, "_", opt$top_num,"_with_dotplotmarks_", opt$filelabel,".pdf"),
					 paste0("Celltypemarker_DotPlot_", opt$filelabel, ".pdf"))
file.copy(files, folder_name ,overwrite = T)
file.remove(files) 

print("Happy~Finished!!")
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)

