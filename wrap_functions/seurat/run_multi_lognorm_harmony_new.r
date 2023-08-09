#!/home/zhengjh/miniconda3/envs/r403/bin/Rscript
# parameter

library(optparse)
library(getopt)

option_list <- list(
	make_option(c("-o", "--output"), type = "character", default = FALSE,
							action = "store", help = "This is the output directory."
	),
	make_option(c( "--seuratobject"), type = "character", default = "10X",
							action = "store", help = "This is absoulte path of seuratobject."
	),
	
	make_option(c("--min_gene"), type = "integer", default = 50,
							action = "store", help = "This is min gene number to filter, default is 50."
	),
	make_option(c("--max_gene"), type = "integer", default = 10000,
							action = "store", help = "This is max gene number to filter, default is 10000."
	),
	make_option(c("--min_ncountRNA"), type = "integer", default = 500,
							action = "store", help = "This is min ncountRNA to filter, default is 5000."
	),
	make_option(c("--max_ncountRNA"), type = "integer", default = 100000,
							action = "store", help = "This is min ncountRNA to filter, default is 100000."
	),
	make_option(c("--max_mt_percent"), type = "integer", default = 60,
							action = "store", help = "This is min mitochondria percentage  to filter, default is 60."
	),
	make_option(c("--max_rb_percent"), type = "integer", default = 60,
							action = "store", help = "This is min ribosome percentage to filter, default is 60."
	),
	make_option(c("--file_prefix"), type = "character",
							action = "store", help = "This is number of file prefix  to used in the following analysis."
	),
	make_option(c("--scale_factor"), type = "integer", default = 10000,
							action = "store", help = "This is scale factor used in Seurat NormalizeData function, default is 10000."
	),
	make_option(c("--nfeatures"), type = "integer", default = 2000,
							action = "store", help = "This is number of features used in Seurat FindVariableFeatures function, default is 2000."
	),
	make_option(c("--npc_used"), type = "integer", default = 20,
							action = "store", help = "This is number of principles to use in RunTSNE and RunUMAP function, default is 20."
	),
	
	make_option(c("--k_parameter"), type = "integer", default = 20,
							action = "store", help = "This is k parameter used in FindNeighbors function which will influence the clusering number, default is 20."
	)
)

# -help 
opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is general processing of single cell exp matrix!"))
print(opt)

source("/home/zhengjh/scripts/seurat/SeuratWrapperFunction.R")

###### step0. subset seurat object ####
timestart<-Sys.time() 

# set the output path
setwd(opt$output)

exp_count.seurat <- readRDS(opt$seuratobject)
print("Step0: Subset Seurat Object")
cat(opt$min_gene, "< gene number <", opt$max_gene, "\n")
cat(opt$min_ncountRNA, "< ncountRNA <", opt$max_ncountRNA, "\n")
cat("max_mt_percent ", opt$max_mt_percent,  "\n")
cat("max_rb_percent ", opt$max_rb_percent,  "\n")
exp_count.seurat <- SubsetSeuratData(exp_count.seurat = exp_count.seurat, 
																		 min_gene = opt$min_gene, max_gene = opt$max_gene, 
																		 min_ncountRNA = opt$min_ncountRNA, max_ncountRNA = opt$max_ncountRNA, 
																		 max_mt_percent = opt$max_mt_percent, max_rb_percent = opt$max_rb_percent)

### step1. normalize data
###### harmony ######
print("step1. normalize data")
DefaultAssay(exp_count.seurat) <- "RNA"
exp_count.seurat <- NormalizeData(exp_count.seurat, normalization.method = "LogNormalize",
																	 scale.factor = opt$scale_factor)
exp_count.seurat <- FindVariableFeatures(exp_count.seurat, selection.method = "vst",  nfeatures = opt$nfeatures)
exp_count.seurat <- ScaleData(exp_count.seurat)

# Eliminate batch effects with harmony and cell classification
exp_count.seurat <- RunPCA(exp_count.seurat, pc.genes = exp_count.seurat@var.genes, 
														npcs = opt$npc_used, verbose = FALSE)
folder_name="2.RunPCA"
if (file.exists(folder_name)){
	print("2.RunPCA existed.")
}else{
	dir.create(folder_name)
} 

selected_pcs_name <- paste0(opt$file_prefix, "_Selected_PCs.pdf")
pdf(file = selected_pcs_name, height = 8, width = 8)
print(ElbowPlot(object = exp_count.seurat, ndims = opt$npc_used, reduction = "pca"))
dev.off()

# Plot 
pca_plot_name <- paste0(opt$file_prefix, "_PCA_plot",".pdf")
pdf(pca_plot_name, 8, 8)
print(DimPlot(exp_count.seurat, reduction="pca", label = TRUE, pt.size = 0.5, group.by = "orig.ident"))
dev.off()

#Copy files to 2.RunPCA
result_files <- c(selected_pcs_name, pca_plot_name)
file.copy(result_files, folder_name ,overwrite = T)#拷贝文件
file.remove(result_files) #移除拷贝完的文件

#### step2. harmony to integrate each dataset and find clusters ####

print("step2. harmony to integrate each dataset and FindNeighbors and clusters and then run tsne and umap reduction based on harmony reduction")

cat("k_param used in Seurat FindNeighbors function is ", opt$k_param, "\n")
cat("resolution_number used in Seurat FindClusters function is ", opt$resolution_number, "\n")
cat("npc_used used in RunTSNE and RUNUMAP is ", opt$npc_used, "\n")
resolution_number <- c(seq(0.2, 2, .2))
exp_count.seurat <- PlotClusterBasedHarmony(exp_count.seurat = exp_count.seurat, 
																 file_prefix = opt$file_prefix, 
																 npc_used = opt$npc_used, 
																 k_param = opt$k_param,
																 resolution_number = resolution_number)

pdf("harmony_correct_clustree.pdf", height = 12, width = 12)
print(clustree(exp_count.seurat@meta.data, prefix = "RNA_snn_res."))
dev.off()

file.copy("harmony_correct_clustree.pdf", folder_name ,overwrite = T)#拷贝文件
file.remove("harmony_correct_clustree.pdf") #移除拷贝完的文件

rds_name <- paste0(opt$file_prefix, "_harmony_seurat.rds")
saveRDS(exp_count.seurat, file = rds_name)

print("Happy~Finished!!")
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)


