#!/home/zhengjh/miniconda3/envs/r403/bin/Rscript
# parameter
library(optparse)
library(getopt)

option_list <- list(
	make_option(c("-o", "--output"), type = "character", default = FALSE,
							action = "store", help = "This is the output directory."
	),
	make_option(c( "--seurat_object_path"), type = "character", default = FALSE,
							action = "store", help = "This is absoulte path of seurat_object_path"
	),
	make_option(c("--tf_path"), type = "character", default = FALSE,
							action = "store", help = "TThis is absoulte path of tf list file download from AnimalTFDB"
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
	make_option(c("--top_num"), type = "integer", default = 10,
							action = "fileprefix", help = "This is top num TF to visualize using heatmap[default 10]."
	),
	make_option(c("--orgdb_name"), type = "character", default = "org.Dr.eg.db",
							action = "fileprefix", help = "This is top num TF to visualize using heatmap[default org.Dr.eg.db]."
	)
)

# -help 
opt = parse_args(OptionParser(option_list = option_list, 
															usage = "calculate_TFSS.R [-o] [--seurat_object_path] [--filelabel]....."))
print(opt)

###### step0.load package, set the output path and read seurat object ####
setwd(opt$output)
cat("Load used R packages and functions\n")
source("/home/zhengjh/scripts/seurat/SeuratWrapperFunction.R")

timestart<-Sys.time() 

seurat_object <- readRDS(opt$seurat_object_path)
tf_list <- read.table(opt$tf_path, sep = "\t", header = T)

###### step1. calculate TF specific score ####
cat("CalculateTFSS\n")
tf_results <- CalculateTFSS(TF_Symbol = tf_list$Symbol, seurat_object = seurat_object, min_pct = opt$min_pct, 
														logfc_threshold = opt$logfc_threshold, p_val_adj = opt$p_val_adj)

###### step2. visualize the top num TF in each cell type #####
tfss_melt <- melt(tf_results$TFSS)
colnames(tfss_melt) <- c("genename", "celltype", "score")
toptf <- tfss_melt %>% group_by(celltype) %>% top_n(n = opt$top_num, wt = score)
tf_visual <- unique(toptf$genename)
# scale函数是按列归一化，我想对TF标化所以需要先转置下
plot_data <- t(tf_results$TFSS[tf_visual, ])
tfss_zscore <- scale(plot_data)
tfss_zscore <- t(tfss_zscore)
tfss_zscore <- as.data.frame(tfss_zscore)

cat("Visualize the top num TF in each cell type using pheatmap\n")
library(pheatmap)
htmap <- pheatmap::pheatmap(tfss_zscore, scale = "none", show_rownames = F, cluster_cols = F, cluster_rows = T)
order <- htmap$tree_row$order
genelist <- row.names(tfss_zscore)[order]
cluster <- htmap$tree_row
#对聚类树进行分簇；
cut <- cutree(cluster, length(unique(tfss_melt$celltype)))
cut <- sort(cut)
genes <- names(cut)
#提取相应的表达量数据；
genesdf <- tfss_zscore[genes, ]

bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
cols_used <- c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), 
							 colorRampPalette(colors = c("white","red"))(length(bk)/2))

dev.off()
pdf(paste0("celltype_specfic_TF_withname_top", opt$top_num ,"_",opt$filelabel,  ".pdf"), 10, 50)
print(pheatmap(genesdf, scale = "none", show_rownames = T, cluster_cols = F, cluster_rows = F,
									 legend_breaks=seq(-2,2,1), colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
									 breaks=seq(-2, 2, length.out = 50)))
dev.off()

genesdf$cluster <- cut
write.csv(genesdf, paste0( "Heatmap_plotdata", opt$filelabel, ".csv"))
write.csv(tf_results$TFSS, paste0( "TFSS_table", opt$filelabel, ".csv"))
write.csv(tf_results$TF_markers, paste0("TF_makrer_", opt$filelabel, ".csv"))

####### step3: conduct enrichment analysis of the TFs ######
cat("Conduct enrichment analysis of the TFs\n")

cat("The enriched pathway of all TFs and cellty enriched TFs!")

all_enrich_pathway_results <- enrichGOByClusterProfiler(genelist = unique(tf_results$TF_markers$gene), GoOrgDbname = get(opt$orgdb_name),
																														pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
all_enrich_pathway_results$gores$celltype <- "all"
all_enrich_pathway_results$genemapid$celltype <- "all"
ego_res <- all_enrich_pathway_results$gores
plotgo_data <- list(all=all_enrich_pathway_results$gores)
gene_map <-  all_enrich_pathway_results$genemapid

for (celltype_i in unique(tf_results$TF_markers$cluster)){
	print(celltype_i)
	celltype_tf <- tf_results$TF_markers$gene[tf_results$TF_markers$cluster==celltype_i]
	print(length(celltype_tf))
	enrich_pathway_results <- enrichGOByClusterProfiler(genelist = celltype_tf, GoOrgDbname = get(opt$orgdb_name),
																											pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
	if(nrow(enrich_pathway_results$gores) != 0){
		enrich_pathway_results$gores$celltype <- celltype_i
		enrich_pathway_results$genemapid$celltype <- celltype_i
		ego_res <- rbind(ego_res,  enrich_pathway_results$gores)
		plotgo_data[[celltype_i]] <- ego_res
		gene_map <- rbind(gene_map, enrich_pathway_results$genemapid)
	}
}

write.csv(ego_res, paste0( "Enriched_go_pathway_info", opt$filelabel, ".csv"))
write.csv(gene_map, paste0( "Enriched_gene_map_to_geneid", opt$filelabel, ".csv"))
df = simplifyGOFromMultipleLists(plotgo_data, verbose=FALSE, ont = "BP", db = opt$orgdb_name, go_id_column = 2,
																 padj_column =  7, show_barplot = F, padj_cutoff = 0.001, draw_word_cloud = F, plot = F)

if(nrow(df) != 0){
	pdf(paste0( "Heatmap_gobp_keyword_enrichment_", opt$filelabel, ".pdf"), width = 10, height = 10)
	print(simplifyGOFromMultipleLists(plotgo_data, verbose=FALSE, ont = "BP", db = opt$orgdb_name, go_id_column = 2,
																		padj_column =  7, show_barplot = F, padj_cutoff = 0.001))
	dev.off()
	
	df$description <- ego_res$Description[match(df$id, ego_res$ID)]
}
write.csv(df, paste0( "simplifyGOFromMultipleLists_plot_go_pathway", opt$filelabel, ".csv"))

rm(plotgo_data)

####### step4: move the saved files to TF folder ###### 

folder_name <- "TFAnalysis"

if (file.exists(folder_name)) {
	print("TFAnalysis existed.")
}else{
	dir.create(folder_name)
}

files <- c(paste0("celltype_specfic_TF_withname_top", opt$top_num ,"_",opt$filelabel,  ".pdf"),
					 paste0( "Heatmap_plotdata", opt$filelabel, ".csv"),
					 paste0( "TFSS_table", opt$filelabel, ".csv"),
					 paste0("TF_makrer_", opt$filelabel, ".csv"),
					 paste0( "Enriched_go_pathway_info", opt$filelabel, ".csv"), 
					 paste0( "Enriched_gene_map_to_geneid", opt$filelabel, ".csv"),
					 paste0( "Heatmap_gobp_keyword_enrichment_", opt$filelabel, ".pdf"),
					 paste0( "simplifyGOFromMultipleLists_plot_go_pathway", opt$filelabel, ".csv"))
file.copy(files, folder_name ,overwrite = T)
file.remove(files) 

print("Happy~Finished!!")
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)



