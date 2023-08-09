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
  make_option(c("--nfeatures"), type = "integer", default = 2000,
              action = "store", help = "This is number of features used in Seurat FindVariableFeatures function, default is 2000."
  ),
  make_option(c("--npc_used"), type = "integer", default = 20,
              action = "store", help = "This is number of principles to use in RunTSNE and RunUMAP function, default is 20."
  ),
  
  make_option(c("--k_parameter"), type = "integer", default = 20,
              action = "store", help = "This is k parameter used in FindNeighbors function which will influence the clusering number, default is 20."
  ),
  make_option(c("--resolution_number"), type = "double", default = 0.5,
              action = "store", help = "This is resolution number used in FindClusters function , default is 0.5"
  ),
  make_option(c("--min_pct"), type = "double", default = 0.25,
              action = "store", help = "This is min.pct used in FindAllMarkers function , default is 0.25"
  ),
  make_option(c("--logfc_threshold"), type = "double", default = 0.25,
              action = "store", help = "This is logfc.threshold used in FindClusters function , default is 0.25"
  ),
  make_option(c("--p_val_adj"), type = "double", default = 0.05,
              action = "store", help = "This is p_val_adj used to filter marker , default is 0.05"
  ),
  make_option(c("--top_num"), type = "integer", default = 10,
              action = "store", help = "This is top num marker of each cluster, default is 10"
  ),
  make_option(c("--mt_rb_pattern"), type = "character", default = "^mt-|^Rpl-|^Rps-",
              action = "store", help = "This is mt_rb_pattern used to filter mt/rb genes found in cluster markers, default is ^mt-|^Rpl-|^Rps-"
  )
  
)

# -help 
opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is general processing of single cell exp matrix!"))
print(opt)

source("/home/zhengjh/scripts/seurat/function/SeuratWrapperFunction.R")


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

#### step1. normalize each dataset ####
print("Step1: LogNormalize Data and Run PCA")
cat("scale factor used in NormalizeData function is ", opt$scale_factor, "\n")
cat("nfeatures used in Seurat FindVariableFeatures function is ", opt$nfeatures, "\n")
cat("npc_plot used in ElbowPlot is ", opt$npc_plot, "\n")

object.list <- SplitObject(exp_count.seurat, split.by = "orig.ident")
temp <- lapply(X = object.list, FUN = SCTransform)

#### step2. integrate each dataset ####
print("Step2: Use cca to integrate data")
features <- SelectIntegrationFeatures(object.list = temp, nfeatures = opt$nfeatures)
temp <- PrepSCTIntegration(object.list = temp, anchor.features = features)
object.anchors <- FindIntegrationAnchors(object.list = temp, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:opt$npc_used)
object.integrated <- IntegrateData(anchorset = object.anchors,  normalization.method = "SCT",
                                   dims = 1:opt$npc_used)
object.integrated <- RunPCA(object.integrated, npcs = opt$npc_used, verbose = FALSE)

folder_name="2.RunPCA"
if (file.exists(folder_name)){
  print("2.RunPCA existed.")
}else{
  dir.create(folder_name)
} 

selected_pcs_name <- paste0(opt$file_prefix, "_Selected_PCs.pdf")
pdf(file = selected_pcs_name, height = 8, width = 8)
print(ElbowPlot(object = object.integrated, ndims = opt$npc_used, reduction = "pca"))
dev.off()

# Plot 
pca_plot_name <- paste0(opt$file_prefix, "_PCA_plot",".pdf")
pdf(pca_plot_name, 8, 8)
print(DimPlot(object.integrated, reduction="pca", label = TRUE, pt.size = 0.5, group.by = "orig.ident"))
dev.off()

#Copy files to 2.RunPCA
result_files <- c(selected_pcs_name, pca_plot_name)
file.copy(result_files, folder_name ,overwrite = T)#拷贝文件
file.remove(result_files) #移除拷贝完的文件


#### step4. find clusters ####

print("Step4: FindNeighbors and clusters and then run tsne and umap reduction")
cat("k_param used in Seurat FindNeighbors function is ", opt$k_param, "\n")
cat("resolution_number used in Seurat FindClusters function is ", opt$resolution_number, "\n")
cat("npc_used used in RunTSNE and RUNUMAP is ", opt$npc_used, "\n")

object.integrated <- PlotCluster(exp_count.seurat = object.integrated, 
                                 file_prefix = opt$file_prefix, 
                                 npc_used = opt$npc_used, 
                                 k_param = opt$k_param,
                                 resolution_number = opt$resolution_number)

rds_name <- paste0(opt$file_prefix, ".pre.cca.integrated.final.rds")
saveRDS(object.integrated, rds_name)

#### step5. find top markers #####
print("Step5: find top markers and then map cluster with top2 marke")
cat("min.pct used in Seurat FindAllMarkers function is ", opt$min.pct, "\n")
cat("logfc.threshold used in Seurat FindAllMarkers function is ", opt$logfc.threshold, "\n")
cat("p_val_adj used to filter markers is ", opt$p_val_adj, "\n")

cluster.markers <- FindmarkerForCluster(exp_count.seurat = object.integrated, 
                                        file_prefix = opt$file_prefix, 
                                        min.pct = opt$min_pct, 
                                        logfc.threshold = opt$logfc_threshold, 
                                        p_val_adj = opt$p_val_adj,
                                        mt_rb_pattern = opt$mt_rb_pattern)

TopMarker <- TopMarkersInCluster(cluster.markers = cluster.markers, 
                                 file_prefix = opt$file_prefix, 
                                 top_num = opt$top_num)
# Top markers heatmap
P_heatmap_top_marker <- DoHeatmap(object.integrated, features = TopMarker$gene)
filename <- paste0(opt$file_prefix, "_Heatmap_plot_top_markers.pdf")
folder_name <- "4.MarkersInCluster"
pdf(filename, 30, 15)
print(P_heatmap_top_marker)
dev.off()

file.copy(filename, folder_name, overwrite = T) #copy files
file.remove(filename)

#### step6. map top2 marker to each cluster #####
print("step6: map top2 marker to each cluster")
object.integrated <- MapTop2MarkerEachCluster(exp_count.seurat = object.integrated, 
                                              cluster.markers = cluster.markers,
                                              file_prefix = opt$file_prefix)
object.integrated$celltype_assign <- object.integrated@active.ident

#### step7. save the clustering plots aftering assign top2 markers ###
print("step7. save the clustering plots aftering assign top2 markers")

p1 <- DimPlot(object.integrated, reduction = "umap",label = TRUE,
              pt.size = 0.5, label.size = 4.5, repel = TRUE, group.by = "ident") + NoLegend()
p1 <- p1 + labs(title="Clustering UMAP Reduction")
p1 <- p1 + theme(plot.title = element_text(hjust = 0.5)) 

p2 <- DimPlot(object.integrated, reduction = "umap",label = TRUE,
              pt.size = 0.5, label.size = 4.5, repel = TRUE, group.by = "ident")
p2 <- p2 + labs(title="Clustering TSNE Reduction")
p2 <- p2 + theme(plot.title = element_text(hjust = 0.5)) 

plots <- plot_grid(p1, p2, ncol=2,rel_widths = c(1.7, 2.2))

cluster_pdfname <- paste0(opt$file_prefix, "_Clustering_TSNE_UMAP.pdf")
folder_name <- "3.PlotCluster"
pdf(cluster_pdfname, width = 18, height = 6)
print(plots)
dev.off()

file.copy(cluster_pdfname, folder_name, overwrite = T)#拷贝文件
file.remove(cluster_pdfname)

rds_name <- paste0(opt$file_prefix, "_final_cca_integrate_seurat.rds")
saveRDS(object.integrated, file = rds_name)

print("Happy~Finished!!")
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)


