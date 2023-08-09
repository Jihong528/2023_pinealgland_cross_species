##------20190724 ZJH 
##------20201205 ZJH update functions 
##------20210103 ZJH update functions
##------20221025 ZJH update functions

## load function 
library(Seurat)
library(tidyr)
library(SingleCellExperiment)
library(scater)
library(MAST)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(Seurat)
library(openxlsx)
library(grid)
library(ggplotify)
library(magrittr)
library(dplyr)
library(ggsignif)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(fgsea)
library(ggplot2)
library(Seurat)
library(Matrix)
library(dplyr)
library(abind)
library(matrixStats)
library(reshape2)
library(viridis)
# library(made4)
library(RColorBrewer)
library(openxlsx)
library(clusterProfiler)
library(openxlsx)
library(org.Dr.eg.db)
library(org.Rn.eg.db)
library(org.Mmu.eg.db)
library(cowplot)
library(DOSE)
library(dplyr)
library(simplifyEnrichment)
library(gridtext)
library(clustree)

# upstream analysis 
CreateFirstFilterSeuratObject <- function(exp_count, pr_name, min_cell, 
                                          min_gene, mt_pattern, rb_pattern){
  # CreateSeuratObject : filter gene and cell by min.cells(Include features detected in at least this many cells) and min.features(Include cells where at least this many features are detected.)
  exp_count.seurat <- CreateSeuratObject(exp_count, 
                                         project = pr_name, 
                                         min.cells = min_cell,
                                         min.features = min_gene)
  Or_ident <- rep(pr_name, ncol(exp_count.seurat))
  exp_count.seurat@meta.data$orig.ident <- as.factor(Or_ident)
  exp_count.seurat@active.ident <- exp_count.seurat$orig.ident
  ## Check the mt gene and rp gene
  exp_count.seurat[["percent.mt"]] <- PercentageFeatureSet(object = exp_count.seurat, pattern = mt_pattern)
  exp_count.seurat[["percent.rb"]] <- PercentageFeatureSet(object = exp_count.seurat, pattern = rb_pattern)
  
  print(levels(exp_count.seurat@active.ident))
  print(table(exp_count.seurat$orig.ident))
  print(head(x = exp_count.seurat@meta.data, 5))
  print(dim(exp_count.seurat@meta.data))
  return(exp_count.seurat)
}

SeuratQC <- function(exp_count.seurat, name_pre){
  folder_name="QC"
  if (file.exists(folder_name)){
    print("QC existed.")
  }
  else{
    dir.create(folder_name)
  } 
  file_name=paste(name_pre,"Violinplot.pdf",sep="_")
  
  reads.drop <- isOutlier(as.numeric(exp_count.seurat$nCount_RNA), nmads = 3, type = "lower")
  feature.drop <- isOutlier(as.numeric(exp_count.seurat$nFeature_RNA), nmads = 3, type = "lower")
  qc.mito2 <- isOutlier(exp_count.seurat$percent.mt, nmads = 3, type="higher")
  qc.ribo2 <- isOutlier(exp_count.seurat$percent.rb, nmads = 3, type="higher")
  
  ### 保留下那些质量差的细胞，看看这些细胞的情况 
  low_q_cells <- (reads.drop | feature.drop | qc.mito2 | qc.ribo2)
  exp_count.seurat$lowquality_cells <- "High_qualitycells"
  exp_count.seurat$lowquality_cells[low_q_cells] <- "Low_qualitycells"
  
  plot1=FeatureScatter(object = exp_count.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "lowquality_cells")
  plot2=FeatureScatter(object = exp_count.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "lowquality_cells")
  plot3=FeatureScatter(object = exp_count.seurat, feature1 = "nCount_RNA", feature2 = "percent.rb", group.by = "lowquality_cells")
  
  pdf(file = file_name, height = 10, width = 20)
  print(VlnPlot(object = exp_count.seurat, features = c("nFeature_RNA", "nCount_RNA", 
                                                        "percent.mt","percent.rb"), ncol = 2, pt.size=0.05)) 
  
  opar <- par(no.readonly = TRUE)
  par(mfrow=c(2,2))
  print(hist(exp_count.seurat$nCount_RNA, breaks = 100, main = "Reads Count Distribution", col = "grey80", 
             xlab = "library size", ylab = "frequency",
             cex.lab = 1.4, cex.axis = 1.4))
  print(abline(v = 10000, col = "red", lwd = 2, lty = 2))
  
  print(hist(exp_count.seurat$nFeature_RNA, breaks = 100, main = "Gene Count Distribution",
             col = "grey80", xlab = "total mRNA", ylab = "frequency", cex.lab = 1.4, cex.axis = 1.4))
  print(abline(v = 1500, col = "red", lwd = 2, lty = 2))
  
  print(hist(exp_count.seurat$percent.mt, breaks = 10, main = "Mitochondrial Percentage Distribution",
             col = "grey80", xlab = "total mRNA", ylab = "frequency", cex.lab = 1.4, cex.axis = 1.4))
  print(abline(v = 20, col = "red", lwd = 2, lty = 2))
  
  print(hist(exp_count.seurat$percent.rb, breaks = 10, main = "Ribosome Percentage Distribution",
             col = "grey80", xlab = "total mRNA", ylab = "frequency", cex.lab = 1.4, cex.axis = 1.4))
  print(abline(v = 20, col = "red", lwd = 2, lty = 2))
  par(opar)
  
  print(plot1 + plot2 + plot3)
  dev.off()
  
  #copy files to QC
  result_files <- c(file_name)
  file.copy(result_files, folder_name ,overwrite = T)
  file.remove(result_files) 
  print("Based on the Violin plot, you can filter the cell by the nFeature_RNA, nCount_RNA and/or  percent.mt rb")
  return(exp_count.seurat)
  }


SubsetSeuratData <- function(exp_count.seurat, min_gene, max_gene,
                             min_ncountRNA,  max_ncountRNA,
                             max_mt_percent, max_rb_percent){
  exp_count.seurat <- subset(exp_count.seurat,  subset = nCount_RNA < max_ncountRNA &  nCount_RNA > min_ncountRNA & nFeature_RNA < max_gene & nFeature_RNA > min_gene  & percent.mt < max_mt_percent & percent.rb < max_rb_percent )
  # QC
  qc_pdf_name <- paste("filtered_by_gene", min_gene, "UMI", min_ncountRNA, "percent.mt", max_mt_percent, "percent.rb", max_rb_percent, sep="_")
  SeuratQC(exp_count.seurat, qc_pdf_name)
  dim(exp_count.seurat)
  print(levels(exp_count.seurat$orig.ident))
  return(exp_count.seurat)
}


SelectPCsByLogNormalize <- function(exp_count.seurat, file_prefix, scale_factor, 
                                    nfeatures, npc_plot){
  
  folder_name="RunPCA"
  if (file.exists(folder_name)){
    print("RunPCA existed.")
  }
  else{
    dir.create(folder_name)
  } 
  
  DefaultAssay(exp_count.seurat) <- "RNA"
  # Normalize Find VariableGene and Scale Data
  exp_count.seurat <- NormalizeData(exp_count.seurat, normalization.method = "LogNormalize",
                                    scale.factor = scale_factor)
  exp_count.seurat <- FindVariableFeatures(exp_count.seurat, selection.method = "vst", nfeatures = nfeatures)
  exp_count.seurat <- ScaleData(exp_count.seurat, vars.to.regress = c("percent.rb","percent.mt"))
  
  # RunPCA
  exp_count.seurat <- RunPCA(exp_count.seurat, pc.genes = exp_count.seurat@var.genes, 
                             npcs = npc_plot, verbose = FALSE,)
  
  # Select PCs:1:30pcs
  selected_pcs_name <- paste0(file_prefix, "_Selected_PCs.pdf")
  pdf(file = selected_pcs_name, height = 8, width = 8)
  print(ElbowPlot(object = exp_count.seurat, ndims = npc_plot, reduction = "pca"))
  dev.off()
  
  # Plot 
  pca_plot_name <- paste0(file_prefix, "_PCA_plot",".pdf")
  pdf(pca_plot_name,8,8)
  print(DimPlot(exp_count.seurat, reduction="pca", label = TRUE, pt.size=0.5, group.by="ident"))
  dev.off()
  
  #Copy files to RunPCA
  result_files <- c(selected_pcs_name,pca_plot_name)
  file.copy(result_files, folder_name ,overwrite = T)#拷贝文件
  file.remove(result_files) #移除拷贝完的文件
  
  return(exp_count.seurat)
}

SelectPCsBySCTransform <- function(exp_count.seurat, file_prefix, npc_plot){
  
  folder_name="RunPCA"
  if (file.exists(folder_name)){
    print("RunPCA existed.")
  }
  else{
    dir.create(folder_name)
  } 
  
  # Normalize Find VariableGene and Scale Data
  exp_count.seurat <- SCTransform(exp_count.seurat, vars.to.regress = c("percent.rb","percent.mt"))
  # RunPCA
  exp_count.seurat <- RunPCA(exp_count.seurat, verbose = FALSE,  npcs = npc_plot)
  
  # Select PCs:1:30pcs
  selected_pcs_name <- paste0(file_prefix, "_Selected_PCs.pdf")
  pdf(file = selected_pcs_name, height = 8, width = 8)
  print(ElbowPlot(object = exp_count.seurat, ndims = npc_plot, reduction = "pca"))
  dev.off()
  
  # Plot 
  pca_plot_name <- paste0(file_prefix, "_PCA_plot",".pdf")
  pdf(pca_plot_name, 8, 8)
  print(DimPlot(exp_count.seurat,reduction="pca", label = TRUE, pt.size= 0.5,group.by="ident",shape.by="orig.ident"))
  dev.off()
  
  #Copy files to RunPCA
  result_files <- c(selected_pcs_name,pca_plot_name)
  file.copy(result_files, folder_name ,overwrite = T)#拷贝文件
  file.remove(result_files) #移除拷贝完的文件
  
  return(exp_count.seurat)
}

# RunUMAP and RunTSNE and Then find cluster
PlotCluster <- function(exp_count.seurat, file_prefix, npc_used, k_param, resolution_number){
  
  folder_name="PlotCluster"
  if (file.exists(folder_name)){
    print("PlotCluster existed.")
  }
  else{
    dir.create(folder_name)
  } 
  
  # FinderNeighbors 
  exp_count.seurat <- FindNeighbors(exp_count.seurat, dims = 1:npc_used, verbose = FALSE, k.param = k_param)
  exp_count.seurat <- FindClusters(exp_count.seurat, verbose = FALSE, resolution = resolution_number)
  
  
  # RunTSNE and RunUMAP
  exp_count.seurat <- RunTSNE(exp_count.seurat, dims = 1:npc_used, verbose = FALSE, check_duplicates = FALSE)
  exp_count.seurat <- RunUMAP(exp_count.seurat, dims = 1:npc_used, verbose = FALSE)
  
  # RenameIdents from zero to one
  levels_define <- as.numeric(levels(exp_count.seurat))
  new.cluster.ids <- levels_define + 1
  names(new.cluster.ids) <- levels(exp_count.seurat)
  exp_count.seurat <- RenameIdents(exp_count.seurat, new.cluster.ids)
  
  combined_cluster_plotname <- paste0(file_prefix, "_combined_cluster_resolution_",resolution_number,".pdf")
  pdf(combined_cluster_plotname,7,7)
  print(DimPlot(exp_count.seurat,reduction="umap",label = TRUE, pt.size = 0.5, 
                label.size = 4.5, repel=TRUE, group.by="ident"))
  print(DimPlot(exp_count.seurat,reduction="tsne", label = TRUE, pt.size = 0.5, 
                label.size = 4.5, repel=TRUE, group.by="ident"))
  dev.off()
  
  split_clustering_plot_name <- paste0(file_prefix, "_split_clustering_resolution_", resolution_number, ".pdf")
  pdf(split_clustering_plot_name,14,7)
  print(DimPlot(exp_count.seurat,reduction="umap",label = TRUE, pt.size = 0.5, 
                label.size = 4.5, group.by="ident", repel=TRUE, split.by="orig.ident"))
  print(DimPlot(exp_count.seurat,reduction="tsne", label = TRUE, pt.size = 0.5, 
                label.size = 4.5, group.by="ident", repel=TRUE, split.by="orig.ident")) 
  dev.off()
  
  # Copy files to 2.Cluster
  file.copy(combined_cluster_plotname, folder_name ,overwrite = T)
  file.remove(combined_cluster_plotname) 
  file.copy(split_clustering_plot_name, folder_name ,overwrite = T)
  file.remove(split_clustering_plot_name) 
  return(exp_count.seurat)
}

PlotClusterBasedHarmony <- function(exp_count.seurat, file_prefix, npc_used, k_param, resolution_number){
  
  folder_name="PlotCluster"
  if (file.exists(folder_name)){
    print("PlotCluster existed.")
  }
  else{
    dir.create(folder_name)
  } 
  
  # remove batch effect
  
  library(harmony)
  exp_count.seurat <- exp_count.seurat %>% RunHarmony("orig.ident", plot_convergence = TRUE)
  
  # RunTSNE and RunUMAP
  exp_count.seurat <- RunTSNE(exp_count.seurat, reduction = "harmony", dims = 1:npc_used, verbose = FALSE, check_duplicates = FALSE)
  exp_count.seurat <- RunUMAP(exp_count.seurat, reduction = "harmony", dims = 1:npc_used, verbose = FALSE)
  
  # FinderNeighbors 
  exp_count.seurat <- FindNeighbors(exp_count.seurat, dims = 1:npc_used, reduction = "harmony", verbose = FALSE, k.param = k_param)
  exp_count.seurat <- FindClusters(exp_count.seurat, verbose = FALSE, resolution = resolution_number)
  

  # RenameIdents from zero to one
  levels_define <- as.numeric(levels(exp_count.seurat))
  new.cluster.ids <- levels_define + 1
  names(new.cluster.ids) <- levels(exp_count.seurat)
  exp_count.seurat <- RenameIdents(exp_count.seurat, new.cluster.ids)
  
  combined_cluster_plotname <- paste0(file_prefix, "_combined_cluster_resolution_",resolution_number,".pdf")
  pdf(combined_cluster_plotname,7,7)
  print(DimPlot(exp_count.seurat,reduction="umap",label = TRUE, pt.size = 0.5, 
                label.size = 4.5, repel=TRUE, group.by="ident"))
  print(DimPlot(exp_count.seurat,reduction="tsne", label = TRUE, pt.size = 0.5, 
                label.size = 4.5, repel=TRUE, group.by="ident"))
  dev.off()
  
  split_clustering_plot_name <- paste0(file_prefix, "_split_clustering_resolution_", resolution_number, ".pdf")
  pdf(split_clustering_plot_name,14,7)
  print(DimPlot(exp_count.seurat,reduction="umap",label = TRUE, pt.size = 0.5, 
                label.size = 4.5, group.by="ident", repel=TRUE, split.by="orig.ident"))
  print(DimPlot(exp_count.seurat,reduction="tsne", label = TRUE, pt.size = 0.5, 
                label.size = 4.5, group.by="ident", repel=TRUE, split.by="orig.ident")) 
  dev.off()
  
  # Copy files to 2.Cluster
  file.copy(combined_cluster_plotname, folder_name ,overwrite = T)
  file.remove(combined_cluster_plotname) 
  file.copy(split_clustering_plot_name, folder_name ,overwrite = T)
  file.remove(split_clustering_plot_name) 
  return(exp_count.seurat)
}

# Find markers for each cluster
FindmarkerForCluster <- function(exp_count.seurat, file_prefix, min.pct, logfc.threshold, p_val_adj, mt_rb_pattern){
  folder_name="MarkersInCluster"
  if (file.exists(folder_name)){
    print("MarkersInCluster file existed.")
  }
  else{
    dir.create(folder_name)
  } 
  
  cluster.markers <- FindAllMarkers(object = exp_count.seurat, only.pos = TRUE, 
                                    min.pct = min.pct, logfc.threshold = logfc.threshold)
  index <- cluster.markers$p_val_adj < p_val_adj
  cluster.markers <- cluster.markers[index,]
  # remove mt and rb gene
  index <- grep(mt_rb_pattern, cluster.markers$gene)
  if (length(index)  > 0){
    cluster.markers <- cluster.markers[-index, ]
  }
  
  save_name <- paste0(file_prefix,"_MarkersInClusters_min.pct_",min.pct, "_logfc.threshold_", logfc.threshold,"_p_val_adj_", ".csv")
  write.csv(cluster.markers,save_name)
  
  file.copy(save_name, folder_name,overwrite = T)
  file.remove(save_name)
  return(cluster.markers)
}

# Top markers for each cluster
TopMarkersInCluster <- function(cluster.markers, file_prefix, top_num){
  library(dplyr)
  folder_name="MarkersInCluster"
  if (file.exists(folder_name)){
    print("MarkersInCluster file existed.")
  }
  else{
    dir.create(folder_name)
  } 
  #将readsCountSM.markers传给group_by，group_by按cluster 排序，再将结果传给top_n，top_n按avg_logFC排序，显示每个类中的前两个
  top_marker <- cluster.markers %>% group_by(cluster) %>% top_n(n = top_num, wt = avg_log2FC)
  file_name=paste(file_prefix, "_top_marker", top_num,".csv",sep="")
  write.csv(top_marker, file =file_name)
  file.copy(file_name, folder_name,overwrite = T)
  file.remove(file_name)
  return(top_marker)
}

# Rename each cluster with top2 markers
MapTop2MarkerEachCluster <- function(exp_count.seurat, cluster.markers, file_prefix){
  library(dplyr)
  library(plyr)
  exp_count.seurat$seurat_clusters <- exp_count.seurat@active.ident
  
  top2 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  top2 <- top2 %>%  group_by(cluster) %>%  
    dplyr::mutate(markers = paste0(gene, collapse = "/")) %>% dplyr::slice(1)  
  marker.names <- top2$markers
  current.cluster.ids <- as.character(1:(length(unique(Idents(exp_count.seurat))))) 
  new.cluster.ids <- marker.names
  Idents(exp_count.seurat) <- plyr::mapvalues(Idents(exp_count.seurat),
                                              from = current.cluster.ids, 
                                              to = new.cluster.ids)
  print(Idents(exp_count.seurat))
  print(exp_count.seurat@active.ident)
  exp_count.seurat$celltype_assign <- exp_count.seurat@active.ident
  return(exp_count.seurat)
}

MapCelltypeByMarkers <- function(exp_count.seurat, markercelltype_df){
  clustavg <- AverageExpression(exp_count.seurat, features = markercelltype_df$markers)
  clustdf <- as.data.frame(clustavg$RNA)
  clustdf <- t(clustdf)
  clustdf <- as_tibble(clustdf, rownames = "cluster")
  print(head(clustdf))
  print(dim(clustdf))
  
  df2 <- clustdf %>% dplyr::select("cluster", all_of(markercelltype_df$markers))
  
  newdf <- df2 %>% gather(markers, max_value, 2:(length(markercelltype_df$markers) + 1)) %>% 
    group_by(cluster) %>% slice(which.max(max_value)) 
  newdf <- left_join(x = newdf, y = df2, by	= "cluster")

  newdf <- left_join(x = newdf, y = markercelltype_df, by	= "markers")
  
  #exp_count.seurat@active.ident <- exp_count.seurat$seurat_clusters
  Idents(exp_count.seurat) <- plyr::mapvalues(Idents(exp_count.seurat), 
                                              from = newdf$cluster, to = newdf$celltype)
  exp_count.seurat$celltype_maxvalue <- Idents(exp_count.seurat)
  return(exp_count.seurat)
}

# harmony used to remove batch effect
RemoveBftByHarmony <- function(exp_count.seurat, npc_used, k_parameter, resolution_number){
  library(harmony)
  exp_count.seurat <- exp_count.seurat %>% RunHarmony("orig.ident", plot_convergence = TRUE)
  
  exp_count.seurat <- exp_count.seurat %>% 
  RunTSNE(reduction = "harmony", dims = 1:npc_used) %>% 
  RunUMAP(reduction = "harmony", dims = 1:npc_used) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:npc_used, k.param = k_parameter) %>% 
  FindClusters(resolution = resolution_number) %>% 
  identity()
  
  new.cluster.ids <- as.numeric(levels(exp_count.seurat)) + 1
  names(new.cluster.ids) <- levels(exp_count.seurat)
  exp_count.seurat <- RenameIdents(exp_count.seurat, new.cluster.ids)
  return(exp_count.seurat)
}

## Identify doublets
IdentifyDoubletsByDoubletFinder <- function(exp_count.seurat, npc_used, nExp_poi_number){
  library(DoubletFinder)
  sweep.res.list <- paramSweep_v3(exp_count.seurat, PCs = 1:npc_used)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  # 0.2
  
  ## Homotypic Doublet Proportion Estimate 
  annotations <- exp_count.seurat@active.ident
  homotypic.prop <- modelHomotypic(annotations)
  
  nExp_poi <- round( nExp_poi_number *ncol(exp_count.seurat@assays$RNA@data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies 
  seurat_filterDouble <- doubletFinder_v3(exp_count.seurat, PCs = 1:npc_used, 
                                          pN = 0.25, pK = mpK, nExp = nExp_poi,
                                          reuse.pANN = FALSE, sct = T)
  index <- grep("^DF.classifications", colnames(seurat_filterDouble@meta.data))
  # Aggregated_seurat <- seurat_filterDouble
  exp_count.seurat$doublets <- seurat_filterDouble@meta.data[ , index]
  return(exp_count.seurat)
}

## parallel function 
SbatchParallel <- function(cl_cores, code_path){
  library(parallel)
  cl <- makeCluster(cl_cores)
  clusterEvalQ(cl, source(file = code_path))
  stopCluster(cl)
}


SeuratIdentifyCelltype <- function(exp_count.seurat, ref.seurat){
  
  exp_count.seurat.anchors <- FindTransferAnchors(reference =ref.seurat, query = exp_count.seurat, 
                                              dims = 1:30, project.query = T, k.filter=150, k.anchor=5)
  predictions <- TransferData(anchorset = exp_count.seurat.anchors, refdata = ref.seurat$celltype, 
                              dims = 1:30)
  exp_count.seurat$celltype_Seurat <- predictions$predicted.id
  return(exp_count.seurat)
}



# Downstream Analysis

BarplotWrap <- function(plot_data, stat_mode, x_group, y_value, 
                                fill_group, fill_cols, max_y, y_steps,
                                title_name, x_axis_name, y_axis_name){
  # plot_data: three columns dataframe(x_group, y_value, fill_group)
  # stat_mode: "identity", stack
library(ggplot2)
theme_used <- theme_bw() + theme(panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), 
                                 axis.line = element_line(colour = "black"),
                                 panel.border = element_blank(), 
                                 aspect.ratio = 1)
theme_used <- theme_used + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, 
                                                            hjust = 0.5, size = 10),
                                 axis.text.y = element_text(size = 10),
                                 axis.title.x = element_text(size = 10),
                                 axis.title.y = element_text(size = 10),
                                 plot.title = element_text(hjust = 0.5, size = 15))

barplot_p <- ggplot(data = plot_data, aes(x = x_group, y = y_value, fill = fill_group)) 
barplot_p <- barplot_p + geom_bar(position = stat_mode, stat = "identity") + theme_used
# change fill colors
barplot_p <- barplot_p + scale_fill_manual(values = fill_cols) 

barplot_p <- barplot_p + scale_y_continuous(limits=c(0, max_y), breaks=seq(0, max_y, y_steps), 
                                  expand = c(0,0))
barplot_p <- barplot_p + ggtitle(title_name) + xlab(x_axis_name) + ylab(y_axis_name)
return(barplot_p)
}

ViolinplotWrap <- function(plot_data, x_group, y_value, 
                        fill_group, fill_cols, max_y, y_steps,
                        title_name, x_axis_name, y_axis_name){
  # plot_data: three columns dataframe(x_group, y_value, fill_group)
  # stat_mode: "identity", stack
  library(ggplot2)
  theme_used <- theme_bw() + theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), 
                                   axis.line = element_line(colour = "black"),
                                   panel.border = element_blank(), 
                                   aspect.ratio = 1)
  theme_used <- theme_used + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, 
                                                              hjust = 0.5, size = 10),
                                   axis.text.y = element_text(size = 10),
                                   axis.title.x = element_text(size = 10),
                                   axis.title.y = element_text(size = 10),
                                   plot.title = element_text(hjust = 0.5, size = 15))
  
  violinplot_p <- ggplot(data = plot_data, aes(x = x_group, y = y_value, fill = fill_group)) 
  violinplot_p <- violinplot_p + geom_violin(show.legend = FALSE)  + theme_used
  # change fill colors
  violinplot_p <- violinplot_p + scale_fill_manual(values = fill_cols) 
  
  violinplot_p <- violinplot_p + scale_y_continuous(limits=c(0, max_y), breaks=seq(0, max_y, y_steps), 
                                              expand = c(0,0))
  violinplot_p <- violinplot_p + ggtitle(title_name) + xlab(x_axis_name) + ylab(y_axis_name)
  return(violinplot_p)
}

GeomcolWrap <- function(plot_data, x_group, y_value, 
                           fill_group, fill_cols,
                           title_name, x_axis_name, y_axis_name){
  # plot_data: three columns dataframe(x_group, y_value, fill_group)
  # stat_mode: "identity", stack
  library(ggplot2)
  theme_used <- theme_bw() + theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), 
                                   axis.line = element_line(colour = "black"),
                                   panel.border = element_blank(), 
                                   aspect.ratio = 1)
  theme_used <- theme_used + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, 
                                                              hjust = 0.5, size = 10),
                                   axis.text.y = element_text(size = 10),
                                   axis.title.x = element_text(size = 10),
                                   axis.title.y = element_text(size = 10),
                                   plot.title = element_text(hjust = 0.5, size = 15))
  
  plot_p <- ggplot(data = plot_data, aes(x = x_group, y = y_value, fill = fill_group)) 
  plot_p <- plot_p + geom_col(position = "fill")  + theme_used
  # change fill colors
  plot_p <- plot_p + scale_fill_manual(values = fill_cols) 
  
  plot_p <- plot_p + scale_y_continuous(expand = c(0,0))
  plot_p <- plot_p + ggtitle(title_name) + xlab(x_axis_name) + ylab(y_axis_name)
  return(plot_p)
}


# DE genes

FindDEGsFromMasttest <- function(exp_count.seurat, cell_type, treat, control,  min_pct_cell, file_prefix){
  library(MAST)
  results_DEGs_mast <- list()
  
  exp_count.seurat$orig.ident <- as.character(exp_count.seurat$orig.ident)
  print("Find DEGs from all cells using MAST test between treat and control")
  cat("treat group is", treat, "\n")
  cat("control group is", control, "\n")
  DefaultAssay(exp_count.seurat) <- "RNA"
  
  Idents(exp_count.seurat) <- exp_count.seurat$group
  print(table(Idents(exp_count.seurat)))
  bulk_MAST<- FindMarkers(exp_count.seurat, slot = "data",
                          logfc.threshold = 0, min.pct = min_pct_cell,
                          ident.1 = treat, ident.2 = control, only.pos = FALSE,
                          verbose = T, test.use = "MAST", latent.vars = c("nCount_RNA", "orig.ident"))
  results_DEGs_mast[["bulk_MAST"]] <- bulk_MAST
  
  ## in order to subset object of different celltype
  Idents(exp_count.seurat) <- exp_count.seurat$celltype_assign
  print("All cells finished!")
  print("Find DEGs from the celltypes you used")
  for (i in 1:length(cell_type)){
    ## find DEGs from treat vs control
    print(cell_type[i])
    print("Find DEGs using MAST test from each celltype between treat and control.....")
    subset_seurat <- subset(exp_count.seurat, idents = c(cell_type[i]))
    print(table(Idents(subset_seurat)))
    
    ## in order to find DEGs from different group
    Idents(subset_seurat) <- subset_seurat$group
    print(table(Idents(subset_seurat)))
    
    dd <- as.data.frame(table(Idents(subset_seurat)))
    if(sum(dd$Freq > 3)){
    DefaultAssay(subset_seurat) <-"RNA"
    cat("treat group is", treat, "\n")
    cat("control group is", control, "\n")
    print("Find DEGs of this celltype from different group")
    print(table(Idents(subset_seurat)))
    DEGs_treat_VS_Control_mast <- FindMarkers(subset_seurat, slot = "data",
                                              logfc.threshold = 0, min.pct = min_pct_cell,
                                              ident.1 = treat, ident.2 = control, test.use = "MAST",
                                              only.pos = FALSE, verbose = T,
                                              latent.vars = c("nCount_RNA", "orig.ident"))
    DEGs_treat_VS_Control_mast$cell_type <- cell_type[i]
    DEGs_treat_VS_Control_mast <- as_tibble(DEGs_treat_VS_Control_mast, rownames = "gene")
    results_DEGs_mast[[cell_type[i]]] <- DEGs_treat_VS_Control_mast
    
    print("finished!")
    }else{
      print("One of the number of this group of this celltype is smaller than 3")
      next
    }
  }
  
  ## save results
  write.xlsx(results_DEGs_mast, file = paste0(file_prefix, "_mast_DEGs_Treat_VS_Control.xlsx"),
             col.names = T, row.names = T)
  saveRDS(results_DEGs_mast, paste0(file_prefix, "_mast_DEGs_Treat_VS_Control.rds"))
}

FindDEGsFromWilcoxtest <- function(exp_count.seurat, cell_type, treat, control,  min_pct_cell, file_prefix){
	library(MAST)
	results_DEGs_mast <- list()
	
	exp_count.seurat$orig.ident <- as.character(exp_count.seurat$orig.ident)
	print("Find DEGs from all cells using wilcox test between treat and control")
	cat("treat group is", treat, "\n")
	cat("control group is", control, "\n")
	DefaultAssay(exp_count.seurat) <- "RNA"
	
	Idents(exp_count.seurat) <- exp_count.seurat$group
	print(table(Idents(exp_count.seurat)))
	bulk_MAST<- FindMarkers(exp_count.seurat, slot = "data",
													logfc.threshold = 0, min.pct = min_pct_cell,
													ident.1 = treat, ident.2 = control, only.pos = FALSE,
													verbose = T, test.use = "wilcox")
	bulk_MAST$gene <- rownames(bulk_MAST)
	results_DEGs_mast[["bulk_MAST"]] <- bulk_MAST
	
	## in order to subset object of different celltype
	Idents(exp_count.seurat) <- exp_count.seurat$celltype_assign
	print("All cells finished!")
	print("Find DEGs from the celltypes you used")
	for (i in 1:length(cell_type)){
		## find DEGs from treat vs control
		print(cell_type[i])
		print("Find DEGs using MAST test from each celltype between treat and control.....")
		subset_seurat <- subset(exp_count.seurat, idents = c(cell_type[i]))
		print(table(Idents(subset_seurat)))
		
		## in order to find DEGs from different group
		Idents(subset_seurat) <- subset_seurat$group
		print(table(Idents(subset_seurat)))
		
		dd <- as.data.frame(table(Idents(subset_seurat)))
		if(sum(dd$Freq > 3)){
			DefaultAssay(subset_seurat) <-"RNA"
			cat("treat group is", treat, "\n")
			cat("control group is", control, "\n")
			print("Find DEGs of this celltype from different group")
			print(table(Idents(subset_seurat)))
			DEGs_treat_VS_Control_mast <- FindMarkers(subset_seurat, slot = "data",
																								logfc.threshold = 0, min.pct = min_pct_cell,
																								ident.1 = treat, ident.2 = control, test.use = "wilcox",
																								only.pos = FALSE, verbose = T)
			DEGs_treat_VS_Control_mast$cell_type <- cell_type[i]
			DEGs_treat_VS_Control_mast <- as_tibble(DEGs_treat_VS_Control_mast, rownames = "gene")
			results_DEGs_mast[[cell_type[i]]] <- DEGs_treat_VS_Control_mast
			
			print("finished!")
		}else{
			print("One of the number of this group of this celltype is smaller than 3")
			next
		}
	}
	
	## save results
	write.xlsx(results_DEGs_mast, file = paste0(file_prefix, "_wilcox_DEGs_Treat_VS_Control.xlsx"),
						 col.names = T, row.names = T)
	saveRDS(results_DEGs_mast, paste0(file_prefix, "_wilcox_DEGs_Treat_VS_Control.rds"))
}


RunFgsea <- function(DE_dat, pathway_used, cell_name, padj_cutoff) {
  res2 <- dplyr::select(DE_dat, gene, avg_log2FC)
  res2 <- arrange(res2, desc(avg_log2FC))
  gene_value <- tibble::deframe(res2) # to vector
  head(gene_value, 20)
  daty <- as_tibble(fgseaMultilevel(pathways = pathway_used, stats = gene_value,
                                    nPermSimple= 1000))
  daty$cell_type<-cell_name
  daty <- dplyr::filter(daty, padj < padj_cutoff)
  daty <- mutate(daty, leading_edge_count = lengths(daty$leadingEdge))
  return(daty)
}

enrichGOKEGGByClusterProfiler <- function(genelist, GoOrgDbname, pvalueCutoff, qvalueCutoff, pAdjustMethod, KEGGorganism){
	gene_setdiff <- setdiff(gene, keys(GoOrgDbname, keytype="SYMBOL"))
	print("Using bitr function to convert gene symbol to ENTREZID:")
	gene <- bitr(genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb=GoOrgDbname) # groupGO默认使用"ENTREZID"作为输入数据，因此首先将“SYMBOL”转化为“ENTREZID”，以防出错
	print("Conducting enrichGO:")
	ego <- enrichGO(gene=gene$ENTREZID, keyType = "ENTREZID", OrgDb = GoOrgDbname, ont = "ALL", 
									pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, readable = TRUE) # 是否将gene ID转换到 gene symbol)
	print("Conducting enrichKEGG:")
	kk <- enrichKEGG(gene = as.numeric(gene$ENTREZID), organism = KEGGorganism, pvalueCutoff = pvalueCutoff, 
									 qvalueCutoff = qvalueCutoff, pAdjustMethod =pAdjustMethod)
	res <- list(gores=as.data.frame(ego), keggres = as.data.frame(kk), gene_unmap = gene_setdiff, genemapid=gene)
	return(res)
}

enrichGOByClusterProfiler <- function(genelist, GoOrgDbname, pvalueCutoff, qvalueCutoff, pAdjustMethod){
	gene_setdiff <- setdiff(genelist, keys(GoOrgDbname, keytype="SYMBOL"))
	print("Using bitr function to convert gene symbol to ENTREZID:")
	gene <- bitr(genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb=GoOrgDbname) # groupGO默认使用"ENTREZID"作为输入数据，因此首先将“SYMBOL”转化为“ENTREZID”，以防出错
	print("Conducting enrichGO:")
	ego <- enrichGO(gene=gene$ENTREZID, keyType = "ENTREZID", OrgDb = GoOrgDbname, ont = "all", 
									pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, readable = TRUE) # 是否将gene ID转换到 gene symbol)
	print("Finished!")
	res <- list(gores=as.data.frame(ego), genemapid=gene)
	return(res)
}

CalculateCelltypeAvgasinhFC <- function(seurat_object, commGenes, treat_label, control_label, file_prefix){
  
  ## note: 
  # seurat_object$orig.ident is the group you want to use
  # seurat_object@active.ident is the celltype 
  ## a. find asinhFC by celltype 
  
  c.avg <- c(NA)
  c.avg.list <- list()
  celltype <- levels(seurat_object@active.ident)
  
  for (i in celltype){
    temp_object <- subset(seurat_object, ident = i)
    c.avg <- AverageExpression(temp_object, features  = commGenes, 
                               group.by = "orig.ident", assays = "RNA", slot = "data")
    c.avg.list[[i]] <- c.avg                         
    setNames(c.avg.list, paste0('c.',i,'.avg'))
    }
  
  ### b. add asinhFC to clustAvg 
  c.avg.asinhFC <- c(NA)
  c.avg.asinhFC.list <- list()
  for (i in celltype){
    treat <- asinh(expm1(c.avg.list[[i]]$RNA[ , treat_label])) ## 个人理解expm1是把平均值放大，好计算
    control <- asinh(expm1(c.avg.list[[i]]$RNA[ , control_label]))
    asinhFC = treat - control
    c.avg.asinhFC <- cbind(c.avg.list[[i]]$RNA, asinhFC)
    #assign(paste0('c',i,'.avg.asinhFC'), c.avg.asinhFC)
    c.avg.asinhFC.list[[i]] <- c.avg.asinhFC                        
    setNames(c.avg.asinhFC.list, paste0('c.',i,'.avg.asinhFC'))
  }
  
  rm(asinhFC, c.avg.asinhFC, i, treat, control, c.avg)
  
  ### c.merge asinhFC for all clusters into a single data.frame and save
  c.clustAvg.asinhFC <- do.call("cbind", lapply(c.avg.asinhFC.list, function(x) (x[,3])))
  c.clustAvg.asinhFC <- cbind(rownames(c.avg.asinhFC.list[[5]]), c.clustAvg.asinhFC )
  c.clustAvg.asinhFC <- as.data.frame(c.clustAvg.asinhFC[,-1], row.names=c.clustAvg.asinhFC[,1])
  
  names <- c(NA)
  names.list <- list()
  for (i in celltype){names <- paste0('c', rep(i),'.asinhFC')
  assign(paste0('c',i,'.asinhFC'), names)
  names.list[[i]] <- names                        
  setNames(names.list, paste0('c',i,'.asinhFC'))}
  
  colnames(c.clustAvg.asinhFC)<- names.list[1:12]
  ## transvert c.clustAvg.asinhFC factor-->numeric
  c.clustAvg.asinhFC.indx <- sapply(c.clustAvg.asinhFC, is.factor) 
  c.clustAvg.asinhFC[c.clustAvg.asinhFC.indx] <- lapply(c.clustAvg.asinhFC[c.clustAvg.asinhFC.indx], 
                                                        function(x) as.numeric(as.character(x)))
  
  saveRDS(c.clustAvg.asinhFC, paste0(file_prefix, '.clustAvg.master.asinhFC.rds'))
  write.table(c.clustAvg.asinhFC, paste0(file_prefix, '.clustAvg.master.asinhFC.txt'), 
                                         sep="\t", col.names=NA)
  return(c.clustAvg.asinhFC)
}


CalculateCelltypeGenePct <- function(seurat_object, commGenes, file_prefix){

# create matrix of percent positive cells per cluster
rawData <- as.matrix(seurat_object@assays$RNA@counts)[commGenes, ]
celltype <- levels(seurat_object@active.ident)

# create list of cell names for each cluster 
allCells.ident <- as.matrix(seurat_object@active.ident)
colnames(allCells.ident) <- "celltypename"

cells <- c(NA)
cells.list <- list()
for(i in celltype)
{cells <- names(allCells.ident[allCells.ident[ , "celltypename"] == i, ] )
#assign(paste0('c',i,'.cells'), cells)
cells.list[[i]] <- cells                         
setNames(cells.list, paste0('c',i,'.cells'))}

# slice rawData by celltype 
cells.rawData <- c(NA)
cells.rawData.list <- list()
for(i in celltype){cells.rawData <- rawData[ , cells.list[[i]]]
#assign(paste0('c',i,'.cells.rawData'), cells.rawData)
cells.rawData.list[[i]] <- cells.rawData                         
setNames(cells.rawData.list, paste0('c',i,'.cells.rawData'))}

# sums by gene per celltype 
sums <- c(NA)
sums.list <- list()
for(i in celltype){sums <- rowSums(cells.rawData.list[[i]] != 0)
# assign(paste0('c',i,'.sums'), sums)
sums.list[[i]] <- sums                         
setNames(sums.list, paste0('c',i,'.sums'))}

rm(rawData)
rm(cells, sums, cells.rawData)
rm(i)

# pct by gene per celltype
pct <- c(NA)
pct.list <- list()
for(i in celltype){pct <- sums.list[[i]]/length(cells.list[[i]])
# assign(paste0('c',i,'.pct'), pct)
pct.list[[i]] <- pct                         
setNames(pct.list, paste0('c',i,'.pct'))}

pct <- data.frame(sapply(pct.list,c))

rm(i)
rm(pct.list)

for(i in 1:length(celltype))
{names(pct)[i] <- celltype[i]}

saveRDS(pct, paste0(file_prefix, ".celltype.pct.rds"))
write.table(pct, paste0(file_prefix, ".celltype.pct.txt"), sep="\t", col.names=NA)

return(pct)
}

GenerateSigDEGsCountBarplot <- function(DEGs_file , commGenes, pct.min, pValadj.max, log2FC.min,
																				celltype, clust.colors, file_prefix){

	#### 计算的是符合pvaladj 小于等于0.05 avgLog2FC >= 0.1 以及各个组中表达细胞数目大于等于10%的基因占总的找出来差异的基因数目的比例
	### 总的基因数目是去除掉mt- Rpl Rps后的基因数目
	
	celltyeSigDEGs <- list()
	pctSig <- c()
	celltype_deg <- c()
	for (i in names(DEGs_file)){
		
		print(i)
		## remove mt and rb gene
		index <- grep("^mt-|^Rpl|^Rps", DEGs_file[[i]]$gene)
		temp <- DEGs_file[[i]][-index, ]
		
		index_a <- temp$p_val_adj <= pValadj.max & abs(temp$avg_log2FC) >= log2FC.min
		index_b <- temp$pct.1 >= pct.min & temp$pct.2 >= pct.min 
		index <- index_a & index_b
		if(sum(index) > 0){
			temp <- temp[index, ]
			temp <- mutate(temp, updown = ifelse(avg_log2FC>=0, "up", "down"))
			print(table(temp$updown))
			pctSig <- rbind(pctSig, table(temp$updown))
			celltyeSigDEGs[[i]] <- temp
			celltype_deg <- c(celltype_deg, i)
		} else{
			next
		}
		}
	
	rownames(pctSig) <- celltype_deg

	print(celltype_deg)
	pctSig <- as.data.frame(pctSig)
	pctSig$celltype <- factor(rownames(pctSig), levels = celltype[celltype%in%celltype_deg])
	pctSig$group <- file_prefix
	
	plot_data<- reshape2::melt(data = pctSig, id.vars = c("celltype", "group"), measure.vars = c("down", "up"), 
											 value.name = "count")
	plot_data$pct <- (plot_data$count/length(commGenes)) * 100
	plot_data$pct[plot_data$variable=="down"] <- -(plot_data$pct[plot_data$variable=="down"] )
	theme1 <- theme_bw()  + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),# 移除bw主题的网格线
																plot.title = element_text(hjust = 0.5)) # 主题居中
	
	max_y <- ceiling(max(abs(plot_data$pct)))
	p1 <- ggplot(plot_data, aes(x=celltype, y=pct, fill = factor(celltype)))
	p1 <- p1 + geom_bar(stat="identity", color="black", width=0.8) 
	p1 <- p1+ scale_fill_manual(values=as.character(clust.colors[celltype%in%celltype_deg]))
	p1 <- p1 + theme1
	p1 <- p1 + theme(legend.position="none") 
	p1 <- p1 + labs(title=paste0(file_prefix, "_pValadj.max_", pValadj.max, "_log2FC.min_", log2FC.min, "_pct.min_", pct.min))
	
	p1 <- p1 + theme(axis.text.x = element_text(size = 8, face = "bold", vjust = 0.1, hjust = 0.1,
																							angle = 30))
	p1 <- p1 + xlab("") + ylab("% Gene")
	p1 <- p1 + scale_y_continuous(limits=c(-max_y, max_y), breaks=seq(-max_y, max_y, 1))
	p1 <- p1 + coord_flip()
	
	ggsave(paste0(file_prefix, "_pValadj.max_", pValadj.max, "_log2FC.min_", log2FC.min, "_pct.min_", pct.min, "_CountPlot.pdf"), 
				 plot  = p1,  width = 14, height = 7)
	res <- list(celltyeSigDEGs = celltyeSigDEGs, pctSig = pctSig, plot = p1)
	saveRDS(res, paste0(file_prefix, "_pctSig_count.rds"))
	return(res)
}

generateStripPlotData <- function(deg_rds, mt_rb_pattern, padj_cutoff,
																	fc_cutoff, top_gene_number, col_used){
	
	# colnames(deg_rds) is "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cell_type"
	
	plot_data <- NULL
	for (i in seq_along(names(deg_rds))){
		print(i)
		print(names(deg_rds)[i])
		temp <- deg_rds[[i]]
		dim(temp)
		
		index <- grep(mt_rb_pattern, temp$gene)
		if(length(index) > 0){temp <- temp[-index, ]}else{temp <- temp}
		dim(temp)
		temp <- mutate(temp, updown = ifelse(p_val_adj <= padj_cutoff & avg_log2FC >= fc_cutoff, "Upregulated",
																				 ifelse(p_val_adj <= padj_cutoff & avg_log2FC <= -fc_cutoff, "Downregulated", "NS")))
		
		sigs <- dplyr::filter(temp, updown != "NS") 
		top_n_gene <- sigs %>% slice_max(order_by = avg_log2FC, n = top_gene_number)
		bot_n_gene <- sigs %>% slice_min(order_by = avg_log2FC, n = top_gene_number)
		top_gene <- dplyr::full_join(top_n_gene, bot_n_gene)
		temp$updown[temp$gene %in% top_gene$gene] <- "top_label"
		
		#set color to use
		temp <- mutate(temp, col_use = ifelse(updown== "Upregulated" | updown== "Downregulated",
																					col_used[i], "dark gray"))
		temp$col_use[temp$updown == "top_label"] <- "red"
		plot_data <- rbind(plot_data, temp)
	}
	return(plot_data)
}

drawStripPlot <- function(plot_data, plot_title_name, file_prefix, file_width, file_height){
	library(ggplot2)
	theme_strip <- theme_bw()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
																		 legend.position = "none", axis.line = element_line(colour = "black"),
																		 panel.border = element_blank()) 
	theme_strip <- theme_strip + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
	theme_strip <- theme_strip + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
																		 axis.text.x = element_text(size = 10, angle = 90, hjust=1))
	
	strip <- ggplot(plot_data, aes( x = cell_type, y = avg_log2FC)) 
	strip <- strip + ggrastr::geom_jitter_rast(data=plot_data[plot_data$updown == "NS",], 
																						 color="dark grey", width = 0.2, height = 0.0, 
																						 alpha = .25, shape = 1, raster.dpi = 300)
	strip <- strip + ggrastr::geom_jitter_rast(data = plot_data[plot_data$updown != "NS", ],
																						 color = plot_data[plot_data$updown != "NS", ]$col_use, 
																						 width = 0.2, height = 0.0, alpha = 1, shape = 1, raster.dpi = 300, size = 0.8)
	strip <- strip + ggrepel::geom_text_repel(data = plot_data[plot_data$updown == "top_label", ], 
																						label = plot_data[plot_data$updown == "top_label", ]$gene, 
																						fontface = "italic", size = 2.5, 
																						max.overlaps = Inf) 
	strip <- strip + theme_strip + ggtitle(plot_title_name)
	save_file_name <- paste0(file_prefix,"strip_plot.pdf")
	ggsave(filename = save_file_name, strip, width = file_width, height = file_height)
	return(strip)}

IdentifyCelltypeByTop2Markers <- function(topn_markers = top, cellmarker = cellmarker){
	
	predict_celltype <- data.frame(cluster = unique(top$cluster))
	predict_celltype$top2gene <- ""
	predict_celltype$celltype <- ""
	
	for (row_i in  predict_celltype$cluster){
		
		cat("Cluster", row_i, "\n")
		search_gene <- paste0(paste0("^", top$gene[top$cluster==row_i]), collapse = "|")
		print(search_gene)
		top_celltype <- unique(cellmarker$cell[grep(search_gene, cellmarker$gene, ignore.case = T)])
		predict_celltype$celltype[row_i] <- paste(top_celltype, collapse = ";")
		predict_celltype$top2gene[row_i] <- paste(search_gene, collapse = ";")
	}
	return(predict_celltype)
}


library(purrr)
GatherData <- function(object, ...) {
	
	UseMethod("GatherData")
	
}


GatherData.Seurat <- function(object,
															assay,
															slot_use,
															...) {
	
	assay <- assay %||% "RNA"
	slot_use <- slot_use %||% "data"
	obj_data <- GetAssayData(
		object = object,
		assay = assay,
		slot = slot_use
	) %>%
		as.matrix()
	return(obj_data)
}


PseudoCell  <- function(object,
												assay_use = NULL,
												slot_use = NULL,
												cluster_use =NULL,
												pseudocell.size  =NULL){
	message("tips: 
  Cluster_use : one col in metadata
  pseudocell.size : how many cell will be pseudo")
	
	Inter<- GatherData(object = object,
										 assay = assay_use,
										 slot_use = slot_use) 
	Inter[Inter<0]=0
	idd<-object@meta.data
	Inter.id<-cbind(rownames(idd),as.vector(idd[,cluster_use]))
	
	rownames(Inter.id)<-rownames(idd)
	colnames(Inter.id)<-c("CellID","Celltype")
	
	Inter.id<-as.data.frame(Inter.id)
	Inter1<-Inter[,Inter.id$CellID]
	Inter<-as.matrix(Inter1)
	pseudocell.size = pseudocell.size ## 10 test
	new_ids_list = list()
	Inter.id$Celltype <- as.factor(Inter.id$Celltype)
	message("generate new_ids_list")
	for (i in 1:length(levels(Inter.id$Celltype))) {
		cluster_id = levels(Inter.id$Celltype)[i]
		cluster_cells <- rownames(Inter.id[Inter.id$Celltype == cluster_id,])
		cluster_size <- length(cluster_cells)       
		pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
		pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
		names(pseudo_ids) <- sample(cluster_cells)  
		new_ids_list[[i]] <- pseudo_ids     
	}
	
	new_ids <- unlist(new_ids_list)
	new_ids <- as.data.frame(new_ids)
	new_ids_length <- table(new_ids)
	
	new_colnames <- rownames(new_ids)  ###add
	all.data<-Inter[,as.character(new_colnames)] ###add
	all.data <- t(all.data)###add
	
	message("aggregate new_ids_list gene expression through mean method")
	new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
											list(name=new_ids[,1]),FUN=mean)
	rownames(new.data)<-new.data$name
	new.data<-new.data[,-1]
	
	new_ids_length<-as.matrix(new_ids_length) ##
	short <- which(new_ids_length < pseudocell.size -1 ) ## 去掉小于9的细胞数目
	new_good_ids<-as.matrix(new_ids_length[-short,]) ##
	result<-t(new.data)[,rownames(new_good_ids)]
	rownames(result)<-rownames(Inter)
	
	message("create seuratobject and save results")
	newobject <- CreateSeuratObject(result)
	newobject@misc$idtrans <- new_ids
	newobject@commands$PseudoCell <- LogSeuratCommand(newobject, return.command = TRUE)
	return(newobject)
}


## celltype similarity
AvgExpSpecies <- function(species1_marker, species2_marker, species1_object, species2_object){
	
	intersect_markers <- intersect(species1_marker$gene, species2_marker$gene) #568
	species1_avgexp <- AverageExpression(species1_object, group.by = "ident", features = intersect_markers, assays = "RNA")
	species2_avgexp <- AverageExpression(species2_object, group.by = "ident", features = intersect_markers, assays = "RNA")
	
	Sp1 = as.data.frame(species1_avgexp$RNA[rowSums(as.matrix(species1_avgexp$RNA))!=0,])
	colnames(Sp1) <- paste0("Sp1_", colnames(Sp1) )
	Sp2 = as.data.frame(species2_avgexp$RNA[rowSums(as.matrix(species2_avgexp$RNA))!=0,])
	colnames(Sp2) <- paste0("Sp2_", colnames(Sp2) )
	
	nDESp1 = length(species1_marker$gene)
	nDESp2 = length(species2_marker$gene)
	avg_list <- list(species1_avgexp = Sp1,
									 species2_avgexp = Sp2,
									 nDESp1 = nDESp1,
									 nDESp2 = nDESp2,
									 intersect_markers = intersect_markers)
	return(avg_list)
}

CorrCompareFunction <- function(Sp1_exp, Sp2_exp, intersect_markers, nDESp1, nDESp2){
	
	#Step5: Scale Expression Tables by gene average
	avg = rowMeans(Sp1)
	Sp1 = sweep(Sp1,1,avg,"/")
	rm(avg)
	
	avg = rowMeans(Sp2)
	Sp2 = sweep(Sp2,1,avg,"/")
	rm(avg)
	
	#Step6: Merge Expression Tables
	geTable = merge(Sp1,Sp2, by='row.names', all=F)
	rownames(geTable) = geTable$Row.names
	geTable = geTable[,2:ncol(geTable)]
	
	#Step7:  Correlation
	#7a:  Correlation
	corr.method <- "spearman"
	Corr.Coeff.Table = cor(geTable,method=corr.method)
	
	#7b:  Shuffle data
	shuffled.cor.list = list()
	pb   <- txtProgressBar(1, 100, style=3)
	nPermutations <- 1000
	for (i in 1:nPermutations){
		shuffled = apply(geTable[,1:ncol(Sp1)],1,sample)
		shuffled2 = apply(geTable[,(ncol(Sp1)+1):ncol(geTable)],1,sample)
		shuffled = cbind(t(shuffled),t(shuffled2))
		shuffled.cor = cor(shuffled,method=corr.method)
		shuffled.cor.list[[i]] = shuffled.cor
		rm(list=c('shuffled','shuffled2','shuffled.cor'))
		if ((i %% 100) ==0){
			setTxtProgressBar(pb, (i*100)/nPermutations)
		}
	}
	
	p.value.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
	rownames(p.value.table) = colnames(geTable)
	colnames(p.value.table) = colnames(geTable)
	
	shuffled.mean.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
	rownames(shuffled.mean.table) = colnames(geTable)
	colnames(shuffled.mean.table) = colnames(geTable)
	
	a = combn(1:ncol(geTable),2)
	for (i in 1:ncol(a)){
		cor.scores = sapply(shuffled.cor.list,"[",a[1,i],a[2,i])
		shuffled.mean.table[a[1,i],a[2,i]] = mean(cor.scores)
		shuffled.mean.table[a[2,i],a[1,i]] = mean(cor.scores)
		p.value = mean(abs(cor.scores)>=abs(Corr.Coeff.Table[a[1,i],a[2,i]]))
		p.value.table[a[1,i],a[2,i]] = p.value
		p.value.table[a[2,i],a[1,i]] = p.value
		rm(list=c('cor.scores','p.value'))
		setTxtProgressBar(pb, (i*100)/ncol(a))
	}
	
	p.value.table[p.value.table==0] <- 0.000001
	neg.log10.p = -log10(p.value.table)
	
	#step8 "Overlap in Markers"
	#for all pairs of cell-types generate list of genes that are at least 1.5x avg in both cells
	
	#from above a = combn(1:ncol(geTable),2)
	marker.overlap.list = list()
	for (i in 1:ncol(a)){
		datasubset = cbind(geTable[,a[1,i]],geTable[,a[2,i]])
		markers = rownames(geTable[datasubset[,1]>1 & datasubset[,2]>1,])
		marker.overlap.list[[i]] = markers
		names(marker.overlap.list)[i] = paste(colnames(geTable)[a[1,i]], colnames(geTable)[a[2,i]],sep='_')
		rm(list=c('datasubset','markers'))
	}
	
	
	allgene.list.to.return = list(Corr.Coeff.Table,shuffled.mean.table,
																p.value.table,neg.log10.p,
																intersect_markers,
																rownames(geTable),
																length(intersect_markers),
																length(rownames(geTable)),nDESp1, nDESp2, geTable,marker.overlap.list)
	names(allgene.list.to.return) = c('corr.coeff','shuffled_correlation_score_means',
																		'p.value', 
																		'negative_log10_p.value',
																		'DEgenes_intersect',
																		'DEgenes_in_analysis',
																		'nDEgenes_intersect',
																		'nDEgenes_in_analysis',
																		'nDEgenes_Sp1','nDEgenes_Sp2',
																		'scaled_table','overlapping_markers')
	return(allgene.list.to.return)
}

cal_scLink <- function(count.norm){
	library(openxlsx)
	library(scLink)
	library(ggplot2)
	library(ggplotify)
	library(cowplot)
	library(dplyr)
	library(ggraph)
	library(igraph)
	library(RColorBrewer)
	
	networks = sclink_net(expr = count.norm, ncores = 1, lda = c(0.5, 0.2, 0.1, 0.05, 0.01))
	corr = sclink_cor(expr = count.norm, ncores = 1)
	res_wt <- c(networks = networks, corr = corr)
	
	net <- networks$summary[[1]]
	aaa <- data.frame(gene1="a", gene2="b", adj=1, Sigma=1, corr=1, abs_corr = 1)
	genenames <- rownames(net$adj)
	nrow(net$adj)
	for (i in 1:(nrow(net$adj)-1)){
		# print(i)
		nextcol <- i + 1
		
		for (j in nextcol:nrow(net$adj)){
			#print(j)
			temp <- data.frame(gene1=genenames[i], gene2=genenames[j], 
												 adj=net$adj[i,j], Sigma=net$Sigma[i,j], corr=networks$cor[i,j])
			temp$abs_corr <- abs(temp$corr)
			aaa <- rbind(aaa, temp)
		} 
	}
	aaa <- aaa[-1, ]
	return(aaa)
}


plot_network <- function(res, color_genename){
	
	
	info <- data.frame(from = res$gene1, to = res$gene2, 
										 combined_score = res$abs_corr)
	
	# 转换stringID为Symbol，只取前两列和最后一列
	links <- info %>%  
		dplyr::select(from, to , last_col()) %>% 
		dplyr::rename(weight = combined_score)
	
	# 去除游离的互作关系
	# 如果links数据框的一个link的from只出现过一次，同时to也只出现一次，则将其去除
	links_2 <- links %>% mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
		mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
		filter(!(from_c == 1 & to_c == 1)) %>%
		dplyr::select(1,2,3)
	# 新的节点数据
	nodes_2 <- links_2 %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
	# 创建网络图
	net_2 <- igraph::graph_from_data_frame(d=links_2,vertices=nodes_2,directed = F)
	# 添加一些参数信息用于后续绘图
	# V和E是igraph包的函数，分别用于修改网络图的节点（nodes）和连线(links)
	# 添加必要的参数
	igraph::V(net_2)$deg <- igraph::degree(net_2)
	igraph::V(net_2)$size <- igraph::degree(net_2)/5
	igraph::E(net_2)$width <- igraph::E(net_2)$weight
	
	p_net <- ggraph(net_2,layout = "linear", circular = TRUE)+
		geom_edge_arc(aes(edge_width=width),  color="lightgrey", show.legend = T)+
		geom_edge_arc(aes(filter=width>0.5, edge_colour="red", edge_width=width)) +
		geom_node_point(aes(size=size,color=factor(group)),  alpha=0.5)+
		#geom_node_text(aes(filter=deg>5, label=name), size = 5, repel = F, family = "Arial")+
		geom_node_text(aes(label=name), size = 5, repel = F, family = "Arial") +
		scale_edge_width(range = c(0.001,0.5), breaks = c(0.2, 0.5, 0.7), limits = c(0.2, 1)) + 
		scale_size_continuous(range = c(1,max(igraph::degree(net_2))) )+
		guides(size=F)+
		theme_graph()
	p_net$data$group <- p_net$data$name %in% color_genename
	p_net <- p_net + geom_node_point(aes(size=size,color=factor(group)),  alpha=0.5)
	return(p_net)
}


hclustBasedTopNumMarkers <- function(seurat_object, file_prefix, min_pct, logfc_threshold, p_val_adj,
																		 mt_rb_pattern, top_num, width_num, height_num, cluster_method){
	
	library(factoextra)
	library(igraph)
	
	Cluster_markers <- FindmarkerForCluster(seurat_object, file_prefix, min.pct = min_pct, logfc.threshold = logfc_threshold,
																					p_val_adj = p_val_adj, mt_rb_pattern = mt_rb_pattern)
	TopMarkersInClusters <- TopMarkersInCluster(Cluster_markers, file_prefix, top_num)
	
	#### hclust
	Variable_gene <- rev(unique(TopMarkersInClusters$gene))
	
	hclust_exp <- t(seurat_object@assays$RNA@data[Variable_gene, ])
	hclust_exp_mean <- aggregate(hclust_exp, 
															 by = list(celltype = as.character(seurat_object@active.ident)), mean)
	
	celltypenames <- hclust_exp_mean$celltype
	hclust_exp_mean <- hclust_exp_mean[, -1]
	rownames(hclust_exp_mean) <- celltypenames
	
	bc.scaled <- scale(hclust_exp_mean)
	d <- dist(bc.scaled)
	fit1 <- hclust(d, method = cluster_method)
	
	p2 <- fviz_dend(fit1, k=4, rect =T, rect_fill = T,horiz = F,
									rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"))
	
	#### top4 marker
	levels(x = seurat_object) <- celltypenames[fit1$order]

	TopMarkersInClusters <- TopMarkersInCluster(Cluster_markers, file_prefix, 4)
	TopMarkersInClusters$cluster <- factor(as.character(TopMarkersInClusters$cluster),
																				 levels = celltypenames[fit1$order])
	TopMarkersInClusters <- arrange(TopMarkersInClusters, cluster)
	markers <- unique(TopMarkersInClusters$gene)
	
	p_dotplot <- DotPlot(seurat_object, features = markers, col.min = 0, col.max = 2)
	
	p_dotplot <- p_dotplot + coord_flip()
	p_dotplot <- p_dotplot + scale_colour_gradient2(low = "lightgrey",high = "Red", midpoint = 1, breaks=c(0, 1, 2), label = c(0, 1, 2))
	p_dotplot <- p_dotplot + theme(axis.text.x = element_text(size = 10, vjust = 1, hjust = 1, angle = 45))
	ggsave(paste0("Celltypemarker_DotPlot_", file_prefix, ".pdf"),  plot = p_dotplot, width = width_num, height = height_num)
	p_dotplot_no <- p_dotplot + NoLegend()
	
	plots <- plot_grid(p2, p_dotplot_no, ncol = 1, rel_heights = c(0.2, 1), align = "v")
	
	pdf(paste0("Cluster_dendrogram_", cluster_method, "_", top_num,"_with_dotplotmarks_", file_prefix,".pdf"), width_num, height_num)
	print(plots)
	dev.off()
}


AvgExpSpecies <- function(species1_marker, species2_marker, species1_object, species2_object){
	
	intersect_markers <- intersect(species1_marker$gene, species2_marker$gene) #568
	species1_avgexp <- AverageExpression(species1_object, group.by = "ident", features = intersect_markers, assays = "RNA")
	species2_avgexp <- AverageExpression(species2_object, group.by = "ident", features = intersect_markers, assays = "RNA")
	
	Sp1 = as.data.frame(species1_avgexp$RNA[rowSums(as.matrix(species1_avgexp$RNA))!=0,])
	colnames(Sp1) <- paste0("Sp1_", colnames(Sp1) )
	Sp2 = as.data.frame(species2_avgexp$RNA[rowSums(as.matrix(species2_avgexp$RNA))!=0,])
	colnames(Sp2) <- paste0("Sp2_", colnames(Sp2) )
	
	nDESp1 = length(species1_marker$gene)
	nDESp2 = length(species2_marker$gene)
	avg_list <- list(species1_avgexp = Sp1,
									 species2_avgexp = Sp2,
									 nDESp1 = nDESp1,
									 nDESp2 = nDESp2,
									 intersect_markers = intersect_markers)
	return(avg_list)
}

CorrCompareFunction <- function(Sp1, Sp2, intersect_markers, nDESp1, nDESp2){
	
	#Step5: Scale Expression Tables by gene average
	avg = rowMeans(Sp1)
	Sp1 = sweep(Sp1,1,avg,"/")
	rm(avg)
	
	avg = rowMeans(Sp2)
	Sp2 = sweep(Sp2,1,avg,"/")
	rm(avg)
	
	#Step6: Merge Expression Tables
	geTable = merge(Sp1,Sp2, by='row.names', all=F)
	rownames(geTable) = geTable$Row.names
	geTable = geTable[,2:ncol(geTable)]
	
	#Step7:  Correlation
	#7a:  Correlation
	corr.method <- "spearman"
	Corr.Coeff.Table = cor(geTable, method=corr.method)
	
	#7b:  Shuffle data
	shuffled.cor.list = list()
	pb   <- txtProgressBar(1, 100, style=3)
	nPermutations <- 1000
	for (i in 1:nPermutations){
		shuffled = apply(geTable[,1:ncol(Sp1)],1,sample)
		shuffled2 = apply(geTable[,(ncol(Sp1)+1):ncol(geTable)],1,sample)
		shuffled = cbind(t(shuffled),t(shuffled2))
		shuffled.cor = cor(shuffled, method=corr.method)
		shuffled.cor.list[[i]] = shuffled.cor
		rm(list=c('shuffled','shuffled2','shuffled.cor'))
		if ((i %% 100) ==0){
			setTxtProgressBar(pb, (i*100)/nPermutations)
		}
	}
	
	p.value.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
	rownames(p.value.table) = colnames(geTable)
	colnames(p.value.table) = colnames(geTable)
	
	shuffled.mean.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
	rownames(shuffled.mean.table) = colnames(geTable)
	colnames(shuffled.mean.table) = colnames(geTable)
	
	a = combn(1:ncol(geTable),2)
	for (i in 1:ncol(a)){
		cor.scores = sapply(shuffled.cor.list,"[",a[1,i],a[2,i])
		shuffled.mean.table[a[1,i],a[2,i]] = mean(cor.scores)
		shuffled.mean.table[a[2,i],a[1,i]] = mean(cor.scores)
		p.value = mean(abs(cor.scores)>=abs(Corr.Coeff.Table[a[1,i],a[2,i]]))
		p.value.table[a[1,i],a[2,i]] = p.value
		p.value.table[a[2,i],a[1,i]] = p.value
		rm(list=c('cor.scores','p.value'))
		setTxtProgressBar(pb, (i*100)/ncol(a))
	}
	
	p.value.table[p.value.table==0] <- 0.000001
	neg.log10.p = -log10(p.value.table)
	
	#step8 "Overlap in Markers"
	#for all pairs of cell-types generate list of genes that are at least 1.5x avg in both cells
	
	#from above a = combn(1:ncol(geTable),2)
	marker.overlap.list = list()
	for (i in 1:ncol(a)){
		datasubset = cbind(geTable[,a[1,i]],geTable[,a[2,i]])
		markers = rownames(geTable[datasubset[,1]>1 & datasubset[,2]>1,])
		marker.overlap.list[[i]] = markers
		names(marker.overlap.list)[i] = paste(colnames(geTable)[a[1,i]], colnames(geTable)[a[2,i]],sep='_')
		rm(list=c('datasubset','markers'))
	}
	
	
	allgene.list.to.return = list(Corr.Coeff.Table,shuffled.mean.table,
																p.value.table,neg.log10.p,
																intersect_markers,
																rownames(geTable),
																length(intersect_markers),
																length(rownames(geTable)),nDESp1, nDESp2, geTable,marker.overlap.list)
	names(allgene.list.to.return) = c('corr.coeff','shuffled_correlation_score_means',
																		'p.value', 
																		'negative_log10_p.value',
																		'DEgenes_intersect',
																		'DEgenes_in_analysis',
																		'nDEgenes_intersect',
																		'nDEgenes_in_analysis',
																		'nDEgenes_Sp1','nDEgenes_Sp2',
																		'scaled_table','overlapping_markers')
	return(allgene.list.to.return)
}

Hfunc <- function(v) {
	v <- v[v > 0]
	return(-sum(v * log(v)))
}

CalculateTFSS <- function(TF_Symbol, seurat_object, min_pct, logfc_threshold, p_val_adj){
	
	### intersect with seurat_object
	TF_list <- intersect(TF_Symbol, rownames(seurat_object))
	
	### narrow down the TF list
	DefaultAssay(seurat_object) <- "RNA"
	seurat_object <- NormalizeData(seurat_object)
	TF_markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_threshold, features = TF_list)
	TF_markers <- TF_markers[TF_markers$p_val_adj < p_val_adj, ]
	
	cat("The length of the TF_markers:")
	print(nrow(TF_markers))
	### generate p_E_table
	p_E_table <- seurat_object@assays$RNA@data[unique(TF_markers$gene), ]
	p_E_table <- p_E_table[rowSums(p_E_table > 0) > 10 , ]
	p_E_table <- as.matrix(p_E_table)
	## normalize the data so that the sum of the vector is 1
	p_E_table <- sweep(p_E_table, 1, rowSums(p_E_table), FUN = "/")
	cat("The row and column of p_E_table:")
	print(dim(p_E_table))
	
	### generate cell_label
	
	cell_label <- matrix(0, length(unique(seurat_object@active.ident)), length(seurat_object@active.ident))
	rownames(cell_label) <- levels(seurat_object@active.ident)
	colnames(cell_label) <- colnames(seurat_object)
	for(i in rownames(cell_label)){
		index <- as.character(seurat_object@active.ident) == i
		cell_label[i, index] <- 1
	}
	## normalize the data so that the sum of the vector is 1
	cell_label <- sweep(cell_label, 1, rowSums(cell_label), FUN = "/")
	cat("The row and column of cell_label:")
	print(dim(cell_label))
	
	### calculate Trancription Factor Specific score
	TF_num <- nrow(p_E_table)
	celltype_num <- nrow(cell_label)
	TFSS <- matrix(0, TF_num, celltype_num)
	rownames(TFSS) <- rownames(p_E_table)
	colnames(TFSS) <- rownames(cell_label)
	
	for(TF_i in 1:TF_num){
		print(TF_i)
		p_E <- p_E_table[TF_i, ]
		for(celltype_i in 1:celltype_num){
			print(celltype_i)
			p_C <- cell_label[celltype_i, ]
			JSD <- Hfunc(0.5 * (p_E+p_C)) -  0.5 * (Hfunc(p_E) + Hfunc(p_C))
			TFSS[TF_i, celltype_i] <- 1-sqrt(JSD)
		}
	}
	res <- list(TFSS = TFSS, TF_markers = TF_markers)
	return(res)
}


generate_sclink_network <- function(count, used_genes){
	
	library(scLink)
	library(reshape)
	# count  行是细胞，列是基因
	count.norm = sclink_norm(count, scale.factor = 1e6, filter.genes = F, gene.names = used_genes)
	
	## filer cells and genes
	row_sum <- apply(count.norm > 0, 1, sum)
	col_sum <- apply(count.norm > 0, 2, sum)
	count.norm <- count.norm[row_sum > 2, col_sum > 2]
	## network
	lda_value <- seq(1, 0.1, -0.05)
	networks = sclink_net(expr = count.norm, ncores = 4, lda = lda_value)
	
	bic_value <- c()
	for (i in 1:length(lda_value)){ bic_value <- c(bic_value, networks$summary[[i]]$bic)}
	
	plot_bic <- data.frame(lda = lda_value, bic = bic_value)
	print(ggplot(plot_bic) +
		geom_point(mapping = aes(x = lda_value, y = bic_value)))
	
	mat <- networks$summary[[which.min(bic_value)]]$adj
	mat[!upper.tri(mat, diag = TRUE)] <- 0
	net_adj <- reshape::melt(mat)
	net_adj <- net_adj[net_adj$value==1, ]
	net_adj <- net_adj[as.character(net_adj$X1)!=as.character(net_adj$X2), ]
	networks$summary[[which.min(bic_value)]]$net_adj <- net_adj
	networks$summary[[which.min(bic_value)]]$corr <- networks$cor
	networks$summary[[which.min(bic_value)]]$count.norm <- count.norm
	return(networks$summary[[which.min(bic_value)]])
}


get_corr_p_permutation <- function(geTable, nPermutations = 500){
	library(scLink)
	#step1:  original correlation
	# geTable: row sample, col gene
	# scLink corr Correlation
	Corr.Coeff.Table = sclink_cor(geTable, ncores = 1)
	
	#step2:  Shuffle data
	shuffled.cor.list = list()
	pb   <- txtProgressBar(1, 100, style=3)
	for (i in 1:nPermutations){
		shuffled = apply(geTable, 2,sample)
		shuffled.cor = sclink_cor(shuffled, ncores = 1)
		shuffled.cor.list[[i]] = shuffled.cor
		rm(list=c('shuffled','shuffled.cor'))
		if ((i %% 100) ==0){
			setTxtProgressBar(pb, (i*100)/nPermutations)
		}
	}
	
	# step3: pvalue 
	p.value.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
	rownames(p.value.table) = colnames(geTable)
	colnames(p.value.table) = colnames(geTable)
	
	shuffled.mean.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
	rownames(shuffled.mean.table) = colnames(geTable)
	colnames(shuffled.mean.table) = colnames(geTable)
	
	a = combn(1:ncol(geTable),2)
	for (i in 1:ncol(a)){
		cor.scores = sapply(shuffled.cor.list,"[",a[1,i],a[2,i])
		shuffled.mean.table[a[1,i],a[2,i]] = mean(cor.scores)
		shuffled.mean.table[a[2,i],a[1,i]] = mean(cor.scores)
		p.value = mean(abs(cor.scores)>=abs(Corr.Coeff.Table[a[1,i],a[2,i]]))
		p.value.table[a[1,i],a[2,i]] = p.value
		p.value.table[a[2,i],a[1,i]] = p.value
		rm(list=c('cor.scores','p.value'))
		setTxtProgressBar(pb, (i*100)/ncol(a))
	}
	
	p.value.table[p.value.table==0] <- 0.000001
	return(p.value.table)
}


single_parallel <- function(func,iterable,...){
	"
  :param func:被并行函数
  :param iteralbe:func的1个动态参数(vector、list)
  :param ...:func的静态参数
  :return list,与iterable等长
  "
	#1.加载包
	library(parallel)
	#2.计算计算机内核数
	cores <- detectCores(logical = FALSE)
	#3.打开并行计算
	cl <- makeCluster(cores)
	#4.给每个单独内核传递变量，函数等
	clusterExport(cl,deparse(substitute(func)))
	#5.开始并行计算（用法与sapply类似）
	result <- parSapply(cl,iterable,func,...)
	#6.关闭并行计算
	stopCluster(cl)
	return(result)
}


select_GENIE3_Links <- function(linkList, corrMat, pMat, weightThreshold = 0.001, topThr = 0.01, 
																nTopTfs = c(5, 10, 20), nTopTargets = 20, aFun = topPerTf, weightCol = "weight", verbose = TRUE) 
{

	if (!all(c("TF", "Target", weightCol) %in% colnames(linkList))) 
		stop("The link list colnames should be \"TF\", \"Target\", \"", 
				 weightCol, "\"")
	cntPairs <- table(table(linkList[, "TF"], linkList[, "Target"]))
	if (any(names(cntPairs) > 1)) 
		stop("There are duplicated regulator-target (gene id/name) pairs in the input link list.")
	msg <- paste0(format(Sys.time(), "%H:%M"), "\tCreating TF modules")

	linkList <- linkList[which(linkList[, weightCol] >= weightThreshold), ]
	if (verbose) {
		message("Number of links between TFs and targets (weight>=", 
						weightThreshold, "): ", nrow(linkList))
	}
	tfModules <- list()
	linkList$TF <- as.character(linkList$TF)
	linkList$Target <- as.character(linkList$Target)
	allName <- paste0("w", format(weightThreshold, scientific = FALSE))
	tfModules[[allName]] <- split(linkList$Target, factor(linkList$TF))
	if (!is.null(topThr)) {
		topThr <- setNames(topThr, paste("w", format(topThr, 
																								 scientific = FALSE), sep = ""))
		for (i in seq_along(topThr)) {
			llminW <- linkList[which(linkList[, weightCol] > 
															 	topThr[i]), ]
			tfModules[[names(topThr)[i]]] <- split(llminW$Target, 
																						 factor(llminW$TF))
		}
	}
	if (!is.null(nTopTargets)) {
		nTopTargets <- setNames(nTopTargets, paste("top", nTopTargets, 
																							 sep = ""))
		for (i in seq_along(nTopTargets)) {
			tfModules[[names(nTopTargets)[i]]] <- lapply(tfModules[[allName]], 
																									 function(x) x[1:(min(length(x), nTopTargets[i]))])
		}
	}
	tfModules.melted <- reshape2::melt(tfModules)
	colnames(tfModules.melted) <- c("Target", "TF", "method")
	topTFsperTarget.asDf <- NULL
	if (!is.null(nTopTfs)) {
		linkList_byTarget <- split(linkList, factor(linkList$Target))
		nTopTfs <- setNames(nTopTfs, paste("top", nTopTfs, "perTarget", 
																			 sep = ""))
		topTFsperTarget <- lapply(linkList_byTarget, function(llt) {
			nTFs <- nTopTfs[which(nTopTfs <= nrow(llt))]
			reshape2::melt(lapply(nTFs, function(x) llt[1:x, 
																									"TF"]))
		})
		topTFsperTarget <- topTFsperTarget[which(!sapply(sapply(topTFsperTarget, 
																														nrow), is.null))]
		topTFsperTarget.asDf <- data.frame(data.table::rbindlist(topTFsperTarget, 
																														 idcol = TRUE))
		colnames(topTFsperTarget.asDf) <- c("Target", "TF", "method")
	}
	byFun <- NULL
	if (!is.null(aFun)) {
		byFun <- aFun(linkList, weightCol)
	}
	tfModules <- rbind(tfModules.melted, topTFsperTarget.asDf, 
										 byFun)
	rm(tfModules.melted)
	rm(topTFsperTarget.asDf)
	tfModules$TF <- as.character(tfModules$TF)
	tfModules$Target <- as.character(tfModules$Target)
	if (verbose) 
		print(rbind(nTFs = length(unique(tfModules$TF)), nTargets = length(unique(tfModules$Target)), 
								nGeneSets = nrow(unique(tfModules[, c("TF", "method")])), 
								nLinks = nrow(tfModules)))

	if (!is.null(corrMat) & !is.null(pMat)) {
		tfs <- unique(tfModules$TF)
		missingTFs <- tfs[which(!tfs %in% rownames(corrMat))]
		if (length(missingTFs) > 0) {
			warning("The following TFs are missing from the correlation matrix: ", 
							paste(missingTFs, collapse = ", "))
			tfs <- tfs[which(tfs %in% rownames(corrMat))]
			corrMat <- corrMat[tfs, ]
		}
		tfModules_byTF <- split(tfModules, as.factor(tfModules$TF))
		tfModules_withCorr_byTF <- lapply(tfModules_byTF[tfs], 
																			function(tfGeneSets) {
																				tf <- as.character(unique(tfGeneSets$TF))
																				targets <- as.character(tfGeneSets$Target)
																				cbind(tfGeneSets, corr = as.numeric(corrMat[tf, targets]),
																							pvalue = as.numeric(pMat[tf, targets]))
																			})
		tfModules_withCorr_byTF <- tfModules_withCorr_byTF[which(lengths(tfModules_withCorr_byTF) > 
																														 	0)]
		tfModules_withCorr <- data.frame(data.table::rbindlist(tfModules_withCorr_byTF))
		if (length(missingTFs) > 0) {
			tfModules_withCorr <- rbind(tfModules_withCorr, data.frame(tfModules[tfModules$TF %in% 
																																					 	missingTFs, ], corr = NA))
		}
	} else {
		tfModules_withCorr <- data.frame(tfModules, corr = NA)
		if (verbose) 
			message("Correlation information not available. It will not be added to the modules.")
	}
	return(tfModules_withCorr)
}

runCorr <- function(exprMatr, method = "spearman") 
{
	res_corr <- psych::corr.test(exprMatr, method = method, ci = F)
	return(res_corr)
}

topPerTf <- function(ll, weightCol) 
{
	ll <- split(ll, factor(as.character(ll$TF)))
	tptf <- setNames(lapply(names(ll), function(tf) {
		tfMean <- mean(ll[[tf]][, weightCol])
		tfSd <- sd(ll[[tf]][, weightCol])
		list(top1sd = ll[[tf]][which(ll[[tf]][, weightCol] >= 
																 	tfMean + tfSd), "Target"], top3sd = ll[[tf]][which(ll[[tf]][, 
																 																															weightCol] >= tfMean + (3 * tfSd)), "Target"])
	}), names(ll))
	tptf <- reshape2::melt(tptf)
	colnames(tptf)[which(colnames(tptf) == "value")] <- "Target"
	colnames(tptf)[which(colnames(tptf) == "L1")] <- "TF"
	colnames(tptf)[which(colnames(tptf) == "L2")] <- "method"
	tptf <- tptf[, c("Target", "TF", "method")]
	return(tptf)
}

# Hclust function
perform_hclust <- function(exp_count_seurat, top_marker, cluster_method = "ward.D2") {
	variable_gene <- rev(unique(top_marker$gene))
	
	hclust_exp <- exp_count_seurat@assays$RNA@data[variable_gene, ]
	hclust_exp <- t(as.matrix(hclust_exp))
	hclust_exp_mean <- aggregate(hclust_exp,
															 by = list(celltype = as.character(exp_count_seurat@active.ident)), mean)
	
	celltypenames <- hclust_exp_mean$celltype
	hclust_exp_mean <- hclust_exp_mean[, -1]
	rownames(hclust_exp_mean) <- celltypenames
	
	hclust_exp_mean_scaled <- scale(hclust_exp_mean)
	dist_matrix <- dist(hclust_exp_mean)
	fit_hclust <- hclust(dist_matrix, method = cluster_method)
	return(list(hclust_exp_mean_scaled = hclust_exp_mean_scaled, fit_hclust = fit_hclust, dist_matrix = dist_matrix))
}

# Compute the WCSS and average silhouette for different values of k
find_optimal_k <- function(hclust_exp_mean_scaled, fit_hclust, dist_matrix) {
	
	max_k <- nrow(hclust_exp_mean_scaled) - 2
	
	k_range <- 2:max_k
	
	# Initialize silhouette and WSS
	avg_sil <- numeric(length(k_range))
	wss <- numeric(length(k_range))
	
	# Calculate silhouette and WSS
	for (i in seq_along(k_range)) {
		k <- k_range[i]
		cluster_assignments <- cutree(fit_hclust, k)
		
		# Silhouette
		silhouette_result <- silhouette(cluster_assignments, dist_matrix)
		avg_sil[i] <- mean(silhouette_result[, 3])
		
		# WSS
		wss[i] <- sum((nrow(hclust_exp_mean_scaled) - 1) * diag(as.matrix(dist(hclust_exp_mean_scaled))) / 2) -
			sum((cluster.stats(dist_matrix, cluster_assignments))$within.cluster.ss)
	}
	
	# Determine the best K values
	best_k_silhouette <- k_range[which.max(avg_sil)]
	best_k_wss <- k_range[which.min(wss[2:length(wss)]) + 1]
	
	
	return(list(k_range = k_range, wcss = wss, best_k_wss = best_k_wss, silhouette_avg = avg_sil, best_k_silhouette = best_k_silhouette))
}

# Visualize the dendrogram with the chosen number of clusters (k)
visualize_dendrogram <- function(fit_hclust, optimal_k, method_name, main_title) {
	p <- fviz_dend(fit_hclust, k = optimal_k, rect = TRUE, rect_fill = TRUE, horiz = FALSE,
								 rect_border = rainbow(optimal_k), main = main_title)
	return(p)
}

# Save and visualize the dotplot
visualize_dotplot <- function(exp_count_seurat, file_prefix, top_marker, fit_hclust, celltypenames) {
	levels(x = exp_count_seurat) <- celltypenames[fit_hclust$order]
	top_marker$cluster <- factor(as.character(top_marker$cluster),
															 levels = celltypenames[fit_hclust$order])
	top_marker <- arrange(top_marker, cluster)
	markers <- unique(top_marker$gene)
	
	p_dotplot <- DotPlot(exp_count_seurat, features = markers)
	p_dotplot <- p_dotplot + coord_flip()
	p_dotplot <- p_dotplot + scale_color_gradient2(low = "lightgrey", high = "red", midpoint = 1,
																								 breaks = c(0, 1, 2), label = c(0, 1, 2))
	
	p_dotplot <- p_dotplot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
																 axis.text.y = element_text(size = 8),
																 legend.position = "right")
	
	return(p_dotplot)
}

# GeomcolWrap function 
GeomcolWrap <- function(plot_data, x_group, y_value, 
												fill_group, fill_cols,
												title_name, x_axis_name, y_axis_name){
	# plot_data: three columns dataframe(x_group, y_value, fill_group)
	# stat_mode: "identity", stack
	
	theme_used <- theme_bw() + theme(panel.grid.major = element_blank(),
																	 panel.grid.minor = element_blank(), 
																	 axis.line = element_line(colour = "black"),
																	 panel.border = element_blank(), 
																	 aspect.ratio = 1)
	theme_used <- theme_used + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, 
																															hjust = 0.5, size = 10),
																	 axis.text.y = element_text(size = 10),
																	 axis.title.x = element_text(size = 10),
																	 axis.title.y = element_text(size = 10),
																	 plot.title = element_text(hjust = 0.5, size = 15))
	
	plot_p <- ggplot(data = plot_data, aes(x = x_group, y = y_value, fill = fill_group)) 
	plot_p <- plot_p + geom_col(position = "fill")  + theme_used
	# change fill colors
	plot_p <- plot_p + scale_fill_manual(values = fill_cols) 
	
	plot_p <- plot_p + scale_y_continuous(expand = c(0,0))
	plot_p <- plot_p + ggtitle(title_name) + xlab(x_axis_name) + ylab(y_axis_name)
	return(plot_p)
}


library(purrr)
GatherData <- function(object, ...) {
	
	UseMethod("GatherData")
	
}


GatherData.Seurat <- function(object,
															assay,
															slot_use,
															...) {
	
	assay <- assay %||% "RNA"
	slot_use <- slot_use %||% "data"
	obj_data <- GetAssayData(
		object = object,
		assay = assay,
		slot = slot_use
	) %>%
		as.matrix()
	return(obj_data)
}

### used function
AvgExpSpecies <- function(species1_marker, species2_marker, species1_object, species2_object, species1_name, species2_name){
	
	intersect_markers <- intersect(species1_marker$gene, species2_marker$gene) #568
	species1_avgexp <- AverageExpression(species1_object, group.by = "ident", features = intersect_markers, assays = "RNA")
	species2_avgexp <- AverageExpression(species2_object, group.by = "ident", features = intersect_markers, assays = "RNA")
	
	Sp1 = as.data.frame(species1_avgexp$RNA[rowSums(as.matrix(species1_avgexp$RNA))!=0,])
	colnames(Sp1) <- paste(species1_name, colnames(Sp1), sep = "_")
	Sp2 = as.data.frame(species2_avgexp$RNA[rowSums(as.matrix(species2_avgexp$RNA))!=0,])
	colnames(Sp2) <- paste(species2_name, colnames(Sp2), sep = "_" )
	
	nDESp1 = length(species1_marker$gene)
	nDESp2 = length(species2_marker$gene)
	avg_list <- list(species1_avgexp = Sp1,
									 species2_avgexp = Sp2,
									 nDESp1 = nDESp1,
									 nDESp2 = nDESp2,
									 intersect_markers = intersect_markers)
	return(avg_list)
}

CorrCompareFunction <- function(Sp1, Sp2, intersect_markers, nDESp1, nDESp2){
	
	#Step5: Scale Expression Tables by gene average
	avg = rowMeans(Sp1)
	Sp1 = sweep(Sp1,1,avg,"/")
	rm(avg)
	
	avg = rowMeans(Sp2)
	Sp2 = sweep(Sp2,1,avg,"/")
	rm(avg)
	
	#Step6: Merge Expression Tables
	geTable = merge(Sp1,Sp2, by='row.names', all=F)
	rownames(geTable) = geTable$Row.names
	geTable = geTable[,2:ncol(geTable)]
	
	#Step7:  Correlation
	#7a:  Correlation
	corr.method <- "spearman"
	Corr.Coeff.Table = cor(geTable, method=corr.method)
	
	#7b:  Shuffle data
	shuffled.cor.list = list()
	pb   <- txtProgressBar(1, 100, style=3)
	nPermutations <- 1000
	for (i in 1:nPermutations){
		shuffled = apply(geTable[,1:ncol(Sp1)],1,sample)
		shuffled2 = apply(geTable[,(ncol(Sp1)+1):ncol(geTable)],1,sample)
		shuffled = cbind(t(shuffled),t(shuffled2))
		shuffled.cor = cor(shuffled, method=corr.method)
		shuffled.cor.list[[i]] = shuffled.cor
		rm(list=c('shuffled','shuffled2','shuffled.cor'))
		if ((i %% 100) ==0){
			setTxtProgressBar(pb, (i*100)/nPermutations)
		}
	}
	
	p.value.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
	rownames(p.value.table) = colnames(geTable)
	colnames(p.value.table) = colnames(geTable)
	
	shuffled.mean.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
	rownames(shuffled.mean.table) = colnames(geTable)
	colnames(shuffled.mean.table) = colnames(geTable)
	
	a = combn(1:ncol(geTable),2)
	for (i in 1:ncol(a)){
		cor.scores = sapply(shuffled.cor.list,"[",a[1,i],a[2,i])
		shuffled.mean.table[a[1,i],a[2,i]] = mean(cor.scores)
		shuffled.mean.table[a[2,i],a[1,i]] = mean(cor.scores)
		p.value = mean(abs(cor.scores)>=abs(Corr.Coeff.Table[a[1,i],a[2,i]]))
		p.value.table[a[1,i],a[2,i]] = p.value
		p.value.table[a[2,i],a[1,i]] = p.value
		rm(list=c('cor.scores','p.value'))
		setTxtProgressBar(pb, (i*100)/ncol(a))
	}
	
	p.value.table[p.value.table==0] <- 0.000001
	neg.log10.p = -log10(p.value.table)
	
	#step8 "Overlap in Markers"
	#for all pairs of cell-types generate list of genes that are at least 1.5x avg in both cells
	
	#from above a = combn(1:ncol(geTable),2)
	marker.overlap.list = list()
	for (i in 1:ncol(a)){
		datasubset = cbind(geTable[,a[1,i]],geTable[,a[2,i]])
		markers = rownames(geTable[datasubset[,1]>1 & datasubset[,2]>1,])
		marker.overlap.list[[i]] = markers
		names(marker.overlap.list)[i] = paste(colnames(geTable)[a[1,i]], colnames(geTable)[a[2,i]],sep='_')
		rm(list=c('datasubset','markers'))
	}
	
	
	allgene.list.to.return = list(Corr.Coeff.Table,shuffled.mean.table,
																p.value.table,neg.log10.p,
																intersect_markers,
																rownames(geTable),
																length(intersect_markers),
																length(rownames(geTable)),nDESp1, nDESp2, geTable,marker.overlap.list)
	names(allgene.list.to.return) = c('corr.coeff','shuffled_correlation_score_means',
																		'p.value', 
																		'negative_log10_p.value',
																		'DEgenes_intersect',
																		'DEgenes_in_analysis',
																		'nDEgenes_intersect',
																		'nDEgenes_in_analysis',
																		'nDEgenes_Sp1','nDEgenes_Sp2',
																		'scaled_table','overlapping_markers')
	return(allgene.list.to.return)
}

select_GENIE3_Links <- function(linkList, corrMat, pMat, weightThreshold = 0.001, topThr = 0.01, 
																nTopTfs = c(5, 10, 20), nTopTargets = c(5, 10, 20), aFun = topPerTf, weightCol = "weight", verbose = TRUE) 
{
	
	if (!all(c("TF", "Target", weightCol) %in% colnames(linkList))) 
		stop("The link list colnames should be \"TF\", \"Target\", \"", 
				 weightCol, "\"")
	cntPairs <- table(table(linkList[, "TF"], linkList[, "Target"]))
	if (any(names(cntPairs) > 1)) 
		stop("There are duplicated regulator-target (gene id/name) pairs in the input link list.")
	msg <- paste0(format(Sys.time(), "%H:%M"), "\tCreating TF modules")
	
	linkList <- linkList[which(linkList[, weightCol] >= weightThreshold), ]
	if (verbose) {
		message("Number of links between TFs and targets (weight>=", 
						weightThreshold, "): ", nrow(linkList))
	}
	tfModules <- list()
	linkList$TF <- as.character(linkList$TF)
	linkList$Target <- as.character(linkList$Target)
	allName <- paste0("w", format(weightThreshold, scientific = FALSE))
	tfModules[[allName]] <- split(linkList$Target, factor(linkList$TF))
	if (!is.null(topThr)) {
		topThr <- setNames(topThr, paste("w", format(topThr, 
																								 scientific = FALSE), sep = ""))
		for (i in seq_along(topThr)) {
			llminW <- linkList[which(linkList[, weightCol] > 
															 	topThr[i]), ]
			tfModules[[names(topThr)[i]]] <- split(llminW$Target, 
																						 factor(llminW$TF))
		}
	}
	if (!is.null(nTopTargets)) {
		nTopTargets <- setNames(nTopTargets, paste("top", nTopTargets, 
																							 sep = ""))
		for (i in seq_along(nTopTargets)) {
			tfModules[[names(nTopTargets)[i]]] <- lapply(tfModules[[allName]], 
																									 function(x) x[1:(min(length(x), nTopTargets[i]))])
		}
	}
	tfModules.melted <- reshape2::melt(tfModules)
	colnames(tfModules.melted) <- c("Target", "TF", "method")
	topTFsperTarget.asDf <- NULL
	if (!is.null(nTopTfs)) {
		linkList_byTarget <- split(linkList, factor(linkList$Target))
		nTopTfs <- setNames(nTopTfs, paste("top", nTopTfs, "perTarget", 
																			 sep = ""))
		topTFsperTarget <- lapply(linkList_byTarget, function(llt) {
			nTFs <- nTopTfs[which(nTopTfs <= nrow(llt))]
			reshape2::melt(lapply(nTFs, function(x) llt[1:x, 
																									"TF"]))
		})
		topTFsperTarget <- topTFsperTarget[which(!sapply(sapply(topTFsperTarget, 
																														nrow), is.null))]
		topTFsperTarget.asDf <- data.frame(data.table::rbindlist(topTFsperTarget, 
																														 idcol = TRUE))
		colnames(topTFsperTarget.asDf) <- c("Target", "TF", "method")
	}
	byFun <- NULL
	if (!is.null(aFun)) {
		byFun <- aFun(linkList, weightCol)
	}
	tfModules <- rbind(tfModules.melted, topTFsperTarget.asDf, 
										 byFun)
	rm(tfModules.melted)
	rm(topTFsperTarget.asDf)
	tfModules$TF <- as.character(tfModules$TF)
	tfModules$Target <- as.character(tfModules$Target)
	if (verbose) 
		print(rbind(nTFs = length(unique(tfModules$TF)), nTargets = length(unique(tfModules$Target)), 
								nGeneSets = nrow(unique(tfModules[, c("TF", "method")])), 
								nLinks = nrow(tfModules)))
	
	if (!is.null(corrMat) & !is.null(pMat)) {
		tfs <- unique(tfModules$TF)
		missingTFs <- tfs[which(!tfs %in% rownames(corrMat))]
		if (length(missingTFs) > 0) {
			warning("The following TFs are missing from the correlation matrix: ", 
							paste(missingTFs, collapse = ", "))
			tfs <- tfs[which(tfs %in% rownames(corrMat))]
			corrMat <- corrMat[tfs, ]
		}
		tfModules_byTF <- split(tfModules, as.factor(tfModules$TF))
		tfModules_withCorr_byTF <- lapply(tfModules_byTF[tfs], 
																			function(tfGeneSets) {
																				tf <- as.character(unique(tfGeneSets$TF))
																				targets <- as.character(tfGeneSets$Target)
																				cbind(tfGeneSets, corr = as.numeric(corrMat[tf, targets]),
																							pvalue = as.numeric(pMat[tf, targets]))
																			})
		tfModules_withCorr_byTF <- tfModules_withCorr_byTF[which(lengths(tfModules_withCorr_byTF) > 
																														 	0)]
		tfModules_withCorr <- data.frame(data.table::rbindlist(tfModules_withCorr_byTF))
		if (length(missingTFs) > 0) {
			tfModules_withCorr <- rbind(tfModules_withCorr, data.frame(tfModules[tfModules$TF %in% 
																																					 	missingTFs, ], corr = NA))
		}
	} else {
		tfModules_withCorr <- data.frame(tfModules, corr = NA)
		if (verbose) 
			message("Correlation information not available. It will not be added to the modules.")
	}
	return(tfModules_withCorr)
}

runCorr <- function(exprMatr, method = "spearman") 
{
	res_corr <- psych::corr.test(exprMatr, method = method, ci = F)
	return(res_corr)
}

topPerTf <- function(ll, weightCol) 
{
	ll <- split(ll, factor(as.character(ll$TF)))
	tptf <- setNames(lapply(names(ll), function(tf) {
		tfMean <- mean(ll[[tf]][, weightCol])
		tfSd <- sd(ll[[tf]][, weightCol])
		list(top1sd = ll[[tf]][which(ll[[tf]][, weightCol] >= 
																 	tfMean + tfSd), "Target"], top3sd = ll[[tf]][which(ll[[tf]][, 
																 																															weightCol] >= tfMean + (3 * tfSd)), "Target"])
	}), names(ll))
	tptf <- reshape2::melt(tptf)
	colnames(tptf)[which(colnames(tptf) == "value")] <- "Target"
	colnames(tptf)[which(colnames(tptf) == "L1")] <- "TF"
	colnames(tptf)[which(colnames(tptf) == "L2")] <- "method"
	tptf <- tptf[, c("Target", "TF", "method")]
	return(tptf)
}

GenerateGCNByGENIE3  <- function(exp_count.seurat, cluster_markers, TF_symbol, targets_geneset = NULL, topTargets = 50, topTFs = 20, nsample_cells = 1000, seed_use = 1234, nCores = 6){
	
	top_targets <- (cluster_markers %>% group_by(cluster) %>% top_n(n = topTargets, wt = avg_log2FC))
	tf_markers <- cluster_markers[cluster_markers$gene %in% TF_symbol, ] 
	top_tfs <- (tf_markers %>% group_by(cluster) %>% top_n(n = topTFs, wt = avg_log2FC))
	
	celltypes <- unique(exp_count.seurat@active.ident)
	celltype_tfModules <- c()
	for (celltype_i in celltypes){
		
		message("Cell type ", celltype_i)
		celltype_object <- subset(exp_count.seurat, ident = celltype_i)
		if(ncol(celltype_object) > nsample_cells){
			## subset cells
			subcell <- sample(colnames(celltype_object), nsample_cells)
			scRNAsub <- celltype_object[, subcell]
		} else {
			scRNAsub <- celltype_object
		}
		
		combined_gene <- unique(c(targets_geneset, top_targets$gene[top_targets$cluster==celltype_i], top_tfs$gene[top_tfs$cluster==celltype_i]))
		exprMatr <- as.matrix(scRNAsub@assays$RNA@data[combined_gene, ])
		
		set.seed(seed_use) # For reproducibility of results
		regulators <-  top_tfs$gene[top_tfs$cluster==celltype_i]
		weightMat <- GENIE3(exprMatr, regulators=regulators, nCores = nCores)
		linkList <- getLinkList(weightMat)
		colnames(linkList) <- c("TF", "Target", "weight")
		# the input matrix of runCorr: row is cell, column is gene
		res_corr <- runCorr(t(exprMatr))
		corrMat <- res_corr$r
		pMat <- res_corr$p
		tfModules <- select_GENIE3_Links(linkList = linkList, corrMat = corrMat, pMat = pMat)
		tfModules <- tfModules[, c( "TF", "Target", "method", "corr", "pvalue")]
		tfModules$celltype <- celltype_i
		celltype_tfModules <- rbind(celltype_tfModules, tfModules) 
		
	}
	return(celltype_tfModules)
}

GenerateGCNByscLink <- function(exp_count.seurat, cluster_markers, TF_symbol, targets_geneset = NULL, topTargets = 50, topTFs = 20, 
																seed_use = 1234, nCores = 6, lda_value = seq(3, 0.01, -0.05)){
	
	top_targets <- (cluster_markers %>% group_by(cluster) %>% top_n(n = topTargets, wt = avg_log2FC))
	tf_markers <- cluster_markers[cluster_markers$gene %in% TF_symbol, ] 
	top_tfs <- (tf_markers %>% group_by(cluster) %>% top_n(n = topTFs, wt = avg_log2FC))
	
	celltypes <- unique(exp_count.seurat@active.ident)
	celltype_tfModules <- c()
	for (celltype_i in celltypes){
		
		message("Cell type ", celltype_i)
		celltype_object <- subset(exp_count.seurat, ident = celltype_i)
		combined_gene <- unique(c(targets_geneset, top_targets$gene[top_targets$cluster==celltype_i], top_tfs$gene[top_tfs$cluster==celltype_i]))
		# count.norm row is cell, column is gene
		count.norm = t(as.matrix(celltype_object@assays$RNA@data[combined_gene, ]))
		
		## filer cells and genes
		row_sum <- apply(count.norm > 0, 1, sum)
		col_sum <- apply(count.norm > 0, 2, sum)
		count.norm <- count.norm[row_sum > 0.01 * ncol(count.norm), col_sum > 0.01 * nrow(count.norm)]
		
		## network
		set.seed(seed_use)
		networks = scLink::sclink_net(expr = count.norm, ncores = nCores, lda = lda_value)
		
		bic_value <- c()
		for (i in 1:length(lda_value)){ bic_value <- c(bic_value, networks$summary[[i]]$bic)}
		#plot_bic <- data.frame(lda = lda_value, bic = bic_value)
		#print(ggplot(plot_bic) + geom_point(mapping = aes(x = lda_value, y = bic_value)))
		# select bic value 
		
		mat <- networks$summary[[which.min(bic_value)]]$adj
		mat[!upper.tri(mat, diag = TRUE)] <- 0
		net_adj <- reshape::melt(mat)
		net_adj <- net_adj[net_adj$value==1, ]
		net_adj <- net_adj[as.character(net_adj$X1)!=as.character(net_adj$X2), ]
		res_net <- list(net_adj = net_adj, corr = networks$cor, count.norm = count.norm)
		
		index1 <- res_net$net_adj$X1 %in% TF_symbol
		index2 <- res_net$net_adj$X2 %in% TF_symbol
		if(sum(index1) > 0 | sum(index2) > 0) {
			aa <- res_net$net_adj[index1, ]
			bb <- res_net$net_adj[index2, 	c(2, 1, 3)]
			colnames(aa) <- c("TF", "Target", "value")
			colnames(bb) <- c("TF", "Target", "value")
			res_net$net_adj <- rbind(aa, bb)
			res_net$net_adj$TF <- as.character(res_net$net_adj$TF)
			res_net$net_adj$Target <- as.character(res_net$net_adj$Target)
			# corr 
			
			net_genes <- unique(c(as.character(res_net$net_adj$TF), as.character(res_net$net_adj$Target)))
			res_net$corr <- res_net$corr[net_genes, net_genes]
			res_net$count.norm <- res_net$count.norm[ , net_genes]
			spearman <- psych::corr.test(res_net$count.norm, method = "spearman", ci = F)
			
			for(i in 1:nrow(res_net$net_adj)){
				res_net$net_adj$corr[i] <- res_net$corr[res_net$net_adj$TF[i], res_net$net_adj$Target[i]]
				res_net$net_adj$pvalue[i] <- spearman$p[res_net$net_adj$TF[i], res_net$net_adj$Target[i]]
			}
			
			res_net$net_adj$celltype <- celltype_i
			celltype_tfModules <- rbind(celltype_tfModules, res_net$net_adj)
		}
	}
	
	return(celltype_tfModules)
}

GetNetFromStringdb <- function(string_version = "11.5", species_id, score_thr = 400,  tf_genes, target_genes, OrgDb_id = "org.Dr.eg.db"){
	
	library(tidyverse)
	library(clusterProfiler)
	library(org.Dr.eg.db)
	library(STRINGdb)
	library(igraph)
	library(ggraph)
	
	# 创建STRINGdb对象
	string_db <- STRINGdb$new(version = string_version, species = species_id, score_threshold = score_thr, 
														network_type = "full", input_directory ="")
	
	combinedgene <- unique(c(tf_genes, target_genes))
	
	# 将Gene Symbol转换为Entrez ID
	gene <- combinedgene %>% bitr(fromType = "SYMBOL", 
																toType = "ENTREZID", 
																OrgDb = OrgDb_id, 
																drop = T)
	
	data_mapped <- string_db$map(gene, "ENTREZID", removeUnmappedRows = TRUE)
	data_links <- data_mapped$STRING_id %>% string_db$get_interactions()
	
	# stringID to Symbol，只取前两列和最后一列
	links <- data_links %>%
		mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
		mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
		dplyr::select(from, to , last_col()) %>% 
		dplyr::rename(weight = combined_score)
	
	if( sum(duplicated(links)) > 1){
		links <- links[duplicated(links), ]
	}
	
	links1 <- links[links$from %in% tf_genes, ]
	colnames(links1) <- c("TF", "Target", "weight")
	links2 <- links[links$to %in% tf_genes, c(2, 1, 3)]
	colnames(links2) <- c("TF", "Target", "weight")
	links_tf <- rbind(links1, links2)
	
	return(links_tf)
}

Classification_Seurat <- function(input_object,ref.seurat){
	
	input_object.anchors <- FindTransferAnchors(reference =ref.seurat, query = input_object, 
																							dims = 1:30,project.query = T,k.filter=150,k.anchor=5)
	predictions <- TransferData(anchorset = input_object.anchors, refdata = ref.seurat$celltype, 
															dims = 1:30)
	input_object <- AddMetaData(object = input_object, metadata = predictions)
	input_object$celltype_Seurat <- input_object$predicted.id
	return(input_object)
}

calc_var_ratio <- function(seurat_obj, species, genes_list) {
	
	# Create a dataframe with the same structure as df_var
	df_var <- data.frame(
		matrix(nrow = length(levels(seurat_obj@active.ident)) + 1, ncol = 1, dimnames = list(c("all cells", levels(seurat_obj@active.ident)), c(species))),
		stringsAsFactors = FALSE
	)
	
	# Calculate variance ratio for the specified species
	df <- apply(seurat_obj@assays$RNA@data, 1, var)
	var_tot <- sum(df)
	var_sub <- sum(df[match(genes_list[[species]], names(df))])
	cat(paste0("explained variance ", species, " - all cells: ", var_sub/var_tot, "\n"))
	df_var["all cells", species] <- var_sub/var_tot
	
	for (group in levels(seurat_obj@active.ident)) {
		temp <- subset(seurat_obj, idents = group)
		df <- apply(temp@assays$RNA@data, 1, var)
		var_tot <- sum(df)
		var_sub <- sum(df[match(genes_list[[species]], names(df))])
		cat(paste0("explained variance ", species, " - ", group, " cells: ", var_sub/var_tot, "\n"))
		df_var[group, species] <- var_sub/var_tot
		rm(temp)
	}
	
	return(df_var)
}


