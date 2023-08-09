##------20190724 ZJH 
##------20201205 ZJH update functions 
##------20210103 ZJH update functions
##------20221025 ZJH update functions

## load function 
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(MAST)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(Seurat)
library(openxlsx)
library("grid")
library("ggplotify")
library(magrittr)
library(dplyr)
library(ggsignif)
library("ggplot2")
library(cowplot)
library(openxlsx)

CreateFirstFilterSeuratObject <- function(exp_count, pr_name, Or_ident, 
                                          min_cell, min_gene, mt_pattern,rb_pattern){
# CreateSeuratObject : filter gene and cell by min.cells(Include features detected in at least this many cells) and min.features(Include cells where at least this many features are detected.)
exp_count.seurat <- CreateSeuratObject(exp_count, 
									project = pr_name, 
									min.cells = min_cell,
									min.features = min_gene)
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
folder_name="1.QC"
if (file.exists(folder_name)){
print("1.QC existed.")
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
                                                      "percent.mt","percent.rb"), ncol = 2, pt.size=0)) 

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

#copy files to 1.QC
result_files <- c(file_name)
file.copy(result_files, folder_name ,overwrite = T)
file.remove(result_files) 
print("Based on the Violin plot, you can filter the cell by the nFeature_RNA, nCount_RNA and/or  percent.mt rb")
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


SelectPCsByLogNormalize <- function(exp_count.seurat, file_prefix, scale_factor, nfeatures, npc_plot){
  
  folder_name="2.RunPCA"
  if (file.exists(folder_name)){
    print("2.RunPCA existed.")
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
  print(DimPlot(exp_count.seurat, reduction="pca", label = TRUE, pt.size=1.1, group.by="ident"))
  dev.off()
  
  #Copy files to 2.RunPCA
  result_files <- c(selected_pcs_name,pca_plot_name)
  file.copy(result_files, folder_name ,overwrite = T)#拷贝文件
  file.remove(result_files) #移除拷贝完的文件
  
  return(exp_count.seurat)
}

SelectPCsBySCTransform <- function(exp_count.seurat,condition,npc_plot){

folder_name="2.RunPCA"
if (file.exists(folder_name)){
print("2.RunPCA existed.")
}
else{
dir.create(folder_name)
} 

# Normalize Find VariableGene and Scale Data
exp_count.seurat <- SCTransform(exp_count.seurat,vars.to.regress = c("percent.rb","percent.mt"))
# RunPCA
exp_count.seurat <- RunPCA(exp_count.seurat,verbose = FALSE,  npcs = npc_plot)

# Select PCs:1:30pcs
selected_pcs_name <- paste0(condition,"_Selected_PCs.pdf")
pdf(file = selected_pcs_name, height = 8, width = 8)
print(ElbowPlot(object = exp_count.seurat, ndims = npc_plot, reduction = "pca"))
dev.off()

# Plot 
pca_plot_name <- paste0(condition, "_PCA_plot",".pdf")
pdf(pca_plot_name,8,8)
print(DimPlot(exp_count.seurat,reduction="pca",label = TRUE,pt.size=1.1,group.by="ident",shape.by="orig.ident"))
dev.off()

#Copy files to 2.RunPCA
result_files <- c(selected_pcs_name,pca_plot_name)
file.copy(result_files, folder_name ,overwrite = T)#拷贝文件
file.remove(result_files) #移除拷贝完的文件

return(exp_count.seurat)
}

# RunUMAP and RunTSNE and Then find cluster
PlotCluster <- function(exp_count.seurat, file_prefix, npc_used, k_param, resolution_number){
  
  folder_name="3.PlotCluster"
  if (file.exists(folder_name)){
    print("3.PlotCluster existed.")
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

# Find markers for each cluster
FindmarkerForCluster <- function(exp_count.seurat, file_prefix, min.pct, logfc.threshold, p_val_adj, mt_rb_pattern){
folder_name="4.MarkersInCluster"
if (file.exists(folder_name)){
  print("4.MarkersInCluster file existed.")
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

save_name <- paste0(file_prefix,"_MarkersInClusters.csv")
write.csv(cluster.markers,save_name)

file.copy(save_name, folder_name,overwrite = T)
file.remove(save_name)
return(cluster.markers)
}

# Top markers for each cluster
TopMarkersInCluster <- function(cluster.markers, file_prefix, top_num){
library(dplyr)
folder_name="4.MarkersInCluster"
if (file.exists(folder_name)){
  print("4.MarkersInCluster file existed.")
}
else{
  dir.create(folder_name)
} 
#将readsCountSM.markers传给group_by，group_by按cluster 排序，再将结果传给top_n，top_n按avg_logFC排序，显示每个类中的前两个
top_marker <- cluster.markers %>% group_by(cluster) %>% top_n(n = top_num, wt = avg_log2FC)
file_name=paste(file_prefix, "_top_marker", top_num,".csv",sep="")
write.csv(top_marker, file =file_name)
file.copy(file_name, folder_name,overwrite = T)
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
  return(exp_count.seurat)
  }

##-------------------------------------------------  8.Plot Feature  ------------------------------------------##
Plot_Features <- function(exp_count.seurat,Features_used,condition,height,width){
folder_name="5.Plot_Features"
if (file.exists(folder_name)){
print("5.Plot_Features file existed")
}
else{
dir.create(folder_name)
}

# VlnPlot
#name1 <- paste("VlnPlot_",condition,".pdf",sep="")
#pdf(file = name1,height = height, width = width)
#print(VlnPlot(object = exp_count.seurat, features = Features_used))
#dev.off()
#file.copy(name1, folder_name,overwrite = T)#拷贝文件
#file.remove(name1)

# FeaturePlot
#name2 <- paste("FeaturePlot_",condition,".pdf",sep="")
#pdf(file = name2,height = height, width = width)
#print(FeaturePlot(object = exp_count.seurat, features = Features_used,cols = c("lightgrey", "red"),min.cutoff = 0, max.cutoff = 1,order=T,reduction="tsne"))
#dev.off()
#file.copy(name2, folder_name,overwrite = T)#拷贝文件
#file.remove(name2)

#HeatMap
name3 <- paste("HeatMap_",condition,".pdf",sep="")
pdf(file = name3,height = height, width = width)
print(DoHeatmap(exp_count.seurat, features = Features_used) + NoLegend())
dev.off()
file.copy(name3, folder_name, overwrite = T)#拷贝文件
file.remove(name3)
}

##-----------------------------------------------------------------9. celltype identification----------------------------------------##
# CHETAH
Classification_CHETAH <- function(input_object,ref_object,n_gene_used)
{
									 
input_object_count <- input_object@assays$RNA@counts
input_object_reduced_tsne <- Embeddings(input_object, reduction = "tsne")#使用embeddings 函数来调用tsne 值								 

all.equal(rownames(input_object_reduced_tsne), colnames(input_object_count))
# 创建input_object
input_cell <- SingleCellExperiment(assays = list(counts = input_object_count),
                                  reducedDims = SimpleList(TSNE = input_object_reduced_tsne))
								  
## Classification
input_cell <- CHETAHclassifier(input = input_cell,
                              ref_cells = ref_object, n_genes = n_gene_used)

pdf("PlotCHETAH_confidence_score_0.1.pdf",18,12)
PlotCHETAH(input_cell)
dev.off()
input_object$celltype_CHETAH_confidence_score_0.1 <- as.factor(input_cell$celltype_CHETAH)
## confidence score =0
input_cell <- Classify(input_cell, 0)
pdf("PlotCHETAH_confidence_score_0.pdf",8,8)
PlotCHETAH(input = input_cell, tree = FALSE)
dev.off()
input_object$celltype_CHETAH_confidence_score_0 <- as.factor(input_cell$celltype_CHETAH)
return(input_object)
}

# Seurat

Make_ref_data_seurat <- function(ref_exp,ref_celltype){

ref.seurat <- CreateSeuratObject(counts = ref_exp)
ref.seurat <- SCTransform(ref.seurat)
ref.seurat <- RunPCA(object = ref.seurat, npcs = 30, verbose = FALSE)
ref.seurat$celltype <- ref_celltype
saveRDS(ref.seurat, file = "ref_seurat.rds")
return(ref.seurat)
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

## SingleR (http://comphealth.ucsf.edu/SingleR/).
#library(SingleR)
#singler = CreateSinglerSeuratObject(counts.file, annot, project.name,
# min.genes = 500, technology, species = "Human" (or "Mouse"), citation,
#  normalize.gene.length = F, min.cells = 2, npca = 10
#  regress.out = "nUMI", reduce.seurat.object = T)


#counts.file = 'GSE74923_L1210_CD8_processed_data.txt'
# This file was probably proccessed with Excel as there are duplicate gene names
                                      # (1-Mar, 2-Mar, etc.). They were removed manually.
#annot.file = 'GSE74923_L1210_CD8_processed_data.txt_types.txt' # a table with two columns 
#                                    # cell name and the original identity (CD8 or L1210)
#singler = CreateSinglerSeuratObject(counts.file, annot.file, 'GSE74923', 
#                                   variable.genes='de', regress.out='nUMI', 
#                                   technology='C1', species='Mouse', 
#                                   citation='Kimmerling et al.', reduce.file.size = F, 
#                                   normalize.gene.length = T)
#save(singler,file='GSE74923.RData'))

# The object can then be saved and uploaded to the SingleR web-app for further analysis and visualization or using functions available in the SingleR package (see vignette).
#save(singler,file=paste0(project.name,'.RData')

##-------------------------------------10. celltype hist between different condition----------------------------------------##
## sample_info <- as.factor(sample_info)
## cell_type <- as.factor(exp_count.seurat$celltype)

Sample_celltype_hist <- function(sample_info,cell_type,save_name1){
print("sample_info and cell_type must be factor")
Sample_celltype <- table(sample_info[cell_type == levels(cell_type)[1]])
col_name <-  levels(cell_type)[1]
len = length(levels(cell_type))
for (i in 2:len)
{
celltype <- table(sample_info[cell_type==levels(cell_type)[i]])
col_name_1 <- levels(cell_type)[i]
Sample_celltype <- cbind(Sample_celltype,celltype)
col_name <- cbind(col_name,col_name_1)
}
colnames(Sample_celltype) <- col_name
save_name1 <- paste0(save_name1,".csv")
write.csv(Sample_celltype,file = save_name1)
return(Sample_celltype)}

Plot_celltype_hist <- function(plot_data,pdf_name){
## ggplot2
library(ggplot2)
p1 <- ggplot(plot_data,aes(x=celltype,y=cell_number,fill=type))+geom_bar(stat="identity",position ="dodge")
p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) #移除网格线背景以及加坐标线
p1 <- p1 + theme(axis.text.x = element_text(size = 10, color = "black", face="bold",
           vjust = 0.8, hjust = 0.6, angle = 30))
pdf(paste0(pdf_name,".pdf"),8,8)
print(p1)
dev.off()
}


##-------------------------12.差异基因数目的barplot----------------##
### computer up/down gene number
gene_number_stats <- function(file, sheetnames, FC_cutoff, sig_method, sig_cutoff){
  
  if(sig_method == "p_val"){
    stats <- c()
    for(i in sheetnames){
      
      DEGs <- read.xlsx(file, sheet = i, colNames = T, rowNames = T)
      
      up_index <- DEGs$avg_log2FC > FC_cutoff & DEGs$p_val < sig_cutoff
      up_num <- nrow(DEGs[up_index, ])
      
      down_index <- DEGs$avg_log2FC < -FC_cutoff & DEGs$p_val < sig_cutoff
      down_num <- nrow(DEGs[down_index, ])
      
      tmp <- c(up_num, down_num, i)
      stats <- rbind(stats,tmp)
    }
  }
  else{
    
    stats <- c()
    for(i in sheetnames){
      
      DEGs <- read.xlsx(file, sheet = i, colNames = T, rowNames = T)
      
      up_index <- DEGs$avg_log2FC > FC_cutoff & DEGs$p_val_adj < sig_cutoff
      up_num <- nrow(DEGs[up_index, ])
      
      down_index <- DEGs$avg_log2FC < -FC_cutoff & DEGs$p_val_adj < sig_cutoff
      down_num <- nrow(DEGs[down_index, ])
      
      tmp <- c(up_num, down_num, i)
      stats <- rbind(stats, tmp)
    }
  }
  stats <- data.frame(stats)
  stats[,1] <- as.numeric(as.character(stats[,1]))
  stats[,2] <- as.numeric(as.character(stats[,2]))
  colnames(stats) <- c("up","down","celltype")
  return(stats)
}

###plot
Plot_function <- function(){
library(ggplot2)

plot_data <- data.frame(number=c(data$up,-1*data$down),
                        celltype=rep(data$celltype,2),
                        ud=c(rep("up",dim(data)[1]),rep("down",dim(data)[1])))

plot_data$celltype <- factor(plot_data$celltype,levels=c("all_cells",levels(Aggregated_seurat$celltype_assign)[-3]))

range(plot_data$number)

p1 <- ggplot(plot_data, aes(x=celltype, y=number, fill = ud))+
  geom_col(position = position_dodge(width = 0), width = 0.6, size = 0.3, colour = "black")+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = "black", fill = "transparent"),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "celltype", y = "number") +
  geom_hline(yintercept = 0, size = 0.3) +
  scale_fill_manual(values=c("#377EB8","#E41A1C"))+
  scale_y_continuous(breaks = seq(-500, 300, 100), labels = as.character(abs(seq(-500, 300, 100))), limits = c(-500, 300))

ggsave("20200802_Figure_4_Female_DEG_number_barplot.pdf", p1, width = 5, height = 4)
}

####-------画特定基因在特定细胞中表达的featureplot图--------------------##
FeaturePlot_SpecificGene <- function(subtype_object, intertest_DEGs, gene_num, pic_name){
  raw_exp <- FetchData(subtype_object, vars = intertest_DEGs$GeneName)
  raw_exp$celltype <- subtype_object$celltype
  new_exp <- matrix(0, nrow(raw_exp), ncol(raw_exp)-1)
  
  j = 1
  for (i in intertest_DEGs$Celltypes)
  {
    if(j <= gene_num){
      print(i)
      index <- raw_exp$celltype==i
      new_exp[index, j] <- raw_exp[index, j]
      j <- j + 1}
    else{
      break
    }
  }
  colnames(new_exp) <- colnames(raw_exp)[1:gene_num]
  rownames(new_exp) <- rownames(raw_exp)
  new_exp <- as.data.frame(new_exp)
  
  subtype_object <- AddMetaData(subtype_object, metadata = new_exp,
                                col.name = colnames(new_exp))
  metrics <- colnames(new_exp)
  p3 <- FeaturePlot(subtype_object, features = metrics, 
                    col = c("lightgrey", "red"),
                    split.by = "orig.ident",
                    reduction = "tsne", min.cutoff = 'q10')
  
  ggsave(pic_name, p3, width = 10, height = 40)
  return(subtype_object)
}

##------------找给定一个seurat中每类细胞类型中的实验组和对照组之间的差异基因--------------##
##--20201205 
FindDEGsFromCelltypes <- function(Aggregated_seurat, cell_type, treat, 
                                  control,logfc_threshold,
                                  min_pct_cell, file_prefix){
results_DEGs_mast <- list()
results_DEGs_bimod <- list()
results_DEGs_roc  <- list()
for (i in 1:length(cell_type)){
  ## 找每一类Obesity 和 Control的差异基因
  subset_seurat <- subset(x = Aggregated_seurat, subset = celltype_assign == cell_type[i])
  
  print(cell_type[i])
  ## wilcox
  print("Mast test")
  DEGs_treat_VS_Control_mast <- FindMarkers(subset_seurat,slot="data",
                                                   logfc.threshold = logfc_threshold,ident.1 = treat, 
                                                   ident.2 = control,min.pct = min_pct_cell,
                                                   group.by = "orig.ident",test.use = "MAST",
                                                   min.cells.group = 1)
  results_DEGs_mast[[i]] <- DEGs_treat_VS_Control_mast
  
  ## bimod
  print("Biomod test")
  DEGs_treat_VS_Control_bimod <- FindMarkers(subset_seurat,slot="data",
                                                 logfc.threshold = logfc_threshold,ident.1 = treat, 
                                                 ident.2 = control,min.pct = min_pct_cell,
                                                 group.by = "orig.ident",test.use = "bimod",
                                                 min.cells.group = 1)
  results_DEGs_bimod[[i]] <- DEGs_treat_VS_Control_bimod
  
  ## roc
  print("Roc test")
  DEGs_treat_VS_Control_roc <- FindMarkers(subset_seurat,slot="data",
                                                 logfc.threshold = logfc_threshold,ident.1 = treat, 
                                                 ident.2 = control,min.pct = min_pct_cell,
                                                 group.by = "orig.ident",test.use = "roc",
                                                 min.cells.group = 1)
  
  results_DEGs_roc[[i]] <- DEGs_treat_VS_Control_roc
  
}
## 把list的列名重命名
names(results_DEGs_mast) <-  as.character(cell_type)#重命名
names(results_DEGs_bimod) <-  as.character(cell_type)#重命名
names(results_DEGs_roc) <-  as.character(cell_type)#重命名

## 将结果写到excel中
write.xlsx(results_DEGs_mast, file = paste0(file_prefix, "_mast_DEGs_Treat_VS_Control.xlsx"),
          col.names = T, row.names = T)
write.xlsx(results_DEGs_bimod, file = paste0(file_prefix, "_biomod_DEGs_Treat_VS_Control.xlsx"),
          col.names = T, row.names = T)
write.xlsx(results_DEGs_roc, file = paste0(file_prefix, "_roc_DEGs_Treat_VS_Control.xlsx"),
           col.names = T, row.names = T)
}


FindDEGsFromCelltypesMasttest <- function(Aggregated_seurat, cell_type, treat, 
                                  control,logfc_threshold,
                                  min_pct_cell, file_prefix){
  results_DEGs_mast <- list()

  for (i in 1:length(cell_type)){
    ## 找每一类Obesity 和 Control的差异基因
    subset_seurat <- subset(x = Aggregated_seurat, subset = celltype_assign == cell_type[i])
    
    print(cell_type[i])
    ## wilcox
    print("MAST test")
    DEGs_treat_VS_Control_mast <- FindMarkers(subset_seurat,slot="data",
                                                logfc.threshold = logfc_threshold,ident.1 = treat, 
                                                ident.2 = control,min.pct = min_pct_cell,
                                                group.by = "orig.ident",test.use = "MAST",
                                                min.diff.pct = 0,
                                                min.cells.group = 1)
    results_DEGs_mast[[i]] <- DEGs_treat_VS_Control_mast
    
    
  }
  ## 把list的列名重命名
  names(results_DEGs_mast) <-  as.character(cell_type)#重命名
  
  ## 将结果写到excel中
  write.xlsx(results_DEGs_mast, file = paste0(file_prefix, "_mast_DEGs_Treat_VS_Control.xlsx"),
             col.names = T, row.names = T)
}


FindDEGsFromCelltypesWilcoxtest <- function(Aggregated_seurat, cell_type, treat, 
                                            control,logfc_threshold,
                                            min_pct_cell, file_prefix){
  results_DEGs_wilcox <- list()
  
  for (i in 1:length(cell_type)){
    ## 找每一类Obesity 和 Control的差异基因
    subset_seurat <- subset(x = Aggregated_seurat, subset = celltype_assign == cell_type[i])
    
    print(cell_type[i])
    ## wilcox
    print("Wilcox test")
    DEGs_treat_VS_Control_wilcox <- FindMarkers(subset_seurat,slot="data",
                                                logfc.threshold = logfc_threshold,ident.1 = treat, 
                                                ident.2 = control,min.pct = min_pct_cell,
                                                group.by = "orig.ident",test.use = "wilcox",
                                                min.cells.group = 1,
                                                min.diff.pct = 0)
                                                #max.cells.per.ident = 1200,
                                                #random.seed = 123)
    results_DEGs_wilcox[[i]] <- DEGs_treat_VS_Control_wilcox
    
    
  }
  ## 把list的列名重命名
  names(results_DEGs_wilcox) <-  as.character(cell_type)#重命名
  
  ## 将结果写到excel中
  write.xlsx(results_DEGs_wilcox, file = paste0(file_prefix, "_wilcox_DEGs_Treat_VS_Control.xlsx"),
             col.names = T, row.names = T)
}


####--------------------神经元分析经常画的一些图------------------------#####
## 聚类图和细胞marker图的参数
#reduction_type <- "umap"
#group_by_type <- "celltype_assign"
#pdf_prefix <- "Ob_female_6wk"
#celltype_markers <- c("Snap25", "Syt1", "Slc17a6", 
                      #"Slc32a1", "Gad1", "Gad2", 
                      #"Apoe", "Hcrt", "Pomc",
                      #"Agrp", "Cartpt", "Npy")
## 每类细胞数目比例图的plot_data
#plot_data <- table(Aggregated_seurat$orig.ident, Aggregated_seurat$celltype_assign)
#plot_data <- data.frame(Percentage = c(plot_data[1, ]/sum(plot_data[1, ]), 
                                       #plot_data[2, ]/sum(plot_data[2, ])),
                        #Treat = c(rep("Control", 4), rep("Ob", 4)),
                        #Type = rep(colnames(plot_data), 2))

#plot_data$Type <- factor(plot_data$Type,levels = c("Apoe", "GABA", "Glu", "Hybrid"))
## 印记基因热图所用的marker
#imprinted_genes <- c("Nap1l5", "Ndn", "Peg3", "Snrpn", "Ube3a", "Gnas")

NeuronSubtypePlots <- function(Aggregated_seurat, reduction_type, group_by_type, 
                              celltype_markers, plot_data, imprinted_genes)
  {
## 1.聚类图：总的聚类图和分开的聚类图
# 总的聚类图
p1 <- DimPlot(Aggregated_seurat, reduction=reduction_type, label = TRUE,
              pt.size = 0.8 , label.size = 4.5, repel=TRUE, 
              group.by = group_by_type)
pdf(paste0("F1A_", pdf_prefix, "_neuron_clustering_assign_celltyp.pdf"), 6, 5)
print(p1)
dev.off()

# 分开的聚类图
p2 <- DimPlot(Aggregated_seurat, reduction=reduction_type, label = TRUE,
              pt.size = 0.8 , label.size = 4.5, repel=TRUE, 
              group.by = group_by_type, split.by = "orig.ident")
pdf(paste0("F1B_", pdf_prefix, "_neuron_clustering_assign_celltyp_split.pdf"), 10, 5)
print(p2)
dev.off()

## 2.鉴定神经元用的marker的featureplot图和小提琴图
# 用的marker的featureplot图
p3 <- FeaturePlot(Aggregated_seurat, features = celltype_markers, ncol = 4,
                  reduction = reduction_type, cols = c("lightgrey", "red"),
                  pt.size = 0.8)
pdf(paste0("F2A_", pdf_prefix, "_featureplot_of_neurons_markers.pdf"), 16, 12)
print(p3)
dev.off()

# 用的marker的小提琴图
p4 <- VlnPlot(Aggregated_seurat, features = celltype_markers, ncol = 4, pt.size = 0)
pdf(paste0("F2B_", pdf_prefix, "_vlnplot_of_neurons_markers.pdf"), 16, 12)
print(p4)
dev.off()

## 3.各个神经元类别在对照组和实验组之间的比例图
p5 <- ggplot(plot_data, aes(x = Treat,y = Percentage, fill = Type)) + geom_bar(position="stack", stat="identity")
p5 <- p5 + theme_classic()
# p <- p + facet_wrap(~Sex)
pdf(paste0("F3_", pdf_prefix, "_Neuron_percentage.pdf"), 4, 4)
print(p5)
dev.off()

## 4.神经元印记基因的热图
cor_len <- length(unique(Aggregated_seurat$celltype_assign))
annotation_col <- data.frame(subtype = Aggregated_seurat$celltype_assign,
                             row.names = rownames(Aggregated_seurat@meta.data))
heatmap_byGroup <- FetchData(object = Aggregated_seurat,
                             vars = imprinted_genes,
                             slot = 'scale.data')
heatmap_byGroup <- cbind(heatmap_byGroup,annotation_col)
heatmap_byGroup <- arrange(heatmap_byGroup, subtype)
heatmap_byGroup_order <- heatmap_byGroup[, 1:length(imprinted_genes)]
colorCount=18
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
col1=getPalette(colorCount)[1:cor_len]
names(col1) <- levels(annotation_col$subtype)
ann_colors = list(group = col1)
pdf(paste0("F4_", pdf_prefix, "_Neuron_subtype_imprinted_DEGs.pdf"),9,3)
print(pheatmap(t(heatmap_byGroup_order), #fontsize_row=3, 
         colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10),
         breaks=seq(-1, 1, length.out = 10),
         treeheight_row=0.5, treeheight_col=1, 
         border_color='grey', cluster_cols = F,cluster_rows = F,
         fontsize_row = 18,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = F,show_rownames = T, scale='none',
         border=TRUE, angle_col = "0",
         fontsize = 10))
dev.off()
}



####------20201226 Pieplot of specific genes----------####
# Barplot
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    axis.text.x = element_blank()
  )

Pieplot <- function(plot_data){
  bp<- ggplot(plot_data, aes(x=gene, y=value, fill=color, label = label))+
    geom_bar(width = 0.5, stat = "identity", position = "stack") + facet_wrap(~treat)
  # bp <- bp + scale_y_continuous(limits=c(0, 100), breaks=seq(0,90,20)) 
  bp <- bp +  blank_theme
  bp <- bp + geom_text(size = 3, colour = "black", fontface = "bold",
                       position=position_stack(0.5))
  pie <- bp + coord_polar("y", start=0)
  pie <- pie + scale_fill_manual(values=c("#E69F00", "#56B4E9"))
  pie <- pie + theme(legend.title = element_blank(), legend.position = "top")
  return(pie)
}

Pieplot_Func <- function(seurat_object, used_features, celltype, pre_fix, width, height){
  
  CellPertage <- DotPlot(seurat_object, features = used_features,
                         group.by = "celltype", split.by = "orig.ident")$data
  id <- as.matrix(as.character(CellPertage$id))
  len <-  nrow(CellPertage)
  CellPertage$celltype <- matrix(unlist(strsplit(id, '_')), len, 2, byrow = TRUE)[1:len, 1]
  CellPertage$treat <- matrix(unlist(strsplit(id, '_')), len, 2, byrow = TRUE)[1:len, 2]
  CellPertage$colors <- NULL
  write.csv(CellPertage, paste0(pre_fix, "_cell_percentage.csv"))
  
  pdf(paste0("Pieplot_of_", pre_fix,"_cell_percentage.pdf"), width, height)
  for (j in celltype){
    print(j)
    CellPertage_temp <- CellPertage[CellPertage$celltype==j, ]
    library(ggplot2)
    plot_data <- data.frame(value = c(CellPertage_temp$pct.exp, 100-CellPertage_temp$pct.exp))
    plot_data$treat <- c(CellPertage_temp$treat,CellPertage_temp$treat)
    plot_data$color <- as.factor(c(rep("Expressed", length(CellPertage_temp$pct.exp)),
                                   rep("Unexpressed", length(CellPertage_temp$pct.exp))))
    genenames <- as.character(CellPertage_temp$features.plot)
    plot_data$gene <- as.factor(c(rep(genenames, 2)))
    plot_data$label <- paste0(round(plot_data$value, 2), "%")
    p_plots <- list()
    for (i in used_features){
      p_plot <- paste0("p_",i)
      index <- plot_data$gene==i
      temp <- plot_data[index, ]
      temp_plot <- Pieplot(temp)
      assign(p_plot,temp_plot)
      p_plots[[i]] <- get(p_plot)
    }
    plots <- plot_grid(plotlist=p_plots, labels = j)
    print(plots)
  }
  dev.off()
}

## Usage
### Pomc Cartpt Gal Trh
## features_1 <- c("Pomc", "Cartpt", "Gal", "Trh")
# Female
## celltype <- c("GABA1",  "GABA2", "GABA3", "GABA4", "GABA5", 
              #"GABA6", "GABA7", "Glu1", "Glu2", "Apoe")
## pre_fix <- "Female_Pomc_Cartpt_Gal_Trh"
## width <- 8
## height <- 5
##Pieplot_Func(female, features_1, celltype, pre_fix, width, height)

##------------20210104 Selected Cells and then reclustering--------------##
SelectCells <- function(Aggregated_seurat, n.cells, seed.use){
  ## select samples
  cellid <- rownames(Aggregated_seurat@meta.data)
  index_control <- Aggregated_seurat$orig.ident=="WT"
  control_cellid <- cellid[index_control] 
  
  index_ob <- Aggregated_seurat$orig.ident=="ob/ob"
  ob_cellid <- cellid[index_ob]
  
  set.seed(seed.use)
  subset_controlcellid <- sample(control_cellid, n.cells, replace=TRUE)
  subset_obcellid <- sample(ob_cellid, n.cells, replace=TRUE)
  select_cellid <- as.vector(c(subset_controlcellid, subset_obcellid))
  Aggregated_seurat$cellnames <- rownames(Aggregated_seurat@meta.data)
  subset_Object <- subset(Aggregated_seurat, 
                          cells = select_cellid)
  return(subset_Object)
}

Reclustering  <- function(Aggregated_seurat, npc_used, resolution_number, k_param){
  
  Aggregated_seurat <- SCTransform(Aggregated_seurat,
                                   vars.to.regress = c("percent.rp","percent.mt"))
  # RunPCA
  set.seed(123)
  Aggregated_seurat <- RunPCA(Aggregated_seurat, dims = 1:npc_used, 
                              verbose = FALSE, check_duplicates = FALSE)
  Aggregated_seurat <- RunTSNE(Aggregated_seurat, dims = 1:npc_used, 
                               verbose = FALSE, check_duplicates = FALSE)
  Aggregated_seurat <- RunUMAP(Aggregated_seurat, dims = 1:npc_used, 
                               verbose = FALSE)
  # FinderNeighbors 
  Aggregated_seurat <- FindNeighbors(Aggregated_seurat, dims = 1:npc_used, 
                                     verbose = FALSE, k.param = k_param)
  Aggregated_seurat <- FindClusters(Aggregated_seurat, verbose = FALSE,
                                    resolution = resolution_number)
  
  # RenameIdents
  levels_define <- as.numeric(levels(Aggregated_seurat))
  new.cluster.ids <- levels_define + 1
  names(new.cluster.ids) <- levels(Aggregated_seurat)
  Aggregated_seurat <- RenameIdents(Aggregated_seurat, new.cluster.ids)
  
  return(Aggregated_seurat)
}


ClusteringPlot <- function(Aggregated_seurat, prefix, pt.size, reduction_type,
                           resolution_number){
  # FIgure 1
  p1 <- DimPlot(Aggregated_seurat, reduction = reduction_type, label = TRUE,pt.size = pt.size, 
                label.size = 4, repel = TRUE, group.by = "ident") + NoLegend()
  p1 <- p1 + labs(title = paste0("Clustering resolution", resolution_number))
  p1 <- p1 + theme(plot.title = element_text(hjust = 0.5)) 
  
  p2 <- DimPlot(Aggregated_seurat, reduction = reduction_type, label = TRUE,
                pt.size = pt.size, label.size = 4, repel = TRUE,group.by = "celltype_assign") + NoLegend()
  p2 <- p2 + labs(title = "Celltype Identification") 
  p2 <- p2 + theme(plot.title = element_text(hjust = 0.5)) 
  
  pdf1_name <- paste(prefix, "F1_clustering", reduction_type, resolution_number,
                     ".pdf", sep = "_")
  
  plots <- plot_grid(p1, p2, ncol = 2, rel_widths = c(1.8,2), labels = c("A","B"), label_size = 20)
  pdf(pdf1_name, 15, 8)
  print(plots)
  dev.off()
  
  # FIgure 2
  p2 <- DimPlot(Aggregated_seurat, reduction = reduction_type, label = TRUE, 
                group.by = "celltype_assign", pt.size = pt.size, label.size = 4, 
                repel=TRUE, split.by = "orig.ident") + NoLegend()
  pdf2_name <- paste(prefix, "F2_clustering_split_by_sample",  reduction_type, resolution_number,
                     ".pdf", sep = "_")
  pdf(pdf2_name, 15, 8)
  print(p2)
  dev.off()
}

MarkersPlot <- function(Aggregated_seurat, markers, names, prefix, reduction_type){
  p_celltype_plots <- list()
  clusters <- levels(Aggregated_seurat@active.ident)
  for (i in 1:length(markers)){
    p_celltype <- paste0("p_celltype_",i)
    temp <- VlnPlot(object = Aggregated_seurat, features = markers[i], ncol = 1, pt.size = 0, ) 
    temp <- temp +  NoLegend()
    temp <- temp + labs(title=paste0(clusters[i], " ",markers[i],"(", names[i], ")"))
    temp <- temp + theme(axis.title.y=element_blank()) + theme(plot.title = element_text(hjust = 0.5))
    assign(p_celltype, temp)
    p_celltype_plots[[i]] <- get(p_celltype)
  }
  
  plots_celltype <- plot_grid(plotlist = p_celltype_plots, ncol = 2)
  
  pdf3_name <- paste0(prefix, "_F3_Vlnplot_marker_each_celltype", ".pdf")
  pdf(pdf3_name, 10, 15)
  print(plots_celltype)
  dev.off()
  
  pdf4_name <- paste(prefix,"F4_Featureplot_marker",".pdf", sep = "_")
  p <- FeaturePlot(Aggregated_seurat, features = markers,
                   reduction = reduction_type, ncol = 4, col = c("lightgrey", "red"))
  ggsave(pdf4_name, p, width = 20, height = 12)
}



# err_prob_list <- c(0.01, 0.05, 0.1, 0.15, 0.2)

CellnumberTest <- function(Aggregated_seurat, prefix, control_name, treat_name, err_prob_list, random_num){
  # 1. calculate cell number
  sample_info <- as.factor(Aggregated_seurat$orig.ident)
  cell_type <- Aggregated_seurat$celltype_assign
  obs.counts <- Sample_celltype_hist(sample_info,cell_type,
                                     paste0(prefix, "condition_cell_type"))
  
  # 2. cellnumber test
  print(obs.counts)
  res.table.ControlVsObesity = c()
  ## Go through a series of error probabilities
  for (err_prob in err_prob_list) {
    tip.exp <- generateNull(obs.counts, n=random_num, p=err_prob);
    ## Control vs Obesity
    res.1 = two.class.test(obs.counts, tip.exp, cond.control=control_name, 
                           cond.treatment=treat_name, to.plot=F)
    res.table.ControlVsObesity = rbind(res.table.ControlVsObesity, res.1)
  }
  
  rownames(res.table.ControlVsObesity) = as.character(paste0("error_p_", err_prob_list))
  res.table.ControlVsObesity <- t(res.table.ControlVsObesity)
  
  # 3.save results
  cellnumber <- t(obs.counts)
  cellnumber_pct <- data.frame(control_percent = round(cellnumber[,1]/sum(cellnumber[,1]),4) * 100,
                               obesity_percent = round(cellnumber[,2]/sum(cellnumber[,2]),4) * 100,
                               control_cellnumber = cellnumber[,1],
                               obesity_cellnumber = cellnumber[,2])
  results <- cbind(cellnumber_pct, res.table.ControlVsObesity)
  write.csv(results, paste(prefix, "cellnumber_test.csv", sep = "_"))
}


HeatmapPlot <- function(Aggregated_seurat, features_1, pdfname){
  annotation_col <- data.frame(Aggregated_seurat@active.ident, 
                               row.names=rownames(Aggregated_seurat@meta.data))
  colnames(annotation_col) <- 'group'
  heatmap_byGroup <- FetchData(object = Aggregated_seurat, vars = features_1, slot='scale.data')
  heatmap_byGroup <- cbind(heatmap_byGroup, annotation_col)
  heatmap_byGroup <- arrange(heatmap_byGroup, group)
  heatmap_byGroup_1 <- as.matrix(t(heatmap_byGroup[ ,1:dim(heatmap_byGroup)[2]-1]))
  anno_col <-  data.frame(group = heatmap_byGroup[ ,dim(heatmap_byGroup)[2]])
  colnames(heatmap_byGroup_1) <- rownames(anno_col)
  all.equal(rownames(anno_col), colnames(heatmap_byGroup_1))
  pdf(pdfname,12,20)
  pheatmap::pheatmap(heatmap_byGroup_1, #fontsize_row=3, 
                     color=colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(100),
                     breaks=seq(-2, 2, length.out = 100),
                     treeheight_row=10, treeheight_col=2, 
                     border_color='grey',cluster_cols = F,cluster_rows = F,
                     fontsize_row = 18,
                     annotation_col = anno_col,
                     #annotation_colors = ann_colors,
                     show_colnames = F,scale='none',
                     border=TRUE)
  dev.off()
}


#####--------Basic function to convert human to mouse gene names--------######
# from  https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  #require("biomaRt")
  #human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  #genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  #humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  #return(humanx)
}

######--------Basic function to convert mouse to human gene names#--------#########
#musGenes <- c("Hmmr", "Tlx3", "Cpeb4")

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  
  #require("biomaRt")
  #human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  #genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  #humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  #return(humanx)
}

genes <- convertMouseGeneList(musGenes)
##### 
#\' 使用R基本绘图函数绘制y轴不连续的柱形图
#\'
#\' 绘制y轴不连续的柱形图，具有误差线添加功能。断点位置通过btm和top参数设置，如果不设置，函数可自动计算合适的断点位置。
#\' @title gap.barplot function
#\' @param df 长格式的data.frame，即数据框中每一列为一组绘图数据。
#\' @param y.cols 用做柱形图y值的数据列（序号或名称），一列为一组。
#\' @param sd.cols 与y值列顺序对应的误差值的数据列（序号或名称）。
#\' @param btm 低位断点。如果btm和top均不设置，程序将自动计算和设置断点位置。
#\' @param top 高位断点。
#\' @param min.range 自动计算断点的阈值：最大值与最小值的最小比值
#\' @param max.fold 自动计算断点时最大值与下方数据最大值的最大倍数比
#\' @param ratio 断裂后上部与下部y轴长度的比例。
#\' @param gap.width y轴断裂位置的相对物理宽度（非坐标轴实际刻度）
#\' @param brk.type 断点类型，可设为normal或zigzag
#\' @param brk.bg 断点处的背景颜色
#\' @param brk.srt 断点标记线旋转角度
#\' @param brk.size 断点标记线的大小（长度）
#\' @param brk.col 断点标记线的颜色
#\' @param brk.lwd 断点标记线的线宽
#\' @param cex.error 误差线相对长度，默认为1
#\' @param ... 其他传递给R基本绘图函数barplot的参数
#\' @return 返回barplot的原始返回值，即柱形图的x坐标
#\' @examples
#\' datax <- na.omit(airquality)[,1:4]
#\' cols <- cm.colors(ncol(datax))
#\' layout(matrix(1:6, ncol=2))
#\' set.seed(0)
#\' for (ndx in 1:6){
#\'     dt <- datax[sample(rownames(datax), 10), ]
#\'     par(mar=c(0.5,2,0.5,0.5))
#\'     brkt <- sample(c(\'normal\', \'zigzag\'), 1)
#\'     gap.barplot(dt, col=cols, brk.type=brkt, max.fold=5, ratio=2)
#\' }
#\' @author ZG Zhao
#\' @export
gap.barplot <- function(df, y.cols=1:ncol(df), sd.cols=NULL, btm=NULL, top=NULL, min.range=10, max.fold=5, ratio=1, gap.width=1,
         brk.type='normal', brk.bg='white', brk.srt=135, brk.size=1, brk.col='black', brk.lwd=1, cex.error=1, ...){
  if (missing(df)) stop('No data provided.')
  if (is.numeric(y.cols)) ycol <- y.cols else ycol <- colnames(df)==y.cols
  if (!is.null(sd.cols))
    if (is.numeric(sd.cols)) scol <- sd.cols else scol <- colnames(df)==sd.cols
    ## Arrange data
    opts <- options()
    options(warn=-1)
    y <- t(df[, ycol])
    colnames(y) <- NULL
    if(missing(sd.cols)) sdx <- 0 else sdx <- t(df[, scol])
    sdu <- y + sdx
    sdd <- y - sdx
    ylim <- c(0, max(sdu) * 1.05)
    ## 如果没有设置btm或top，自动计算
    if (is.null(btm) | is.null(top)){
      autox <- .auto.breaks(dt=sdu, min.range=min.range, max.fold=max.fold)
      if (autox$flag){
        btm <- autox$btm
        top <- autox$top
      } else {
        xx <- barplot(y, beside=TRUE, ylim=ylim, ...)
        if (!missing(sd.cols)) errorbar(xx, y, sdu - y, horiz=FALSE, cex=cex.error)
        box()
        return(invisible(xx))
      }
    }
    ## Set up virtual y limits
    halflen <- btm - ylim[1]
    xlen <- halflen * 0.1 * gap.width
    v_tps1 <- btm + xlen        # virtual top positions
    v_tps2 <- v_tps1 + halflen * ratio
    v_ylim <- c(ylim[1], v_tps2)
    r_tps1 <- top               # real top positions
    r_tps2 <- ylim[2]
    ## Rescale data
    lmx <- summary(lm(c(v_tps1, v_tps2)~c(r_tps1, r_tps2)))
    lmx <- lmx$coefficients
    sel1 <- y > top
    sel2 <- y >=btm & y <=top
    y[sel1] <- y[sel1] * lmx[2] + lmx[1]
    y[sel2] <- btm + xlen/2
    sel1 <- sdd > top
    sel2 <- sdd >=btm & sdd <=top
    sdd[sel1] <- sdd[sel1] * lmx[2] + lmx[1]
    sdd[sel2] <- btm + xlen/2
    sel1 <- sdu > top
    sel2 <- sdu >=btm & sdu <=top
    sdu[sel1] <- sdu[sel1] * lmx[2] + lmx[1]
    sdu[sel2] <- btm + xlen/2
    ## bar plot
    xx <- barplot(y, beside=TRUE, ylim=v_ylim, axes = FALSE, names.arg=NULL, ...)
    ## error bars
    if(!missing(sd.cols)) errorbar(xx, y, sdu - y, horiz=FALSE, cex=cex.error)
    ## Real ticks and labels    
    brks1 <- pretty(seq(0, btm, length=10), n=4)
    brks1 <- brks1[brks1 >= 0 & brks1 < btm]
    brks2 <- pretty(seq(top, r_tps2, length=10), n=4)
    brks2 <- brks2[brks2 > top & brks2 <= r_tps2]
    labx <- c(brks1, brks2)
    ## Virtual ticks
    brks <- c(brks1, brks2 * lmx[2] + lmx[1])
    axis(2, at=brks, labels=labx)
    box()
    ## break marks
    pos <- par("usr")
    xyratio <- (pos[2] - pos[1])/(pos[4] - pos[3])
    xlen <- (pos[2] - pos[1])/50 * brk.size
    px1 <- pos[1] - xlen
    px2 <- pos[1] + xlen
    px3 <- pos[2] - xlen
    px4 <- pos[2] + xlen
    py1 <- btm
    py2 <- v_tps1
    rect(px1, py1, px4, py2, col=brk.bg, xpd=TRUE, border=brk.bg)
    x1 <- c(px1, px1, px3, px3)
    x2 <- c(px2, px2, px4, px4)
    y1 <- c(py1, py2, py1, py2)
    y2 <- c(py1, py2, py1, py2)
    px <- .xy.adjust(x1, x2, y1, y2, xlen, xyratio, angle=brk.srt*pi/90)
    if (brk.type=='zigzag'){
      x1 <- c(x1, px1, px3)
      x2 <- c(x2, px2, px4)
      if (brk.srt > 90){
        y1 <- c(y1, py2, py2)
        y2 <- c(y2, py1, py1)
      } else {
        y1 <- c(y1, py1, py1)
        y2 <- c(y2, py2, py2)
      }
    }
    if (brk.type=='zigzag') {
      px$x1 <- c(pos[1], px2, px1, pos[2], px4, px3)
      px$x2 <- c(px2, px1, pos[1], px4, px3, pos[2])
      mm <- (v_tps1 - btm)/3
      px$y1 <- rep(c(v_tps1, v_tps1 - mm, v_tps1 - 2 * mm), 2)
      px$y2 <- rep(c(v_tps1 - mm, v_tps1 - 2 * mm, btm), 2)
    }
    par(xpd=TRUE)
    segments(px$x1, px$y1, px$x2, px$y2, lty=1, col=brk.col, lwd=brk.lwd)
    options(opts)
    par(xpd=FALSE)
    invisible(xx)
}



######----FindConditionDEGs------####
FindConditionDEGs <- function(SeuratObject, FC.cutoff, ident1, 
                              ident2, minpct, grouptoby, imprinting_genename){
  DEGs_unfilter <- FindMarkers(SeuratObject,slot="data",
                               logfc.threshold = FC.cutoff,
                               ident.1 = ident1, ident.2 = ident2,
                               min.pct = minpct,group.by = grouptoby,test.use="wilcox")
  DEGs_unfilter$avg_log2FC <- log2(exp(DEGs_unfilter$avg_logFC))
  DEGs_unfilter$pct.ob <-  DEGs_unfilter$pct.1
  DEGs_unfilter$pct.wt <-  DEGs_unfilter$pct.2
  DEGs_unfilter$pct.1 <- NULL
  DEGs_unfilter$pct.2 <- NULL
  Up_sig_DEGs <- DEGs_unfilter[DEGs_unfilter$avg_log2FC >=0 & DEGs_unfilter$p_val_adj <= 0.05, ]
  Down_sig_DEGs <- DEGs_unfilter[DEGs_unfilter$avg_log2FC < 0 & DEGs_unfilter$p_val_adj <= 0.05, ]
  non_sig_DEGs <- DEGs_unfilter[DEGs_unfilter$p_val_adj > 0.05, ]
  non_sig_DEGs$FC_rank <- rep(0, nrow(non_sig_DEGs)[1])
  non_sig_DEGs$direction <- rep("notSig", nrow(non_sig_DEGs)[1])
  
  ## order
  if(nrow(Up_sig_DEGs)==0 & nrow(Down_sig_DEGs)==0){
    DEGs_order <- non_sig_DEGs
    index <- intersect(imprinting_genename, rownames(DEGs_order))
    DEGs_order$imprinting <- rownames(DEGs_order) %in% index
    Imprinting_DEGs <- DEGs_order[index, ]
    Imprinting_DEGs <- arrange(Imprinting_DEGs, direction, FC_rank)
  } else if(nrow(Up_sig_DEGs)==0){
    o <- order(Down_sig_DEGs$avg_logFC, decreasing = F)
    Down_sig_DEGs <- Down_sig_DEGs[o, ]
    Down_sig_DEGs$FC_rank <- 1:nrow(Down_sig_DEGs)[1]
    Down_sig_DEGs$direction <- rep("Down", nrow(Down_sig_DEGs)[1])
    DEGs_order <- rbind(Down_sig_DEGs, non_sig_DEGs)
    index <- intersect(imprinting_genename, rownames(DEGs_order))
    DEGs_order$imprinting <- rownames(DEGs_order) %in% index
    Imprinting_DEGs <- DEGs_order[index, ]
  } else if(nrow(Down_sig_DEGs)==0){
    o <- order(Up_sig_DEGs$avg_logFC, decreasing = T)
    Up_sig_DEGs <- Up_sig_DEGs[o, ]
    Up_sig_DEGs$FC_rank <- 1:nrow(Up_sig_DEGs)[1]
    Up_sig_DEGs$direction <- rep("Up", nrow(Up_sig_DEGs)[1])
    DEGs_order <- rbind(Up_sig_DEGs, non_sig_DEGs)
    index <- intersect(imprinting_genename, rownames(DEGs_order))
    DEGs_order$imprinting <- rownames(DEGs_order) %in% index
    Imprinting_DEGs <- DEGs_order[index, ]
  } else {
    o <- order(Up_sig_DEGs$avg_logFC, decreasing = T)
    Up_sig_DEGs <- Up_sig_DEGs[o, ]
    Up_sig_DEGs$FC_rank <- 1:nrow(Up_sig_DEGs)[1]
    Up_sig_DEGs$direction <- rep("Up", nrow(Up_sig_DEGs)[1])
    
    o <- order(Down_sig_DEGs$avg_logFC, decreasing = F)
    Down_sig_DEGs <- Down_sig_DEGs[o, ]
    Down_sig_DEGs$FC_rank <- 1:nrow(Down_sig_DEGs)[1]
    Down_sig_DEGs$direction <- rep("Down", nrow(Down_sig_DEGs)[1])
    
    DEGs_order <- rbind(Up_sig_DEGs, Down_sig_DEGs)
    DEGs_order <- rbind(DEGs_order, non_sig_DEGs)
    index <- intersect(imprinting_genename, rownames(DEGs_order))
    DEGs_order$imprinting <- rownames(DEGs_order) %in% index
    index1 <- intersect(imprinting_genename, rownames(Up_sig_DEGs))
    index2 <- intersect(imprinting_genename, rownames(Down_sig_DEGs))
    a <- Up_sig_DEGs[index1, ]
    b <- Down_sig_DEGs[index2, ]
    Imprinting_DEGs <- rbind(a[order(a$FC_rank), ], b[order(b$FC_rank), ])
  }
  
  results <- list(DEGs_order = DEGs_order,
                  Imprinting_DEGs = Imprinting_DEGs)
  return(results)
}

###########-----20210316 构建GSEA genelist and genesets------------------########

# library(ReactomePA)
# library(reactome.db)
library(openxlsx)
library(clusterProfiler)

GenerateGSEAGeneSymbolList <- function(filename, i){
  # 1.读入数据
  temp <- read.xlsx(filename, colNames = T, rowNames = T, sheet = i)
  head(temp)
  gene <- rownames(temp)
  print(length(gene))
  print(head(gene))
  
  # 2.添加log2FC
  gene_df <- data.frame(avg_log2FC=temp[gene, ]$avg_log2FC, #可以是foldchange
                        SYMBOL = gene) #记住你的基因表头名字
  head(gene_df)

  # 3.按avg_log2FC降序排
  geneList <- gene_df$avg_log2FC #第二列可以是folodchange，也可以是logFC
  names(geneList)=gene_df$SYMBOL #使用转换好的ID
  geneList=sort(geneList,decreasing = T) #从高到低排序
  head(geneList)
  return(geneList)
}

GenerateGSEAGeneList <- function(filename, i){
  # 1.读入数据
  temp <- read.xlsx(filename, colNames = T, rowNames = T, sheet = i)
  head(temp)
  gene <- rownames(temp)
  print(length(gene))
  print(head(gene))
  
  # 2.开始ID转换
  gene=bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  length(gene$ENTREZID)
  head(gene)
  
  # 3.去重
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  gene_df <- data.frame(avg_log2FC=temp[gene$SYMBOL, 4], #可以是foldchange
                        SYMBOL = gene$SYMBOL) #记住你的基因表头名字
  head(gene_df)
  gene_df <- merge(gene_df, gene, by="SYMBOL")
  head(gene_df)
  # 4.按avg_log2FC降序排
  geneList <- gene_df$avg_log2FC #第二列可以是folodchange，也可以是logFC
  names(geneList)=gene_df$ENTREZID #使用转换好的ID
  geneList=sort(geneList,decreasing = T) #从高到低排序
  head(geneList)
  return(geneList)
}


GenerateMmMsigDB <- function(rds_object){
  MmMsigDB <- data.frame(term = "term", gene = "NA")
  for (i in 1:length(rds_object)){
    gene_id <- rds_object[[i]]
    gene_df <- bitr(gene_id, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db") #会有部分基因数据丢失，或者ENSEMBL
    gene_name <- gene_df$SYMBOL
    temp_name <- names(rds_object)[i]
    temp_MsigDB <- data.frame(term = temp_name, gene = gene_name)
    MmMsigDB <- rbind(MmMsigDB, temp_MsigDB)
  }
  return(MmMsigDB)
}



#######----------20210324 Change FC rank------------###
ChangeFCrank <- function(DEGs_table, log2FC_cutoff, padj_cutoff){
  index <- abs(DEGs_table$avg_log2FC) >= log2FC_cutoff & DEGs_table$p_val < padj_cutoff
  allDEGs <- DEGs_table[index, ]
  # Up
  allDEGs_up <- allDEGs[allDEGs$avg_log2FC >= 0, ]
  allDEGs_up <- arrange(allDEGs_up, desc(avg_log2FC))
  allDEGs_up$FC_rank <- 1:nrow(allDEGs_up)
  allDEGs_up$direction <- "Up"
  #Down
  allDEGs_down <- allDEGs[allDEGs$avg_log2FC < 0, ]
  allDEGs_down <- arrange(allDEGs_down, avg_log2FC)
  allDEGs_down$FC_rank <- 1:nrow(allDEGs_down)
  allDEGs_down$direction <- "Down"
  # notSig
  notsigrows <- setdiff(rownames(DEGs_table), rownames(allDEGs))
  notsig <- DEGs_table[notsigrows, ]
  DEGs_table_reorder <- rbind(allDEGs_up, allDEGs_down)
  DEGs_table_reorder <- rbind(DEGs_table_reorder, notsig)
  return(DEGs_table_reorder)
}


#############------cellphonedb related function----------------------#######
GenerateDotPlotTable <- function(mypvals, mymeans, gene_pattern,column_selected_pattern){
  costimulatory <- grep(gene_pattern, mymeans$interacting_pair,value = T)
  print(costimulatory)
  mymeans %>% dplyr::filter(interacting_pair %in% costimulatory)%>%
    dplyr::select("interacting_pair",starts_with(column_selected_pattern),)  %>%  
    reshape2::melt() -> meansdf
  colnames(meansdf)<- c("interacting_pair","CC","means")
  print(dim(meansdf))
  mypvals %>% dplyr::filter(interacting_pair %in% costimulatory)%>%  
    dplyr::select("interacting_pair",starts_with(column_selected_pattern))  %>%
    reshape2::melt()-> pvalsdf
  colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
  print(dim(pvalsdf))
  pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
  meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
  pldf <- merge(pvalsdf,meansdf,by = "joinlab")
  return(pldf)
}

heatmat_plot <- function(count_network_file){
  count_net <- read.table(count_network_file, header = T)
  celltype <- c("GABA1", "GABA2", "GABA3", "GABA4", "GABA5", "GABA6",
                "Glu1", "Glu2", "Glu3")
  len <- length(celltype)
  heatmap_table <- matrix(0, len, len)
  rownames(heatmap_table) <- celltype
  colnames(heatmap_table) <- celltype
  heatmap_table
  for(i in celltype){
    index <- which(count_net$SOURCE %in% i)
    index_col <- count_net$TARGET[index]
    heatmap_table[i, index_col] <- count_net$count[index]
  }
  p <- pheatmap::pheatmap(heatmap_table, cluster_cols = F, cluster_rows = F)
  p <- as.ggplot(p)
  return(p)
}

#############------20210727 CalculateReadCountBarcode----------------------#######
CalculateReadCountBarcode <- function(mol.info.file,  row.prefix){
  
  molinfo <- read10xMolInfo(mol.info.file)
  selected_info <- molinfo$data[1:5]
  selected_info$cell <- factor(selected_info$cell)
  cell_read_sum <- tapply(selected_info$reads,selected_info$cell,sum)
  
  cell_read_sum <- as.matrix(cell_read_sum)
  rownames(cell_read_sum) <- paste0(row.prefix,  rownames(cell_read_sum))
  colnames(cell_read_sum) <- "reads_sum"
  print(dim(cell_read_sum)) 
  print(head(cell_read_sum))
  return(cell_read_sum)
}
##########------------20210805  Volcono plot function-------------------########
GenerateVolcanoPlotData <- function(filename, sheetnum, log2FC_cutoff, padj_value){
  
  library(openxlsx)
  DEG_data <- read.xlsx(filename, sheet = sheetnum, rowNames = T, colNames = T)
  DEG_data$log2FoldChange <- log2(exp(DEG_data$avg_logFC))
  DEG_data$padj <- DEG_data$p_val_adj
  
  # 2.padj是0的基因 因为在画火山图的时候 -log10(DEG_data$padj2) 会是无穷大 
  # 所以随机赋值3.81405601117434E-120, 3.81405601117434E-100 之间的数值
  sum_0 <- sum(DEG_data$padj == 0)
  if (sum_0 > 0){
    x1 <- runif(sum_0, 3.81405601117434E-120, 3.81405601117434E-100)
    DEG_data$padj[DEG_data$padj == 0] <- x1 
    DEG_data$logP <- -log10(DEG_data$padj)
  } else {
    DEG_data$logP <- -log10(DEG_data$padj) # 对差异基因矫正后p-value进行log10()转换
  }
  
  # 3. 将基因分为三类：not-siginficant，up，dowm
  
  DEG_data$Group <- "notSig"
  up_index <- which((DEG_data$padj < padj_value) & DEG_data$log2FoldChange > log2FC_cutoff)
  down_index <- which((DEG_data$padj < padj_value) & DEG_data$log2FoldChange < -log2FC_cutoff)
  DEG_data$Group[up_index] = "up"
  DEG_data$Group[down_index] = "down"
  print(table(DEG_data$Group))
  return(DEG_data)
}

GenerateVolcanoPlotData_DESeq2 <- function(filename, sheetnum, log2FC_cutoff, padj_value){
  
  library(openxlsx)
  DEG_data <- read.xlsx(filename, sheet = sheetnum, rowNames = T, colNames = T)
  # DEG_data$log2FoldChange <- log2(exp(DEG_data$avg_logFC))
  # DEG_data$padj <- DEG_data$p_val_adj
  
  # 2.padj是0的基因 因为在画火山图的时候 -log10(DEG_data$padj2) 会是无穷大 
  # 所以随机赋值3.81405601117434E-120, 3.81405601117434E-100 之间的数值
  sum_0 <- sum(DEG_data$padj == 0)
  if (sum_0 > 0){
    x1 <- runif(sum_0, 3.81405601117434E-120, 3.81405601117434E-100)
    DEG_data$padj[DEG_data$padj == 0] <- x1 
    DEG_data$logP <- -log10(DEG_data$padj)
  } else {
    DEG_data$logP <- -log10(DEG_data$padj) # 对差异基因矫正后p-value进行log10()转换
  }
  
  # 3. 将基因分为三类：not-siginficant，up，dowm
  
  DEG_data$Group <- "notSig"
  up_index <- which((DEG_data$padj < padj_value) & DEG_data$log2FoldChange > log2FC_cutoff)
  down_index <- which((DEG_data$padj < padj_value) & DEG_data$log2FoldChange < -log2FC_cutoff)
  DEG_data$Group[up_index] = "up"
  DEG_data$Group[down_index] = "down"
  print(table(DEG_data$Group))
  return(DEG_data)
}

DrawVolcanoPlot <- function(DEG_data, Interest_genelist, title_name){
  library(ggpubr)
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
  DEG_data$SpecificGene <- rownames(DEG_data) %in% Interest_genelist
  up_gene <- rownames(DEG_data)[DEG_data$SpecificGene=="TRUE" & DEG_data$Group == "up"]
  down_gene <- rownames(DEG_data)[DEG_data$SpecificGene=="TRUE" & DEG_data$Group == "down"]
  
  theme <- theme(panel.grid=element_blank(),
                        panel.background = element_rect(fill="white",color="black"),
                        legend.key = element_rect(fill = "white"),
                        legend.text = element_text(size = 15),
                        axis.title.x = element_text(size = 15),
                        axis.title.y = element_text(size = 15),
                        axis.text.x = element_text(size = 13),
                        axis.text.y = element_text(size = 13))
  
  x_lim <- ceiling(max(abs(DEG_data$log2FoldChange)) + 0.5)
  y_lim <- ceiling(max(DEG_data$logP) + 50)
  
  
  p <- ggplot(DEG_data,aes(x=log2FoldChange, y=logP, color = Group))+geom_point(size=0.9)+
    labs(x="log2(fold change)",y="-log10(adjusted P-value)",title=title_name) +
    scale_color_manual("",values =c(down="#0072B5",notSig="grey",up="#BC3C28"))+
    scale_x_continuous(limits = c(-x_lim,x_lim),breaks=seq(-x_lim, x_lim,0.5))+
    scale_y_continuous(limits = c(0, y_lim),breaks=seq(0, y_lim, 50))+
    geom_hline(yintercept=-log10(padj_value), linetype=3)+
    geom_vline(xintercept=c(-log2FC_cutoff,log2FC_cutoff),linetype=3) + theme+
    theme(legend.position=c(0.85,0.85),
          plot.title = element_text(hjust = 0.5,size = 15))+
    #nudge_x nudge_y 更改标签与点的距离
    # up gene label
    geom_text_repel(data = DEG_data[up_gene, ],
                    aes(x=DEG_data[up_gene, ]$log2FoldChange+0.01,
                        y=-log10(DEG_data[up_gene,]$padj)+0.01,
                        label = rownames(DEG_data[up_gene, ])),
                    size = 3.5,color="black",nudge_x=+0.2,nudge_y=1) +
    ## down gene label
    geom_text_repel(data = DEG_data[down_gene,],
                    aes(x=DEG_data[down_gene,]$log2FoldChange-0.01,
                        y=-log10(DEG_data[down_gene,]$padj)-0.01,label = rownames(DEG_data[down_gene,])),
                    size = 3.5,color="black",nudge_x=-0.3,nudge_y=10)
  return(p)
}


##########---20210821 热图--------------###

HeatmapPlotSpecific <- function(plot_object, features_1, pdfname_1, pdfname_2){
  annotation_col <- data.frame( group = plot_object$orig.ident,
                                sex = plot_object$group,
                                smallsubtype = plot_object$small_celltype,
                                row.names = rownames(plot_object@meta.data))
  
  heatmap_byGroup <- FetchData(object=plot_object,
                               vars=features_1,
                               slot='data')
  
  heatmap_byGroup <- cbind(heatmap_byGroup, annotation_col)
  heatmap_byGroup <- arrange(heatmap_byGroup, group, sex, smallsubtype)
  
  heatmap_byGroup_order <- heatmap_byGroup[, 1:length(features_1)]
  anno_col <- data.frame( smallsubtype = factor(heatmap_byGroup$smallsubtype),
                          sex = factor(heatmap_byGroup$sex),
                          group = factor(heatmap_byGroup$group))
  
  rownames(anno_col) <- rownames(heatmap_byGroup_order)
  colorCount=18
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  cols <- getPalette(colorCount)[1:14]
  col1 <- cols[1:2]
  col2 <- c("#FF9289", "#00DAE0")
  col3 <- cols[3:11]
  col4 <- cols[12:13]
  
  names(col1) <- levels(anno_col$sex)
  names(col3) <- levels(anno_col$smallsubtype)
  names(col4) <- levels(anno_col$group)
  ann_colors = list( smallsubtype = col3, 
                     sex = col1,
                     group = col4)
  
  
  library(pheatmap)
  p1 <- pheatmap(t(heatmap_byGroup_order), #fontsize_row=3, 
                 colorRampPalette(c("#CD00CD","black","yellow"))(10),
                 breaks=seq(-1.5, 1.5, length.out = 10),
                 treeheight_row=0.5, treeheight_col=1, 
                 border_color='grey', cluster_cols = F,cluster_rows = T,
                 fontsize_row = 18,
                 annotation_col = anno_col,
                 annotation_colors = ann_colors,
                 show_colnames = F,show_rownames = T, scale='row',
                 border=TRUE, angle_col = "0",
                 fontsize = 10) #, gaps_col = 7104)
  
  pdf(pdfname_1, 12, 50)
  print(p1)
  dev.off()
  
  p2 <- pheatmap(t(heatmap_byGroup_order), #fontsize_row=3, 
                 colorRampPalette(c("#CD00CD","black","yellow"))(10),
                 breaks=seq(-1.5, 1.5, length.out = 10),
                 treeheight_row=0.5, treeheight_col=1, 
                 border_color='grey', cluster_cols = F,cluster_rows = T,
                 fontsize_row = 18,
                 annotation_col = anno_col,
                 annotation_colors = ann_colors,
                 show_colnames = F,show_rownames = F, scale='row',
                 border=TRUE, angle_col = "0",
                 fontsize = 10) #, gaps_col = 7104)
  
  pdf(pdfname_2, 12, 8)
  print(p2)
  dev.off()
  
  res <- list(p_genename = p1,
              p_nogenename = p2)
  return(res)
}
