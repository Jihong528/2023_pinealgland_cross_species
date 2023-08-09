#!/home/zhengjh/miniconda3/envs/r403/bin/Rscript
# parameter

library(optparse)
library(getopt)

option_list <- list(
	make_option(c("-o", "--output"), type = "character", default = FALSE,
							action = "store", help = "This is the output directory."
	),
	make_option(c("--species1_object"), type = "character", default = FALSE,
							action = "store", help = "This is absoulte path of species1_object"
	),
	make_option(c("--species2_object"), type = "character", default = FALSE,
							action = "store", help = "This is absoulte path of species2_object"
	),
	make_option(c("--species1_marker"), type = "character", default = FALSE,
							action = "store", help = "This is absoulte path of species1_marker"
	),
	make_option(c("--species2_marker"), type = "character", default = FALSE,
							action = "store", help = "This is absoulte path of species2_marker"
	),
	make_option(c("--species1_name"), type = "character", default = FALSE,
							action = "fileprefix", help = "This is species1_name."),
	make_option(c("--species2_name"), type = "character", default = FALSE,
									action = "fileprefix", help = "This is species2_name.")
)

# -help
opt = parse_args(OptionParser(option_list = option_list, usage = "cal_cell_similarity.R [-o] [--species1_object] [--species2_object] [--species1_marker] [--species2_marker]......"))
print(opt)

timestart <- Sys.time() 

cat("Load used R package\n")
library(openxlsx)
library(corrplot)
library(Seurat)
library(SeuratObject)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
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

##### folder name 
cat("Set the output path\n")

folder_name <- "CelltyeSimilarityFigs"

if (file.exists(folder_name)) {
	print("CelltyeSimilarityFigs existed.")
}else{
	print("Create folder CelltyeSimilarityFigs.")
	dir.create(folder_name)
} 

##### step1 
setwd(opt$output)
cat("Step1: Read seurat objects and cell type markers of each speices\n")

species1_object <- readRDS(opt$species1_object)
species2_object <- readRDS(opt$species2_object)

cat("Read markers for species1_object\n")
species1_marker <- read.csv(opt$species1_marker,  row.names = 1, header = T)

cat("Read markers for species1_object\n")
species2_marker <- read.csv(opt$species2_marker,  row.names = 1, header = T)

cat("The number of markers of each cell type")
print(table(species1_marker$cluster))
print(table(species2_marker$cluster))

##### step2
cat("Step2: calculate the spearman correlation of cell type pairs\n")
allgene_avg_res <- AvgExpSpecies(species1_marker, species2_marker, species1_object, species2_object,  opt$species1_name, opt$species2_name)
Sp1 <- allgene_avg_res$species1_avgexp
Sp2 <- allgene_avg_res$species2_avgexp
intersect_markers <- allgene_avg_res$intersect_markers
nDESp1 <- allgene_avg_res$nDESp1
nDESp2 <- allgene_avg_res$nDESp2
allgene_corr_res <- CorrCompareFunction(Sp1, Sp2, intersect_markers, nDESp1, nDESp2)
saveRDS(allgene_corr_res, paste(opt$species1_name, opt$species2_name, "celltype_pair_spearman_correlation.rds", sep = "_"))

##### step3
cat("Step3: use heatmap to visualize the spearman correlation of cell type pairs\n")
allgene_corr_res$p.value[is.na(allgene_corr_res$p.value)] <- 0
plotdata <- as.matrix(allgene_corr_res$corr.coeff * (allgene_corr_res$p.value < 0.001))
plotdata <- plotdata[grep(opt$species1_name, rownames(plotdata)), grep(opt$species2_name, rownames(plotdata))]


col1 <- colorRampPalette(c("darkblue", "white","darkred"))

write.csv(plotdata, paste( opt$species1_name, opt$species2_name, "celltype_pair_plotdata.csv", sep = "_"))

pdf(paste(opt$species1_name, opt$species2_name, "celltype_pair_correlation_heatmap.pdf", sep = "_"), 12, 12)
print(pheatmap(plotdata, 
				 colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(50),
				 breaks=c(seq(-1, 1, length.out = 50)),
				 border_color='black', cluster_cols = T,cluster_rows = T))
dev.off()

##### step4

print("try")
plotdata <- as.data.frame(plotdata)
plotdata$row <- rownames(plotdata)
dd <- reshape2::melt(plotdata, id.vars = "row")
colnames(dd) <- c("species1", "species2", "value")
dd <- dd[dd$value !=0 , ]
overlap_index <- which(names(allgene_corr_res$overlapping_markers) %in% paste(dd$species1, dd$species2, sep = "_"))
overlapgene <- allgene_corr_res$overlapping_markers[overlap_index]
overlapgene <- do.call(cbind, lapply(lapply(overlapgene, unlist), 
																		 `length<-`, max(lengths(overlapgene))))
write.xlsx(overlapgene, paste( opt$species1_name, opt$species2_name, "sig_celltype_pair_overlapgene.xlsx", sep = "_"), 
					 colNames = T, rowNames = T)


files <- c(paste(opt$species1_name, opt$species2_name, "celltype_pair_spearman_correlation.rds", sep = "_"),
					 paste(opt$species1_name, opt$species2_name, "celltype_pair_plotdata.csv", sep = "_"),
					 paste(opt$species1_name, opt$species2_name, "celltype_pair_correlation_heatmap.pdf", sep = "_"),
					 paste( opt$species1_name, opt$species2_name, "sig_celltype_pair_overlapgene.xlsx", sep = "_"))
file.copy(files, folder_name ,overwrite = T)
file.remove(files)


print("Happy~Finished!!")
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)



