#!/home/zhengjh/miniconda3/envs/r403/bin/Rscript
# parameter

library(optparse)
library(getopt)

option_list <- list(
	make_option(c("-o", "--output"), type = "character", default = FALSE,
							action = "store", help = "This is the output directory."
	),
	make_option(c( "--seurat_object_path"), type = "character", default = FALSE,
							action = "store", help = "This is absoulte path of seuratobject."
	),
	make_option(c("--lrdb_path"), type = "character", default = FALSE,
							action = "store", help = "This is path of ligand-receptors db object."
	),
	make_option(c("--filelabel"), type = "character", default = FALSE,
							action = "store", help = "This is filelabel"),
	make_option(c("--thresh_p"), type = "double", default = FALSE,
							action = "store", help = "This is filelabel")
)

# -help
opt = parse_args(OptionParser(option_list = option_list, usage = "This script is general processing of single cell exp matrix!"))
print(opt)

###### step0. subset seurat object ####
timestart<-Sys.time() 

source("/home/zhengjh/scripts/cellchat/cellchat_func.r")

setwd(opt$output)

## load funtion
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
library(Seurat)
options(stringsAsFactors = FALSE)

# singledata_path <- list( zebrafish = "/home/zhengjh/projects/PinealGland/analysis/20221110_ZebrafishDataAnalysis/20221108_logtransform_harmony_Zebrafish_pineal_gland.rds",
# 									 rat = "/home/zhengjh/projects/PinealGland/analysis/20221108_RatDataAnalysis/20221108_logtransform_harmony_rat_day_pineal_gland.rds",
# 									 monkey = "/home/zhengjh/projects/PinealGland/analysis/20220602_MonkeyDataAnalysis/Seurat_Analysis/20220606_monkey_pineal_gland.rds")

#### load data ####
print("load data")

seurat_rds <- readRDS(opt$seurat_object_path)
CellChatDB <- readRDS(opt$lrdb_path)
DefaultAssay(seurat_rds) <- "RNA"
seurat_rds <- NormalizeData(seurat_rds)

future::plan("multiprocess", workers = 4) # do parallel

##### Part I: Data input & processing and initialization of CellChat object ####
print("Part I: Data input & processing and initialization of CellChat object")

cellchat <- CreateCellChatObject(seurat_rds = seurat_rds,  CellChatDB = CellChatDB, thresh_p = opt$thresh_p)
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


##### Part II: Inference of cell-cell communication network #####

print("Part II: Inference of cell-cell communication network")
### Compute the communication probability and infer cellular communication network ####
cellchat <- AggregateNetWrap(cellchat, thresh_p = opt$thresh_p)
saveRDS(cellchat, file = paste0("cellchat", opt$filelabel, sep = "_", ".rds"))
opar <- par(no.readonly = TRUE)

##### Part III: Visualization of cell-cell communication network ######
print("Part III: Visualization of cell-cell communication network")

print("Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram")
pathways.show <- cellchat@netP$pathways
VisualizeSignalPathway(cellchat, pathways.show)

print("Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways")
VisualizeLRs(cellchat = cellchat, folder_name = "2.VisualizationCellCellCommunicationNetwork", pathways.show = pathways.show)

###### Part IV: Systems analysis of cell-cell communication network ######

#### Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling ###
# print("Part IV: Systems analysis of cell-cell communication network")
# cellchat <- SystemsAnalysisCellCellNetwork(cellchat = cellchat, nPatterns = 3, pathways.show = pathways.show)


print("Happy~Finished!!")
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)




