#!/home/zhengjh/miniconda3/envs/r403/bin/Rscript
# parameter

library(optparse)
library(getopt)

option_list <- list(
	make_option(c("-o", "--output"), type = "character", default = FALSE,
							action = "store", help = "This is the output directory."
	),
	make_option(c("--seurat_object_path"), type = "character", default = FALSE,
							action = "store", help = "This is the path of seuratobject."
	),
	make_option(c("--assay_use"), type = "character", default = "RNA",
							action = "store", help = "This is assay type you want to use[default:RNA]."
	),
	make_option(c("--slot_use"), type = "character", default = "data",
							action = "store", help = "This is a slot type you want to use[default:data]."
	),
	
	make_option(c("--cluster_use"), type = "character", default = "celltype_assign",
							action = "store", help = "This is vector name you want to use to do pseudocell in the meta.data[default:celltype_assign]."
	),
	
	make_option(c("--pseudocell_size"), type = "integer", default = 10,
							action = "store", help = "This is the pseudocell size you want to use[default:10]"
							
	),
	make_option(c("--filelabel"), type = "character", default = FALSE,
							action = "store", help = "This is the filelabel of saved file you want to use."
							
	)
)

# -help
opt = parse_args(OptionParser(option_list = option_list,
															usage = "This Script is used to generate pseudocell seurat object."))
print(opt)

setwd(opt$output)

## load funtion
source("/home/zhengjh/scripts/seurat/SeuratWrapperFunction.R")

## time start
timestart<-Sys.time()

exp_count.seurat <- readRDS(opt$seurat_object_path)
PseudoCell.seurat <- PseudoCell(exp_count.seurat, opt$assay_use, opt$slot_use, opt$cluster_use, opt$pseudocell_size)
saveRDS(PseudoCell.seurat, paste0(opt$filelabel, ".PseudoCell.seurat.rds"))
print("Happy~Finished!!!")
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)


