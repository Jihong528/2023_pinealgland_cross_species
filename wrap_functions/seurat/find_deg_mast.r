#!/home/zhengjh/miniconda3/envs/r403/bin/Rscript
# parameter

library(optparse)
library(getopt)

option_list <- list(
  make_option(c("-o", "--output"), type = "character", default = FALSE,
              action = "store", help = "This is the output directory."
  ),
  make_option(c("--seuratobject"), type = "character", default = FALSE,
              action = "store", help = "This is the path of exp matrix."
  ),
  make_option(c("--file_prefix"), type = "character", default = FALSE,
              action = "store", help = "This is number of file prefix of results."
  ),
  make_option(c("--cell_type_filepath"), type = "character", default = FALSE,
              action = "store", help = "This is a datapath of celltypes you want to find DEGs."
  ),
  
  make_option(c("--treat"), type = "character", default = FALSE,
              action = "store", help = "This is the name of your treat group."
  ),
  
  make_option(c("--control"), type = "character", default = FALSE,
              action = "store", help = "This is the name of your control group."
              
  ),
  make_option(c("--min_pct_cell"), type = "integer", default = 0.1,
              action = "store", help = "This is min_pct_cell used in MAST test."
              
  )
)

# -help 
opt = parse_args(OptionParser(option_list = option_list, 
                              usage = "This Script is used to find DEGs by MAST test!"))
print(opt)

setwd(opt$output)

## load funtion
source("/home/zhengjh/scripts/seurat/function/SeuratWrapperFunction.R")

## time start
timestart<-Sys.time() 

cell_type <- read.csv(opt$cell_type_filepath)
cell_type <- cell_type$celltype
exp_count.seurat <- readRDS(opt$seuratobject)
FindDEGsFromMasttest(exp_count.seurat = exp_count.seurat,
                     cell_type = cell_type, treat = opt$treat, control = opt$control, min_pct_cell = opt$min_pct_cell, file_prefix = opt$file_prefix)

print("Happy~Finished!!!")
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)

