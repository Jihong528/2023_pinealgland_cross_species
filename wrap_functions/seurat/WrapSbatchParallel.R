#!/home/zhengjh/miniconda3/envs/r403/bin/Rscript
# parameter

library(optparse)
library(getopt)

option_list <- list(
  make_option(c( "--cl_cores"), type = "integer", default = FALSE,
              action = "store", help = "This is cluster cores you want to used."
  ),
  make_option(c("--code_path"), type = "character", default = FALSE,
              action = "store", help = "This is the path of code to run.")
)

# -help 
opt = parse_args(OptionParser(option_list = option_list, 
                              usage = "This Script is used to Sbatch Parallel!"))
print(opt)

timestart<-Sys.time() 

source("/home/zhengjh/scripts/seurat/function/SeuratWrapperFunction.R")

SbatchParallel(opt$cl_cores, opt$code_path)
print("Happy~Finished!!!")
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)

