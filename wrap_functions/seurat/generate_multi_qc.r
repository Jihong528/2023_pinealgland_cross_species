#!/home/zhengjh/miniconda3/envs/r403/bin/Rscript
# parameter

library(optparse)
library(getopt)

option_list <- list(
  make_option(c("-o", "--output"), type = "character", default = FALSE,
              action = "store", help = "This is the output directory."
  ),
  make_option(c("--datapath"), type = "character", default = FALSE,
              action = "store", help = "This is the path of cellranger out which each one folder is one sample."
  ),
  make_option(c("--mt_pattern"), type = "character", default = "^mt-",
              action = "store", help = "This is mitochondria pattern used to calculate its percentage, default is the ^mt-."
  ),
  make_option(c("--rb_pattern"), type = "character", default = "^Rpl|^Rps",
              action = "store", help = "This is ribosome pattern used to calculate its percentage, default is ^Rpl|^Rps"
  ),
  
  make_option(c("--min_cell"), type = "integer", default = 3,
              action = "store", help = "This is min cells number to filter, default is 3."
  ),
  make_option(c("--min_gene"), type = "integer", default = 0,
              action = "store", help = "This is min gene number to filter, default is 0."
  ),
  make_option(c("--file_prefix"), type = "character", default = FALSE,
              action = "store", help = "This is file_prefix of files saved."
  )
)

# -help 
opt = parse_args(OptionParser(option_list = option_list, 
                              usage = "This Script is general overview data quality of multiple single cell exp matrix!"))
print(opt)

#### step0.load data ####

source("/home/zhengjh/scripts/seurat/SeuratWrapperFunction.R")

print("Step0: Load Data")

##### step00. load functions ####
# read in the functions
timestart <- Sys.time() 
# set the output path
cat("Current working dir: ", opt$output, "\n")

setwd(opt$output)
# load ref object

#### step01. load exp matrix ####

print("Step0: Load exp matrix data from cellranger output, \n")
print("File used")
print(list.files(opt$datapath, recursive = T))
file_name <- list.files(opt$datapath)
print("Project name used in Seurat object")
print(file_name)

for (i in 1:length(file_name)){
  matrix.file <- paste(opt$datapath, file_name[i], sep = "/")
  print(matrix.file)
  temp <- Read10X(data.dir = matrix.file)
  newcolname <- paste(file_name[i], colnames(temp), sep = "_")
  colnames(temp) <- newcolname
  assign(file_name[i], temp)
  print(file_name[i])
  print(dim(get(file_name[i])))  
  rm(temp)
}
rm(i, matrix.file, newcolname)

#### step1: Create_seurat_filter_cell ######
object.list <- list()

print("Step1: Creat Seuart Object with filtering lower genes and lower cells")
for (i in 1:length(file_name)){
  temp <- CreateFirstFilterSeuratObject(exp_count = get(file_name[i]),
                                        pr_name = file_name[i],
                                        min_cell = opt$min_cell, 
                                        min_gene = opt$min_gene,
                                        mt_pattern = opt$mt_pattern,
                                        rb_pattern = opt$rb_pattern)
  object.list[[file_name[i]]] <- temp
  rm(temp)
}

#### step2: merge object #####
print("Step2: Megre each seurat object")
mergeobject <- merge(x= object.list[[1]], y = object.list[2:length(object.list)])

name_pre <- paste(opt$file_prefix, "filter_gene_expressed", opt$min_cell, "cells", opt$min_gene, "genes", sep = "_") 
mergeobject <- SeuratQC(mergeobject, name_pre)

## cellnumber
plot_data <- as.data.frame(table(mergeobject$orig.ident))
colnames(plot_data) <- c("group", 'Freq')
plot_data$group <- factor(plot_data$group, levels = file_name )
p <- ggplot(data =plot_data, aes(x = group, y = Freq, fill = group)) + geom_bar(stat="identity") 
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p <- p + geom_text(aes(label = Freq), position = position_dodge(0.9), vjust = -0.8)

pdfname <- paste0(opt$file_prefix, "_Cellnumber.pdf")
pdf(pdfname, 12, 6)
print(p)
dev.off()

folder_name <- "1.QC"
file.copy(pdfname, folder_name, overwrite = T)#拷贝文件
file.remove(pdfname)

rds_name <- paste0(opt$file_prefix, "_multiple_dataquality_overview_seurat.rds")
saveRDS(mergeobject, file = rds_name)

print("Happy~Finished!!!")
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)





