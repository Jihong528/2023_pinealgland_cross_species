# Cross-species single-cell landscape of vertebrate pineal gland

This repository contains scripts to reproduce the results of the single-cell data from: "Cross-species single-cell landscape of vertebrate pineal gland"


The files contain code for the following analyses:

## pinealgland_clustering_analysis 
- **clustering_analysis.Rmd** --> main scripts: QC, preprocessing, clusteirng, cell type identification and batch correction (input data are raw count matrices); generating figure1,2 and related supplemental figures
- **20221103_ZebrafishPinealGlandPipelineAnalysis.sh** --> the script used to remove the batch effect of zebrafish data
- **20221106_RatPinealGlandPipelineAnalysis**  --> the script used to remove the batch effect of rat data
- **20221115_PinealGland_IntegrateThreeSpecies.sh** --> the script used to remove the batch effect of three species data
- **ortholog_gene_mapping.xlsx** --> the one-to-one orthologs tables

## phototransduction_melatonin_circadian_gene_analysis
- **phototransduction_melatonin_circadian_gene_analysis.Rmd** --> the script used to generate the figures related to phototransduction, melatonin & circadian gene expression
- **GeneSets.csv** --> the gene list of phototransduction, melatonin & circadian

## GPCR_analysis
- **GPCR_analysis.Rmd** --> the script used to generate the figures related to GPCR expression
- **GPCRList.xlsx** --> the GPCR list used

## TF_analysis   
- **TF_analysis.Rmd** --> the script used to generate the figures related to TF expression
- **TFList.xlsx** --> the TF list used

## CellChat_analysis
- **CellChat_analysis.Rmd** --> the script used to generate the figures related to cell-cell communication

## retina_analysis
- **retina_data_analysis.Rmd** --> the script used to process the retina data and generate the related figures

## wrap_functions
- **cellchat** --> the functions used in cellchat analysis
- **celltypesimilarity** -> the functions used to calculate the celltype similarity of two species
- **gpcr** --> the functions used to generate the figures in GPCR analysis
- **seurat** --> the functions used in Seurat analysis
- **tf** --> the functions used  to generate the figures in TF analysis


If you have any questions about the data or analysis feel free to contact us. :)

