##  MultiDataLogNormalizeCCAintegrate.R: use lognorm to normalize data ####

echo 'outputpath="/home/zhengjh/projects/PinealGland/analysis/IntegratedAnalysis/SeuratAnalysis"
project_name="Pinealgland_threespecies_integrated_analysis" # is the same in the MultipleDataQualityOverview.R file_prefix
inputdatapath=${outputpath}/20221115_PinealGland_Integrated_Aggregated_seurat.rds
cd ${outputpath}
file_prefix="Pinealgland_threespecies_integrated_lognorm" '  > run_threespecies_pinealgland_scRNA_seq_lognormccaintegrate.sh


echo ' Rscript /home/zhengjh/scripts/seurat/pipeline/MultiDataLogNormalizeCCAintegrate.R -o $outputpath --seuratobject $inputdatapath \
--min_gene 200 --max_gene 6500 --min_ncountRNA 0 --max_ncountRNA 25000 \
--max_mt_percent 30 --max_rb_percent 20 --file_prefix $file_prefix \
--scale_factor 10000 --nfeatures 3000 \
--nfeatures 3000  --npc_used 30 \
--k_parameter 25 --resolution_number 0.8 --min_pct 0.25 --logfc_threshold 0.25  \
--p_val_adj 0.05 --top_num 10 --mt_rb_pattern "^Mt-|^Rpl|^Rps" ' >> run_threespecies_pinealgland_scRNA_seq_lognormccaintegrate.sh

nohup sh run_threespecies_pinealgland_scRNA_seq_lognormccaintegrate.sh > step2.lognormalize.cca.output.txt 2>&1 &

## MultiDataLogNormalizeCCAintegrate.R: use lognorm to normalize data ####
echo 'outputpath="/home/zhengjh/projects/PinealGland/analysis/IntegratedAnalysis/SeuratAnalysis"
project_name="Pinealgland_threespecies_integrated_analysis" # is the same in the MultipleDataQualityOverview.R file_prefix
inputdatapath=${outputpath}/20221115_PinealGland_Integrated_Aggregated_seurat.rds
cd ${outputpath}
file_prefix="Pinealgland_threespecies_integrated_sctnorm" '  > run_threespecies_pinealgland_scRNA_seq_sctnormccaintegrate.sh


echo ' Rscript /home/zhengjh/scripts/seurat/pipeline/MultiDataSCTransformCCAintegrate.R -o $outputpath --seuratobject $inputdatapath \
--min_gene 200 --max_gene 6500 --min_ncountRNA 0 --max_ncountRNA 25000 \
--max_mt_percent 30 --max_rb_percent 20 --file_prefix $file_prefix \
--nfeatures 3000 \
--nfeatures 3000  --npc_used 30 \
--k_parameter 25 --resolution_number 0.8 --min_pct 0.25 --logfc_threshold 0.25  \
--p_val_adj 0.05 --top_num 10 --mt_rb_pattern "^Mt-|^Rpl|^Rps" ' >> run_threespecies_pinealgland_scRNA_seq_sctnormccaintegrate.sh

nohup sh run_threespecies_pinealgland_scRNA_seq_sctnormccaintegrate.sh > step2.sctnormalize.cca.output.txt 2>&1 &