source activate r403

## 1) MultipleDataQualityOverview.R
outputpath=/home/zhengjh/projects/PinealGland/analysis/Zebrafish
cd $outputpath
inputdatapath=/home/zhengjh/projects/PinealGland/data/Zebrafish_GSE123778_RAW
file_prefix="zebrafish_scRNAseq_analysis"
cd /home/zhengjh/projects/PinealGland/analysis/Zebrafish
nohup Rscript /home/zhengjh/scripts/seurat/pipeline/MultipleDataQualityOverview.R -o $outputpath --datapath $inputdatapath --mt_pattern "^mt-" \
--rb_pattern "^rpl|^rps" --min_cell 1 --min_gene 5 --file_prefix $file_prefix > step1_MultipleDataQuality.output.txt &


##2 ) MutilDataLogHarmonyintegrate.R harmony---use Log clustering

cd /home/zhengjh/projects/PinealGland/analysis/Zebrafish
echo '
outputpath=/home/zhengjh/projects/PinealGland/analysis/Zebrafish
project_name="zebrafish_scRNAseq_analysis" # is the same in the MultipleDataQualityOverview.R file_prefix
inputdatapath=${outputpath}/${project_name}_multiple_dataquality_overview_seurat_rm_potencial_pituitary_cell.rds
file_prefix="zebrafish_scRNAseq_logtransform"
cd ${outputpath} ' > run_zebrafish_pinealgland_scRNA_seq_LogNormalizeHarmonyintegrate.sh

echo ' Rscript /home/zhengjh/scripts/seurat/pipeline/MutilDataLogHarmonyintegrate.R -o $outputpath --seuratobject $inputdatapath \
--min_gene 200 --max_gene 6500 --min_ncountRNA 0 --max_ncountRNA 50000 \
--max_mt_percent 30 --max_rb_percent 20 --file_prefix $file_prefix \
--scale_factor 10000 --nfeatures 3000  --npc_used 20 \
--k_parameter 20 --resolution_number 0.5 --min_pct 0.25 --logfc_threshold 0.25  \
--p_val_adj 0.05 --top_num 10 --mt_rb_pattern "^mt-|^Rpl|^Rps" ' >>  run_zebrafish_pinealgland_scRNA_seq_LogNormalizeHarmonyintegrate.sh


nohup sh run_zebrafish_pinealgland_scRNA_seq_LogNormalizeHarmonyintegrate.sh > step2.lognormalize.harmony.output.txt &

