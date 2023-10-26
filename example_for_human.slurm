#!/bin/bash
#SBATCH --partition=queue0,hpcz-1,pool
#SBATCH --job-name=am_talon
#SBATCH	--mem=100GB
#module load python/3.8

source ~/.bashrc

sample="HUMAN"
sam_file_path="/net/seq/data2/projects/amuravyova/nf-long-reads-align/MICE/results/HUMAN.sam"
sample_description="fetal_sample"  # dRNA, cDNA, PCR-cDNA 
platform="PromethIon"
ref_fasta="/home/amuravyova/ref/human_GRCh38/GRCh38_no_alt_analysis_set_GCA_000001405.15_short_names_chr_for__SJ_chrALL.fasta"
cov=0.9
# db="/home/amuravyova/ref/mouse/mouse_mm10.db" # to do
build="human_GRCh38_no_alt"
annotation_v="human_gencode.v29"  #  Which annotation version to use. Will determine which annotation transcripts are considered known or novel relative to. Note: must be in the TALON database. 

mkdir $sample
cd $sample
mkdir labeled  #check exist

echo $sample

talon_initialize_database \
     --f /home/amuravyova/ref/human_GRCh38/gencode.v29.primary_assembly.annotation_UCSC_names.gtf \
     --a human_gencode.v29 \
     --g human_GRCh38_no_alt \
     --o $sample

db="$sample.db"


talon_label_reads --f $sam_file_path\
    --g $ref_fasta \
    --t 10 \
    --ar 20 \
    --deleteTmp \
    --o labeled/$sample

module load bedtools
touch config_$sample.csv
echo "$sample,$sample_description,$platform,labeled/${sample}_labeled.sam"  > config_$sample.csv



talon \
       --f config_$sample.csv\
       --db ${db}  \
       --build $build  \
       --cov $cov \
       --o $sample \
       --t 40   


touch  dataset_$sample.txt 
echo "$sample"  > dataset_$sample.txt 

talon_abundance \
       --db ${db} \
       -a $annotation_v \
       --build $build  \
       -d dataset_$sample.txt  \
       --o $sample


talon_filter_transcripts \
       --db ${db} \
       --datasets $sample \
       -a $annotation_v \
       --maxFracA 0.5 \
       --minCount 5 \
       --o ${sample}_filtered_transcripts.csv

talon_abundance \
       --db ${db} \
       --whitelist ${sample}_filtered_transcripts.csv \
       -a $annotation_v \
       --build $build \
       -d dataset_$sample.txt  \
       --o $sample

talon_create_GTF \
       --db ${db}  \
       --whitelist ${sample}_filtered_transcripts.csv \
       -a $annotation_v \
       --build $build \
       -d dataset_$sample.txt \
       --o $sample  
