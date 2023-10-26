#!/bin/bash
#SBATCH --partition=queue0,hpcz-1,pool
#SBATCH --job-name=am_talon
#SBATCH	--mem=100GB
#module load python/3.8

source ~/.bashrc

sample="MICE"
sam_file_path="/net/seq/data2/projects/amuravyova/nf-long-reads-align/MICE/results/MICE.sam"
sample_description="mice_sample"  # dRNA, cDNA, PCR-cDNA 
platform="PromethIon"
ref_fasta="/home/amuravyova/ref/mouse/mm10_no_alt_analysis_set_ENCODE.fasta"
cov=0.9
db="/home/amuravyova/ref/mouse/mouse_mm10.db"
build="mouse_mm10_no_alt"
annotation_v="mouse_gencode.vM21"  #  Which annotation version to use. Will determine which annotation transcripts are considered known or novel relative to. Note: must be in the TALON database. 

mkdir $sample
cd $sample
mkdir labeled  #check exist

echo $sample

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
