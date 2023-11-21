#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process talon {
  publishDir "${params.outdir}"
  input:
//    tuple  
      path sample_trim_sam_onlyMD            
  output:
    path "${params.sample}/${params.sample}.db"
    path "${params.sample}/${params.sample}_talon_abundance_filtered.tsv"
    path "${params.sample}/${params.sample}_talon.gtf"
  script:

  sample="${params.sample}"
  sample_description="${params.description}"
  platform="${params.platform}"
  ref_fasta="${params.genome_fasta}"
  cov=0.9
  build="human_GRCh38_no_alt"
  annotation_v="human_gencode.v29"  
  db="${sample}.db"

  """
//  module load python/3.8
  mkdir ${sample}
  cd ${sample}
  mkdir labeled  #check exist  
  echo ${sample}
  
  talon_initialize_database \
     --f "${params.genome_gtf}" \
     --a "${annotation_v}" \
     --g "${build}" \
     --o "${sample}"

  talon_label_reads --f ../"${sample_trim_sam_onlyMD}"\
    --g "${ref_fasta}" \
    --t 10 \
    --ar 20 \
    --deleteTmp \
    --o labeled/"${sample}"

  module load bedtools
  touch config_"${sample}".csv
  echo "${sample},${sample_description},${platform},labeled/${sample}_labeled.sam"  > config_"${sample}".csv
  
  talon \
       --f config_"${sample}".csv\
       --db "${db}"  \
       --build "${build}"  \
       --cov "${cov}" \
       --o "${sample}" \
       --t 40   

  touch  dataset_"${sample}".txt 
  echo "${sample}"  > dataset_"${sample}".txt 

  talon_abundance \
       --db "${db}" \
       -a "${annotation_v}" \
       --build "${build}"  \
       -d dataset_"${sample}".txt  \
       --o "${sample}"
  
  talon_filter_transcripts \
       --db "${db}" \
       --datasets "${sample}" \
       -a "${annotation_v}" \
       --maxFracA 0.5 \
       --minCount 5 \
       --o "${sample}"_filtered_transcripts.csv

  talon_abundance \
       --db "${db}" \
       --whitelist "${sample}"_filtered_transcripts.csv \
       -a "${annotation_v}" \
       --build "${build}" \
       -d dataset_"${sample}".txt  \
       --o "${sample}"
  
  talon_create_GTF \
       --db "${db}"  \
       --whitelist "${sample}"_filtered_transcripts.csv \
       -a "${annotation_v}" \
       --build "${build}" \
       -d dataset_"${sample}".txt \
       --o "${sample}" 
  
  """
}

workflow {
  Channel.fromPath( '/net/seq/data2/projects/amuravyova/nf-long-reads-align/FETAL/wo_soft_clipped/23_03_27_PAK79776_barcode03_wo_soft_clipped_noQ_onlyMD.sam' ) 
  | talon
}
