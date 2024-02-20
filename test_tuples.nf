#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


def set_key_for_group_tuple(ch) {
  ch.groupTuple()
  .map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
  .transpose()
}


process align_reads {
    container "/home/amuravyova/nextflow/minimap2_v2.15dfsg-1-deb_cv1.sif"  //   container "docker://biocontainers/minimap2:v2.15dfsg-1-deb_cv1"
    containerOptions "${get_container(params.genome_fasta)}"
    publishDir "${params.outdir}/bams/${ln}"
    tag "${ln}:${fastq}"
    cpus 5
//    scratch true

    input:
       tuple val(row_id), val(ln), path(fastq)
    output:
//       tuple val(row_id), val(ln), path(name_sam), path(log_file)
       tuple val(row_id), val(ln), path(name_sam)    
    script:
    name_sam = "${ln}_${row_id}.sam"
    log_file = "${ln}_${row_id}_minimap2.log"
    publishDir = "${params.outdir}/bams/${ln}"
    """
    mkdir -p ${publishDir}
    minimap2 ${params.genome_fasta} ${fastq} -ax splice  --MD -t ${task.cpus} --secondary=no > ${name_sam} 2> ${publishDir}/${log_file}  
    """    
}

process soft_clip_trimming {
   publishDir "${params.outdir}/bams/${ln}"
   tag "${ln}:${name_sam}"
   conda "${params.conda}"
   input:
     tuple val(row_id), val(ln), path(name_sam)
   output:
     tuple val(row_id), val(ln), path(sample_trim_sam), path(sample_trim_QC), path(sample_skipped_reads)
   script: 
   name_bam = "${name_sam.simpleName}_sorted.bam"
   sample_trim_sam="${name_sam.simpleName}_trimmed.sam"
   sample_trim_bam="${name_sam.simpleName}_trimmed.bam"
   sample_trim_QC="${name_sam.simpleName}_trim_QC.tsv"
   sample_skipped_reads="${name_sam.simpleName}_skipped_reads"
   """
   samtools view  -b  ${name_sam} >  ${name_sam}_unsorted.bam
   samtools sort  -O bam ${name_sam}_unsorted.bam > ${name_bam}
   samtools index ${name_bam}
   python3 /net/seq/data2/projects/amuravyova/nf-long-reads-align/long-read-RNAseq/remove_soft_clipping_part_v2_modified.py ${name_bam} ${sample_trim_bam} ${sample_trim_QC} ${sample_skipped_reads}
   samtools view -h ${sample_trim_bam} > ${sample_trim_sam}
   """
}



process transcriptclean {
    container "/home/amuravyova/nextflow/transcriptclean_v2.0.2_cv1.sif"  //   container "biocontainers/transcriptclean:v2.0.2_cv1"
    containerOptions "${get_container(params.genome_fasta)} ${get_container(params.spl_jnk)} ${get_container(params.known_variants_vcf)}"
    publishDir "${params.outdir}/transcriptclean/${ln}"
    tag "${ln}:${sample_trim_sam}"
    cpus 5 

    input:
        tuple val(row_id), val(ln), path(sample_trim_sam)
//        tuple val(row_id), val(ln), path(sample_trim_sam_onlyMD)
    output:
        tuple val(row_id), val(ln), path(sample_trim_sam_clean), path(clean_log), path(clean_TE_log)

    script:
    sample_trim_sam_clean = "${ln}_${row_id}_clean.sam"
    clean_log = "${ln}_${row_id}_clean.log"
    clean_TE_log = "${ln}_${row_id}_clean.TE.log"
    """ 
    head ${params.spl_jnk}
    head ${params.known_variants_vcf}
    TranscriptClean -s ${sample_trim_sam} -g ${params.genome_fasta}  --spliceJns ${params.spl_jnk} --variants ${params.known_variants_vcf} -t ${task.cpus} -o "${ln}_${row_id}"  --primaryOnly  
    """
}


process take_only_MD {
   publishDir "${params.outdir}/bams/${ln}"
   conda "${params.conda}"
   input:
     tuple val(row_id), val(ln), path(sample_trim_sam_clean)
   output:
     tuple val(row_id), val(ln), path(sample_trim_sam_clean_onlyMD), path(sample_noMD_reads)
   script:
   sample_trim_sam_clean_onlyMD = "${sample_trim_sam_clean.simpleName}_onlyMD.sam"
   sample_noMD_reads = "${ln}_${row_id}_noMD.sam"
   """
   samtools view -H ${sample_trim_sam_clean}  >> ${sample_trim_sam_clean_onlyMD}
   samtools view ${sample_trim_sam_clean} | grep  "MD:"  >> ${sample_trim_sam_clean_onlyMD} 
   touch ${sample_noMD_reads}
   samtools view ${sample_trim_sam_clean} | grep -v "MD:" | wc -l >> ${sample_noMD_reads}
   if (( (samtools view ${sample_trim_sam_clean} | grep -v "MD:" | wc -l) > 0 )) ; then samtools view ${sample_trim_sam_clean} | grep -v "MD:" >> ${sample_noMD_reads} ; fi
   """
}


process merge_files {
    scratch true
    publishDir "${params.outdir}/merged_bams"
    tag "${ln}"
    conda params.conda

    input:
        tuple val(row_ids), val(ln), path(bams)
    output:
        tuple val(ln), path(name), path("${name}.bai"), path(name_sam)
    
    script:
    name = "${ln}_merged.bam"
    name_sam = "${ln}_merged.sam"
    """
    samtools merge  -p -u merged.bam ${bams}
    samtools sort -O bam -o ${name} merged.bam
    samtools index ${name}
    samtools view -h ${name} > ${name_sam} 
    """
}

process talon {
  publishDir "${params.outdir}/TALON"
  input:
    tuple  val(ln), path(name), path("${name}.bai"), path(name_sam)  
//      path sample_trim_sam_onlyMD            
  output:
    tuple path(db) , path("${ln}_talon_abundance_filtered.tsv"), path("${ln}_talon.gtf"), path("${ln}_QC.log")

  script:

  sample="${ln}"
  sample_description="${params.description}"
  platform="${params.platform}"
  ref_fasta="${params.genome_fasta}"
  cov=0.9
  build="human_GRCh38_no_alt"
  annotation_v="human_gencode.v29"  
  db="${ln}.db"

  """
//  module load python/3.8
//  mkdir ${sample}
//  cd ${sample}
  mkdir labeled  #check exist  
  echo ${sample}
  
  talon_initialize_database \
     --f "${params.genome_gtf}" \
     --a "${annotation_v}" \
     --g "${build}" \
     --o "${sample}"

  talon_label_reads --f "${name_sam}"\
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

def get_container(file_name) {
  parent = file(file_name).parent
  container = "--bind ${parent}"
  if (file(file_name).exists()) {
        old_parent = file(file_name).toRealPath().parent
        if (old_parent != parent) {
                container += ",${old_parent}"
        }
  } 
  return container
}

workflow tuple {
     metadata_ch  = Channel.fromPath(params.samples_file)
        .splitCsv(header:true, sep:',')
        .map(row -> tuple(row.row_id, row.ln, row.pathway))
        .groupTuple(by: 1)
       	.map{ it -> tuple(it[0],groupKey(it[1], it[2].size()), it[2]) }
        .transpose()
        .collate(4)
        .flatten()
        .collate(3)
        |align_reads
        |soft_clip_trimming
        |transcriptclean
        |take_only_MD
        |groupTuple(by: 1)
        |view()
        |merge_files
        |talon
}

workflow tuple2 {
     metadata_ch  = Channel.fromPath(params.samples_file)
        .splitCsv(header:true, sep:',')
        .map(row -> tuple(row.ln, (row.pathway)))
        .view()
    sam_files = metadata_ch
        .groupTuple(by:0)
        .view()
        .map{ it -> tuple(groupKey(it[0], it[1].size()), it[1]) }
        .view()
        .transpose()
        .view()
}

workflow old {
    fasta_ch = Channel.fromPath("/home/amuravyova/LR_RNAseq/RNA-Seq_LR/cDNA/SS_pipeline/2_test_fetal_samples_file_list_fastqs.tsv")
        .splitCsv(header:true, sep:'\t')
        .map(row -> tuple(row.sample, file(row.fastq_file)))
        .view()
    sam_files = fasta_ch
        .groupTuple()
        .view()
        .map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
        .view()
        .transpose()
        .view()
}
