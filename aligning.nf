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
    // publishDir "${params.outdir}/minimap_logs", pattern: "log.*"
    tag "${cell_line}:${name}"
    cpus 5
//    scratch true

    input:
        tuple val(cell_line), path(fastq)
    output:
        tuple val(cell_line), path(name)
//      path log_files
    
    script:
    name = "${fastq.simpleName}.sam"
    log_files = "log.${fastq.simpleName}.txt"
    """
    minimap2 ${params.genome_fasta} ${fastq} -ax splice --MD -t ${task.cpus} > ${name}     
    """
}

process sam2bam_sort_indx {
    tag "${cell_line}:${name}"
    conda "${params.conda}"

    input:
        tuple val(cell_line), path(name)
    output:
        tuple val(cell_line), path(name_bam)

    script:
    name_sam = "${name.simpleName}.sam"
    name_bam = "${name.simpleName}.bam"
    """
    samtools view  -b  ${name_sam} >  ${name.simpleName}_unsorted.bam
    samtools sort  -O bam ${name.simpleName}_unsorted.bam > ${name_bam}
    samtools index ${name_bam}
    """
}

process soft_clip_trimming {
//   publishDir "${params.dir_result}"
   tag "${cell_line}:${name}"
   conda "${params.conda}"
   input:
     tuple val(cell_line), path(sample)
//     path dir_result
   output:
     tuple val(cell_line), path(sample_trim_sam)
   script: 
   sample_trim_sam="${sample.simpleName}_wo_soft_clipped_noQ.sam"
   sample_trim_bam="${sample.simpleName}_wo_soft_clipped_noQ.bam"
   """
   echo ${sample}
   echo ${sample_trim_bam}
   samtools index ${sample}
   python3 /net/seq/data2/projects/amuravyova/nf-long-reads-align/remove_soft_clipping_part_v2_modified.py ${sample} ${sample_trim_bam}
   samtools view -h ${sample_trim_bam} > ${sample_trim_sam}
   """
}

process take_only_MD {
//   publishDir "${params.dir_result}"
   conda "${params.conda}"
   input:
     tuple val(cell_line), path(sample_trim_sam)
   output:
     tuple val(cell_line), path(sample_trim_sam_onlyMD)
   script: 
   sample_trim_sam_onlyMD = "${sample_trim_sam.simpleName}_onlyMD.sam"
   """
   samtools view -H ${sample_trim_sam}  >> ${sample_trim_sam_onlyMD} 
   samtools view ${sample_trim_sam} | grep  "MD:"  >> ${sample_trim_sam_onlyMD}  
   """
}


process transcriptclean {
    container "/home/amuravyova/nextflow/transcriptclean_v2.0.2_cv1.sif"  //   container "biocontainers/transcriptclean:v2.0.2_cv1"
    containerOptions "${get_container(params.genome_fasta)} ${get_container(params.spl_jnk)} ${get_container(params.known_variants_vcf)}"
    tag "${cell_line}:${name}"
    cpus 5 

    input:
//        tuple val(cell_line), path(name)
        tuple val(cell_line), path(sample_trim_sam_onlyMD)
    output:
//        tuple val(cell_line), path(name_sam_clean)
        tuple val(cell_line), path(sample_trim_sam_onlyMD_clean)

    script:
//    name_bam = "${name.simpleName}.bam"
//    name_sam_sort = "${name.simpleName}_sorted.sam"
    sample_trim_sam_onlyMD_clean = "${sample_trim_sam_onlyMD.simpleName}_clean.sam"
    """ 
    head ${params.spl_jnk}
    head ${params.known_variants_vcf}
    TranscriptClean -s ${sample_trim_sam_onlyMD} -g ${params.genome_fasta}  --spliceJns ${params.spl_jnk} --variants ${params.known_variants_vcf} -t ${task.cpus} -o "${sample_trim_sam_onlyMD.simpleName}"  --primaryOnly  
    """
}

process merge_files {
    scratch true
    publishDir "${params.outdir}"
    tag "${cell_line}"
    conda params.conda

    input:
        tuple val(cell_line), path(bams)
    output:
        tuple val(cell_line), path(name), path("${name}.bai"), path(name_sam)
    
    script:
    name = "${cell_line}.bam"
    name_sam = "${cell_line}.sam"
    """
    samtools merge  -p -u merged.bam ${bams}
    samtools sort -O bam -o ${name} merged.bam
    samtools index ${name}
    samtools view -h ${name} > ${name_sam} 
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

workflow {
    fasta_ch = Channel.fromPath(params.samples_file)
        .splitCsv(header:true, sep:'\t')
        .map(row -> tuple(row.sample, file(row.fastq_file)))
    sam_files = fasta_ch
        | set_key_for_group_tuple
        | align_reads
        | sam2bam_sort_indx
        | soft_clip_trimming
        | take_only_MD 
        | transcriptclean
        | groupTuple() 
        | merge_files
//    bam_files = merge_files(al_reads.groupTuple())
}


