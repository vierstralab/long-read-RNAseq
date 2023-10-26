#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// Script parameters

params.bams_dir="/net/seq/data2/projects/amuravyova/nf-long-reads-align/test_files/*.bam"
params.dir_result="/net/seq/data2/projects/amuravyova/nf-long-reads-align/test_files/result"
params.description="fetal_sample"
params.platform="PromethIon"
params.example_slurm="/net/seq/data2/projects/amuravyova/nf-long-reads-align/example_for_human.slurm"


process soft_clip_trimming {
   publishDir "${params.dir_result}"
   conda "${params.conda}"
   input:
     path sample
     path dir_result
   output:
     path sample_trim_sam
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
   publishDir "${params.dir_result}"
   conda "${params.conda}"
   input:
     path sample_trim_sam
   output:
     path sample_trim_sam_onlyMD
   script: 
   sample_trim_sam_onlyMD = "${sample_trim_sam.simpleName}_onlyMD.sam"
   """
   samtools view -H ${sample_trim_sam}  >> ${sample_trim_sam_onlyMD} 
   samtools view ${sample_trim_sam} | grep  "MD:"  >> ${sample_trim_sam_onlyMD}  
   """
}

process talon {
  publishDir "${params.dir_result}"
  input:
    path sample_trim_sam_onlyMD           
  output:
    path sample_slurm
  script:
  sample="${sample_trim_sam_onlyMD.simpleName}"
  sample_slurm = "${sample_trim_sam_onlyMD.simpleName}.slurm" 
  sample_path = "${params.dir_result}/${sample_trim_sam_onlyMD}"
  """  
  echo ${sample} ${sample_path} ${params.description} ${params.platform}
  cp ${params.example_slurm} ${sample_slurm}
  sed -i  's/sample=.*/sample="'${sample}'"/' ${sample_slurm}
  sed -i  's|sam_file_path=.*|sam_file_path="'${sample_path}'"|' ${sample_slurm}
  sed -i  's/sample_description=.*/sample_description="'${params.description}'"/' ${sample_slurm}
  sed -i  's/platform=.*/platform="'${params.platform}'"/' ${sample_slurm}
  sbatch ${sample_slurm} 
  """
}


workflow {
   def bams_ch = Channel.fromPath(params.bams_dir)
   bams_ch.view { "path: ${it}" }
   soft_clip_trimming(bams_ch, params.dir_result)|take_only_MD|talon
}

