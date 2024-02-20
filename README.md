# long-read-RNAseq
PIPELINE for identifying and quantifying known and novel genes/isoforms in long-read RNA-seq data

[link to description with pictures](https://docs.google.com/document/d/1mj8DaMMQsriclH1m1FKiJzLcJz12DW3rtCmNr7nHnaA/edit#heading=h.dawyqbpfox7p)


## OVERVIEW
Nextflow pipeline for identifying and quantifying known and novel genes/isoforms in long-read RNA-seq data. Now it works with human data from Oxford Nanopore platforms (in the future - PacBio and mice). 

**alignment** [Minimap2](https://github.com/lh3/minimap2)  mapping Oxford Nanopore reads to the genome

**cleaning** [Transcriptclean](https://github.com/dewyman/TranscriptClean)  corrects mismatches, microindels, and noncanonical splice junctions in long reads that have been mapped to the genome (to fix artifactual noncanonical splice junctions) 

**gene / isoform searching**  [TALON](https://github.com/dewyman/TALON)  identifying and quantifying known and novel genes / isoforms in long-read transcriptome data sets

![image text](https://github.com/vierstralab/long-read-RNAseq/blob/develop/Screen%20Shot%202024-02-12%20at%2010.15.20%20AM.png)

## INPUTS
In file `params.config`:

**Variables to be changed** 

* `samples_file`  - csv-file with row_id (number),  IDs of sample  ( we use LN  as ID) and full pathway to every file with reads for every sample.
Header : row_id,ln,pathway
1 row - 1 pathway to read file 
Example - _/net/seq/data2/projects/amuravyova/nf-long-reads-align/FETAL/11_20_fetal_with_pathways.csv_

* `outdir` - directory where you want to put the results
* `description` - description of the date (does not affect the analysis)
* `platform` - platform that was used for generating the data  (does not affect the analysis)

**Variables could to be changed (please don’t touch them now)**

* `genome_fasta` - fasta file containing the reference genome used in mapping
* `genome_gtf` - gtf file containing the reference annotation 
* `spl_jnk` - high-confidence splice junction file This file is necessary if you want to correct noncanonical splice junctions
* `known_variants_vcf` - vcf file containing variants
* `conda`

## HOW TO RUN
1. create `samples_file`
2. set the required variable values (`samples_file`, `outdir`, `description`, `platform`)  in file _/net/seq/data2/projects/amuravyova/nf-long-reads-align/long-read-RNAseq/**params.config**_ 
> [!CAUTION]
> save file changes !
3. run in tmux : 
```
module load nextflow/22.04.3 
nextflow run test_tuples.nf  -profile Altius -entry tuple
```
> [!IMPORTANT]
> please check that nobody else runs it now !

**Results** will be in the folder you set as `outdir` in file _params.config_

## OUTPUTS
![image text](https://github.com/vierstralab/long-read-RNAseq/blob/develop/Screen%20Shot%202024-02-14%20at%208.57.33%20AM.png)

[QC description](https://docs.google.com/document/d/1cqXZL64vZTV6GQu9x4Il9D0XDdYHMc8J0nMVN1uZ_p8/edit)




 

 

