#!/bin/bash

ln="$1"
fastq_folder=$(grep $ln /net/seq/data2/projects/amuravyova/nf-long-reads-align/FETAL/ALL_fetal_with_pathways.tsv | cut -f9)
fastq_fail_folder=$( echo "$fastq_folder" | sed "s/fastq_pass/fastq_fail/g")

folder="/net/seq/data2/projects/amuravyova/nf-long-reads-align/FETAL_new"
result="$folder/QC_NEW/"$ln"_QC.csv"



# how many reads we have
all_reads=$(for fastq in $fastq_folder/*.fastq.gz ; do zcat $fastq | wc -l ; done | awk '{sum += $1/4} END {printf "%.0f\n",sum}')
echo "All reads,"$all_reads >> $result

# how many failed reads we have
for fastq in $fastq_fail_folder/*.fastq.gz ; do zcat $fastq | wc -l ; done | awk '{sum += $1/4} END {print "Failed reads,",sum}' >> $result

# how many mapped reads (Log minimap2)
mapped_reads=$(for logfile in $folder/bams/$ln/*_minimap2.log ; do grep -Eo "mapped [0-9]+ sequences" $logfile | grep -oE '[0-9]+'; done | awk '{sum += $1} END {print sum}')
echo "mapped reads,"$mapped_reads >> $result

# how many Unmapped reads 
unmapped_reads=$((all_reads - mapped_reads))
echo "Unmapped reads,"$unmapped_reads >> $result

# how many reads don't contain MD-tag 
for noMD in $folder/bams/$ln/*_noMD.sam; do  grep -E '[0-9]+' $noMD ; done | awk '{sum += $1} END {print "noMD reads,",sum}' >> $result


# how many trimmed reads
for trim in $folder/bams/$ln/*_trim_QC.tsv; do wc -l $trim ; done | awk '{sum += $1} END {print "trimmed reads,",sum}' >> $result


# How many skipped_reads by trimming (if read.is_unmapped or read.query_sequence is None or read.query_qualities is None)
for skip in $folder/bams/$ln/*_skipped_reads; do wc -l $skip ; done | awk '{sum += $1} END {print "skipped_reads,",sum}' >> $result


# How many "points" were corrected and uncorrected by Transcriptclean
for clean in $folder/transcriptclean/$ln/*_clean.TE.log; do grep -c "Uncorrected" $clean ; done | awk '{sum += $1} END {print "Uncorrected points,",sum}' >> $result

# header contains “Corrected” - we need (sum-1) for every file
for clean in $folder/transcriptclean/$ln/*_clean.TE.log; do grep -c "Corrected" $clean | awk '{ print  $1 - 1 }' ; done | awk '{sum += $1} END {print "Corrected points,",sum}' >> $result


# How many reads passed and didn't pass TALON QC  ( first 6 lines is header)
# passed_QC
# 	3
# value in column:  1 - if QC passed, 0 - if QC failed
tail -n +7 "$folder/TALON/$ln"_QC.log | awk '{sum += $3} END {print "Reads passed TALON QC,",sum}' >> $result
tail -n +7 "$folder/TALON/$ln"_QC.log | awk '{sum += ($3 == 0)} END {print "Reads did NOT passed TALON QC,",sum}' >> $result


# Mean, SD and median (just for lenght) -  for  lenght, fraction, identity
# read_length     fraction_aligned        identity
# 	5 		6			7 

tail -n +7 "$folder/TALON/$ln"_QC.log | awk '{ sum += $5 } END { print "Mean read_length,",sum/NR }'  >> $result
tail -n +7 "$folder/TALON/$ln"_QC.log | awk '{ sum += $5; sumsq += ($5)^2 } END { mean = sum/NR; print "SD read_length,",sqrt(sumsq/NR - (mean)^2) }'  >> $result
tail -n +7 "$folder/TALON/$ln"_QC.log | awk '{print $5}' | sort -n | awk '{a[NR]=$1} END {if (NR%2) print "Median read_length,",a[int(NR/2)+1]; else print "Median read_length,",(a[NR/2]+a[NR/2+1])/2}'  >> $result

tail -n +7 "$folder/TALON/$ln"_QC.log | awk '{ sum += $6 } END { print "Mean fraction_aligned,",sum/NR }'  >> $result
tail -n +7 "$folder/TALON/$ln"_QC.log | awk '{ sum += $6; sumsq += ($6)^2 } END { mean = sum/NR; print "SD fraction_aligned,",sqrt(sumsq/NR - (mean)^2) }'  >> $result

tail -n +7 "$folder/TALON/$ln"_QC.log | awk '{ sum += $7 } END { print "Mean identity,",sum/NR }'  >> $result
tail -n +7 "$folder/TALON/$ln"_QC.log | awk '{ sum += $7; sumsq += ($7)^2 } END { mean = sum/NR; print "SD identity,",sqrt(sumsq/NR - (mean)^2) }'  >> $result



# How many genes and transcripts were found
# annot_gene_id		annot_transcript_id
# 	3			4

tail -n +2 "$folder/TALON/$ln"_talon_abundance_filtered.tsv | awk '{print $3}' | sort | uniq | wc -l | awk '{print "Genes found,",$1}' >> $result
tail -n +2 "$folder/TALON/$ln"_talon_abundance_filtered.tsv | awk '{print $4}' | sort | uniq | wc -l | awk '{print "Transcripts found,",$1}' >> $result

