#!/bin/bash

samples=('a_2' 'a_3' 'a_4' 'a_5' 'a_6' 'a_7' 'a_8' 'a_9' 'a_10' 'a_11' 'a_12')
#samples=('a_1')

for s in "${samples[@]}";
   do
   mkdir ${LBROMICS_BULK_DATA:-.}/raw_data/${s}/cutadapt_out/bowtie2/
   python2.7 rrna_remove.py trimmed_${s}_ ${LBROMICS_BULK_DATA:-.}/bowtie2_analysis/bowtie_index/mouse_rrnas ${LBROMICS_BULK_DATA:-.}/raw_data/${s}/cutadapt_out/trimmed_${s}_1.fq.gz ${LBROMICS_BULK_DATA:-.}/raw_data/${s}/cutadapt_out/trimmed_${s}_2.fq.gz ${LBROMICS_BULK_DATA:-.}/raw_data/${s}/cutadapt_out/bowtie2/
   gzip ${LBROMICS_BULK_DATA:-.}/raw_data/${s}/cutadapt_out/bowtie2/trimmed_${s}__1_unmapped_sorted_stillpaired.fastq
   gzip ${LBROMICS_BULK_DATA:-.}/raw_data/${s}/cutadapt_out/bowtie2/trimmed_${s}__2_unmapped_sorted_stillpaired.fastq
   STAR --readFilesCommand zcat --genomeDir ${LBROMICS_BULK_DATA:-.}/GRCm38.98_genomeIndexes --peOverlapNbasesMin 10  --runThreadN 12 --readFilesIn ${LBROMICS_BULK_DATA:-.}/raw_data/${s}/cutadapt_out/bowtie2/trimmed_${s}__1_unmapped_sorted_stillpaired.fastq.gz ${LBROMICS_BULK_DATA:-.}/raw_data/${s}/cutadapt_out/bowtie2/trimmed_${s}__2_unmapped_sorted_stillpaired.fastq.gz --outFileNamePrefix trimmed_${s} --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted SortedByCoordinate
done
