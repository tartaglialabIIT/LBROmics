#!/bin/bash

samples=('a_2' 'a_3' 'a_4' 'a_5' 'a_6' 'a_7' 'a_8' 'a_9' 'a_10' 'a_11' 'a_12')
#samples=('a_1')

for s in "${samples[@]}";
   do
   mkdir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/bowtie2/
   python2.7 rrna_remove.py trimmed_${s}_ /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/bowtie2_analysis/bowtie_index/mouse_rrnas /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/trimmed_${s}_1.fq.gz /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/trimmed_${s}_2.fq.gz /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/bowtie2/
   gzip /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/bowtie2/trimmed_${s}__1_unmapped_sorted_stillpaired.fastq
   gzip /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/bowtie2/trimmed_${s}__2_unmapped_sorted_stillpaired.fastq
   /home/jonny/miniconda3/bin/STAR --readFilesCommand zcat --genomeDir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/GRCm38.98_genomeIndexes --peOverlapNbasesMin 10  --runThreadN 12 --readFilesIn /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/bowtie2/trimmed_${s}__1_unmapped_sorted_stillpaired.fastq.gz /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/bowtie2/trimmed_${s}__2_unmapped_sorted_stillpaired.fastq.gz --outFileNamePrefix trimmed_${s} --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted SortedByCoordinate
done
