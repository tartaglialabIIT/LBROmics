#!/bin/bash

mkdir STAR_paired_prova2
cd STAR_paired_prova2/
/home/jonny/miniconda3/bin/STAR --readFilesCommand zcat --genomeDir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/GRCm38.98_genomeIndexes  --runThreadN 12 --alignEndsProtrude 10 ConcordantPair --readFilesIn /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/cutadapt_out/trimmed_a_1_1.fq.gz /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/cutadapt_out/trimmed_a_1_2.fq.gz --outFileNamePrefix trimmed_a_1 --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted SortedByCoordinate

mkdir STAR_prova_notrim
cd STAR_prova_notrim/
/home/jonny/miniconda3/bin/STAR --readFilesCommand zcat --genomeDir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/GRCm38.98_genomeIndexes  --runThreadN 12 --readFilesIn /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/a_1_EKDL210004904-1a_HHJMHDSX2_L3_1.fq.gz /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/a_1_EKDL210004904-1a_HHJMHDSX2_L3_1.fq.gz --outFileNamePrefix a_1 --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted SortedByCoordinate
