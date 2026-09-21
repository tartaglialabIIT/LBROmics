#!/bin/bash


mkdir STAR_prova_save_unmapped
cd STAR_prova_save_unmapped/
/home/jonny/miniconda3/bin/STAR --readFilesCommand zcat --genomeDir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/GRCm38.98_genomeIndexes  --runThreadN 12 --readFilesIn /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/trimmed_a_1_1.fq.gz /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/trimmed_a_1_2.fq.gz --outFileNamePrefix trimmed_a1 --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted SortedByCoordinate --outReadsUnmapped Fastx
