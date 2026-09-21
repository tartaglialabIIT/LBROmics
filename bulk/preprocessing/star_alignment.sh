#!/bin/bash

samples=('a_1' 'a_2' 'a_3' 'a_4' 'a_5' 'a_6' 'a_7' 'a_8' 'a_9' 'a_10' 'a_11' 'a_12')

for s in "${samples[@]}";
   do
   /home/jonny/miniconda3/bin/STAR --genomeDir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/GRCm38.98_genomeIndexes  --runThreadN 25 --readFilesIn /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/trimmed_${s}_1.fq.gz /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/trimmed_{$s}_2.fq.gz --outFileNamePrefix trimmed_${s} --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted SortedByCoordinate
done
