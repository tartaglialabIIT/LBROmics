#!/bin/bash

samples=('a_1' 'a_2' 'a_3' 'a_4' 'a_5' 'a_6' 'a_7' 'a_8' 'a_9' 'a_10' 'a_11' 'a_12')

for s in "${samples[@]}";
   do
   /home/jonny/miniconda3/bin/samtools index /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/bowtie2_analysis/full_scripts/sorted_bam/trimmed_${s}Aligned.sortedByCoord.out.bam
done
