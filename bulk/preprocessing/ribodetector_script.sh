#!/bin/bash

samples=('a_1' 'a_2' 'a_3' 'a_4' 'a_5' 'a_6' 'a_7' 'a_8' 'a_9' 'a_10' 'a_11' 'a_12')

conda activate ribodetector

for s in "${samples[@]}";
   do
   mkdir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/ribodetector/
   mkdir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/ribodetector/nonrrna/
   mkdir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/ribodetector/rrna/
   ribodetector_cpu -t 20 \
     -l 150 \
     -i /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/trimmed_${s}_1.fq.gz /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/trimmed_${s}_2.fq.gz \
     -e rrna \
     --chunk_size 256 \
     -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/ribodetector/nonrrna/trimmed_${s}.nonrrna.1.fq.gz /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/ribodetector/nonrrna/trimmed_${s}.nonrrna.2.fq.gz \
     -r /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/ribodetector/rrna/trimmed_${s}.rrna.1.fq.gz /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/ribodetector/rrna/trimmed_${s}.rrna.2.fq.gz

done


