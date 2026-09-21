#!/bin/bash

samples=('a_1' 'a_2' 'a_3' 'a_4' 'a_5' 'a_6' 'a_7' 'a_8' 'a_9' 'a_10' 'a_11' 'a_12')

for s in "${samples[@]}";
   do
   mkdir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/sortmerna/
   /home/jonny/miniconda3/bin/sortmerna -threads 25 --ref /database/smr_v4.3_default_db.fasta --reads /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/trimmed_${s}_1.fq.gz --reads /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/trimmed_{$s}_2.fq.gz --fastx --paired_out --out2 --aligned rRNA-reads --other non-rRNA-reads --workdir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/${s}/cutadapt_out/sortmerna/
done
