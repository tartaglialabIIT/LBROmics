#!/bin/bash

cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/cutadapt_out/trimmed_a_1_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/cutadapt_out/trimmed_a_1_2.fq.gz *_1.fq.gz *_2.fq.gz

cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_2/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_2/cutadapt_out/trimmed_a_2_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_2/cutadapt_out/trimmed_a_2_2.fq.gz *_1.fq.gz *_2.fq.gz

cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_3/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_3/cutadapt_out/trimmed_a_3_HHHLCDSX2_L2_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_3/cutadapt_out/trimmed_a_3_HHHLCDSX2_L2_2.fq.gz a_3_EKDL210004906-1a_HHHLCDSX2_L2_1.fq.gz a_3_EKDL210004906-1a_HHHLCDSX2_L2_2.fq.gz
/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_3/cutadapt_out/trimmed_a_3_HHHLHDSX2_L2_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_3/cutadapt_out/trimmed_a_3_HHHLHDSX2_L2_2.fq.gz a_3_EKDL210004906-1a_HHHLHDSX2_L2_1.fq.gz a_3_EKDL210004906-1a_HHHLHDSX2_L2_2.fq.gz

cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_4/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_4/cutadapt_out/trimmed_a_4_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_4/cutadapt_out/trimmed_a_4_2.fq.gz *_1.fq.gz *_2.fq.gz

cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_5/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_5/cutadapt_out/trimmed_a_5_HHHLHDSX2_L2_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_5/cutadapt_out/trimmed_a_5_HHHLHDSX2_L2_2.fq.gz a_5_EKDL210004908-1a_HHHLHDSX2_L2_1.fq.gz a_5_EKDL210004908-1a_HHHLHDSX2_L2_2.fq.gz
/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_5/cutadapt_out/trimmed_a_5_HHJMHDSX2_L3_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_5/cutadapt_out/trimmed_a_5_HHJMHDSX2_L3_2.fq.gz a_5_EKDL210004908-1a_HHJMHDSX2_L3_1.fq.gz a_5_EKDL210004908-1a_HHJMHDSX2_L3_2.fq.gz


cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_6/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_6/cutadapt_out/trimmed_a_6_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_6/cutadapt_out/trimmed_a_6_2.fq.gz *_1.fq.gz *_2.fq.gz


cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_7/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_7/cutadapt_out/trimmed_a_7_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_7/cutadapt_out/trimmed_a_7_2.fq.gz *_1.fq.gz *_2.fq.gz


cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_8/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_8/cutadapt_out/trimmed_a_8_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_8/cutadapt_out/trimmed_a_8_2.fq.gz *_1.fq.gz *_2.fq.gz


cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_9/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_9/cutadapt_out/trimmed_a_9_HHHLHDSX2_L2_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_9/cutadapt_out/trimmed_a_9_HHHLHDSX2_L2_2.fq.gz a_9_EKDL210004912-1a_HHHLHDSX2_L2_1.fq.gz a_9_EKDL210004912-1a_HHHLHDSX2_L2_2.fq.gz
/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_9/cutadapt_out/trimmed_a_9_HHJMHDSX2_L4_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_9/cutadapt_out/trimmed_a_9_HHJMHDSX2_L4_2.fq.gz a_9_EKDL210004912-1a_HHJMHDSX2_L4_1.fq.gz a_9_EKDL210004912-1a_HHJMHDSX2_L4_2.fq.gz

cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_10/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_10/cutadapt_out/trimmed_a_10_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_10/cutadapt_out/trimmed_a_10_2.fq.gz *_1.fq.gz *_2.fq.gz

cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_11/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_11/cutadapt_out/trimmed_a_11_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_11/cutadapt_out/trimmed_a_11_2.fq.gz *_1.fq.gz *_2.fq.gz

cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_12/ 

mkdir cutadapt_out

/home/jonny/miniconda3/envs/cutadapt/bin/cutadapt -j 0 -a GATCGGAAGAGCACACGT -A GATCGGAAGAGCGTCGTG -m 18 -o /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_12/cutadapt_out/trimmed_a_12_1.fq.gz -p /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_12/cutadapt_out/trimmed_a_12_2.fq.gz *_1.fq.gz *_2.fq.gz


