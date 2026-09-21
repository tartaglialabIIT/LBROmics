#!/bin/bash

source /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/velo-venv/bin/activate

#Generate index
#You should have the files Mus_musculus.GRCm38.dna.primary_assembly.fa and Mus_musculus.GRCm38.94.gtf

#wget -q http://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
#wget -q http://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz

kb ref -i ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno Mus_musculus.GRCm38.dna.primary_assembly.fa.gz Mus_musculus.GRCm38.98.gtf.gz


#Generate count matrix from fastq
kb count --h5ad -i index.idx -g t2g.txt -x 10xv3 -o /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/RNA_velocity/DR1_WT/ -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno --filter bustools -t 16 --verbose /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s1_DR1_GATAACCTGC/s1_DR1_GATAACCTGC_S1_L002_R1_001.fastq.gz /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s1_DR1_GATAACCTGC/s1_DR1_GATAACCTGC_S1_L002_R2_001.fastq.gz /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s1_DR1_GATAACCTGC/s1_DR1_GATAACCTGC_S21_L004_R1_001.fastq.gz /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s1_DR1_GATAACCTGC/s1_DR1_GATAACCTGC_S21_L004_R2_001.fastq.gz /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s1_DR1_GATAACCTGC/s1_DR1_GATAACCTGC_S31_L003_R1_001.fastq.gz /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s1_DR1_GATAACCTGC/s1_DR1_GATAACCTGC_S31_L003_R2_001.fastq.gz 


kb count --h5ad -i index.idx -g t2g.txt -x 10xv3 -o /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/RNA_velocity/A8_mutant/ -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno --filter bustools -t 16 --verbose /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s2_A8_TTCACACCTT/s2_A8_TTCACACCTT_S2_L002_R1_001.fastq.gz /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s2_A8_TTCACACCTT/s2_A8_TTCACACCTT_S2_L002_R2_001.fastq.gz /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s2_A8_TTCACACCTT/s2_A8_TTCACACCTT_S32_L003_R1_001.fastq.gz /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s2_A8_TTCACACCTT/s2_A8_TTCACACCTT_S32_L003_R2_001.fastq.gz 
 

