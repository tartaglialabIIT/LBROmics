#!/bin/bash
wget http://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.98.gtf.gz

mkdir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/GRCm38.98_genomeIndexes
chmod 777 /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/GRCm38.98_genomeIndexes
/home/jonny/miniconda3/bin/STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/GRCm38.98_genomeIndexes --genomeFastaFiles ./Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile ./Mus_musculus.GRCm38.98.gtf --sjdbOverhang 149 --limitGenomeGenerateRAM 100000000000
