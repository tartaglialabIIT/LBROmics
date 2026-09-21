#!/bin/bash

export PATH=/mnt/large/jfiorentino/cellranger/cellranger-6.1.2:$PATH

cd /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/cellranger_count/

cellranger count —-id s1_DR1_GATAACCTGC  --fastqs /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s1_DR1_GATAACCTGC --transcriptome /mnt/large/jfiorentino/cellranger/refdata-gex-mm10-2020-A

cellranger count —-id s2_A8_TTCACACCTT --fastqs /mnt/large/jfiorentino/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s2_A8_TTCACACCTT --transcriptome /mnt/large/jfiorentino/cellranger/refdata-gex-mm10-2020-A
