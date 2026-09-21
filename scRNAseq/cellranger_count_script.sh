#!/bin/bash

export PATH=${LBROMICS_DATA_ROOT:-.}/cellranger/cellranger-6.1.2:$PATH

cd ${LBROMICS_DATA_ROOT:-.}/Cerase_Data/scRNA_seq/cellranger_count/

cellranger count --id s1_DR1_GATAACCTGC  --fastqs ${LBROMICS_DATA_ROOT:-.}/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s1_DR1_GATAACCTGC --transcriptome ${LBROMICS_DATA_ROOT:-.}/cellranger/refdata-gex-mm10-2020-A

cellranger count --id s2_A8_TTCACACCTT --fastqs ${LBROMICS_DATA_ROOT:-.}/Cerase_Data/scRNA_seq/X201SC19110138-Z01-F029/raw_data/s2_A8_TTCACACCTT --transcriptome ${LBROMICS_DATA_ROOT:-.}/cellranger/refdata-gex-mm10-2020-A
