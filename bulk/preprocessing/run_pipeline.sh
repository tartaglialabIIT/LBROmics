#!/bin/bash

cd ${LBROMICS_BULK_DATA:-.}/

chmod +x run_cutadapt.sh
./run_cutadapt.sh

cd ${LBROMICS_BULK_DATA:-.}/

chmod +x run_fastqc.sh
./run_fastqc.sh

cd ${LBROMICS_BULK_DATA:-.}/

chmod +x merge_fastq.sh
./merge_fastq.sh

