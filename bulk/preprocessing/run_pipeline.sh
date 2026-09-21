#!/bin/bash

cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/

chmod +x run_cutadapt.sh
./run_cutadapt.sh

cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/

chmod +x run_fastqc.sh
./run_fastqc.sh

cd /mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/

chmod +x merge_fastq.sh
./merge_fastq.sh

