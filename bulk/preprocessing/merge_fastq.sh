#!/bin/bash

cd ${LBROMICS_BULK_DATA:-.}/raw_data/a_3/cutadapt_out/

cat trimmed_a_3_HHHLCDSX2_L2_1.fq.gz trimmed_a_3_HHHLHDSX2_L2_1.fq.gz > trimmed_a_3_1.fq.gz
cat trimmed_a_3_HHHLCDSX2_L2_2.fq.gz trimmed_a_3_HHHLHDSX2_L2_2.fq.gz > trimmed_a_3_2.fq.gz

cd ${LBROMICS_BULK_DATA:-.}/raw_data/a_5/cutadapt_out/

cat trimmed_a_5_HHHLHDSX2_L2_1.fq.gz trimmed_a_5_HHJMHDSX2_L3_1.fq.gz > trimmed_a_5_1.fq.gz
cat trimmed_a_5_HHHLHDSX2_L2_2.fq.gz trimmed_a_5_HHJMHDSX2_L3_2.fq.gz > trimmed_a_5_2.fq.gz

cd ${LBROMICS_BULK_DATA:-.}/raw_data/a_9/cutadapt_out/

cat trimmed_a_9_HHHLHDSX2_L2_1.fq.gz trimmed_a_9_HHJMHDSX2_L4_1.fq.gz > trimmed_a_9_1.fq.gz
cat trimmed_a_9_HHHLHDSX2_L2_2.fq.gz trimmed_a_9_HHJMHDSX2_L4_2.fq.gz > trimmed_a_9_2.fq.gz
