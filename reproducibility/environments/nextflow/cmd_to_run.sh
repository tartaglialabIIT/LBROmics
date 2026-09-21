#!/usr/bin/env bash
# nf-core/sammyseq launch used for differential solubility analysis
# Pin: see ../nextflow_pins.env (Nextflow 24.10.4; sammyseq commit fa6f6ffeb3)
# Adapt genome/blacklist/GTF paths to your system.
# Samplesheet: 4fSAMMYseq/pipeline_launch/samplesheet_18_11_2025.csv

nextflow run nf-core/sammyseq -r fa6f6ffeb3 \
  -profile docker \
  --fasta /path/to/mm10.fa \
  --input 4fSAMMYseq/pipeline_launch/samplesheet_18_11_2025.csv \
  --outdir results_sammyseq \
  --blacklist /path/to/mm10-blacklist.v2.bed \
  --keep_regions_bed /path/to/mm10_canonical_chromosomes.bed \
  --comparison S2SvsS3,S2LvsS3,S2SvsS4,S2LvsS4,S4vsS3 \
  --compare_groups MUT_NPCvsWT_NPC,WT_NPCvsWT_ESC,MUT_ESCvsWT_ESC \
  --differential_solubility true \
  --gtf /path/to/mm10.gtf
