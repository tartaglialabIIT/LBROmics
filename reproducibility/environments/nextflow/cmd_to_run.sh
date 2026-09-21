nextflow run nf-core/sammyseq -r dev \
    -profile docker \
    --fasta /mnt/storage/jonni/UCSC_mouse_genome/mm10.fa\
    --input ./samplesheet_11_09_2025.csv \
    --blacklist ./mm10-blacklist.v2.bed \
    --save_reference true \
    --outdir ./results_11_09_2025


nextflow run nf-core/sammyseq -r dev \
    -profile test,docker \
    --outdir ./results


nextflow run nf-core/sammyseq -r dev \
    -resume \
    -profile docker \
    --fasta /mnt/storage/jonni/UCSC_mouse_genome/mm10.fa\
    --input ./samplesheet_11_09_2025.csv \
    --blacklist ./mm10-blacklist.v2.bed \
    --save_reference true \
    --outdir ./results_11_09_2025/comparisons \
    --comparison S2SvsS3,S2LvsS3,S2SvsS4,S2LvsS4,S4vsS3 \
    --comparison_maker spp

nextflow run nf-core/sammyseq -r dev \
  -resume \
  -profile docker \
  --fasta /mnt/storage/jonni/UCSC_mouse_genome/mm10.fa \
  --input ./samplesheet_11_09_2025.csv \
  --outdir ./results_11_09_2025 \
  --blacklist ./mm10-blacklist.v2.bed \
  --comparison S2SvsS3,S2LvsS3,S2SvsS4,S2LvsS4,S4vsS3 \
  --compare_groups WT_NPCvsMUT_NPC,WT_NPCvsWT_ESC,WT_ESCvsMUT_ESC \
  --differential_solubility true \
  --gtf /mnt/storage/jonni/UCSC_mouse_genome/mm10.gtf

nextflow run nf-core/sammyseq -r dev \
  -resume \
  -profile docker \
  --fasta /mnt/storage/jonni/UCSC_mouse_genome/mm10.fa \
  --input ./samplesheet_11_09_2025.csv \
  --outdir ./results_11_09_2025 \
  --blacklist ./mm10-blacklist.v2.bed \
  --keep_regions_bed /mnt/storage/jonni/UCSC_mouse_genome/mm10_canonical_chromosomes.bed \
  --comparison S2SvsS3,S2LvsS3,S2SvsS4,S2LvsS4,S4vsS3 \
  --compare_groups WT_NPCvsMUT_NPC,WT_NPCvsWT_ESC,WT_ESCvsMUT_ESC \
  --differential_solubility true \
  --gtf /mnt/storage/jonni/UCSC_mouse_genome/mm10.gtf

nextflow run nf-core/sammyseq -r dev \
  -resume \
  -profile docker \
  --fasta /mnt/storage/jonni/UCSC_mouse_genome/mm10.fa \
  --input ./samplesheet_18_11_2025.csv \
  --outdir ./results_17_11_2025 \
  --blacklist ./mm10-blacklist.v2.bed \
  --keep_regions_bed /mnt/storage/jonni/UCSC_mouse_genome/mm10_canonical_chromosomes.bed \
  --comparison S2SvsS3,S2LvsS3,S2SvsS4,S2LvsS4,S4vsS3 \
  --compare_groups MUT_NPCvsWT_NPC,WT_NPCvsWT_ESC,MUT_ESCvsWT_ESC \
  --differential_solubility true \
  --gtf /mnt/storage/jonni/UCSC_mouse_genome/mm10.gtf



  nextflow run nf-core/sammyseq -r compartments_subworkflow \
  -resume \
  -profile docker \
  --fasta /mnt/storage/jonni/UCSC_mouse_genome/mm10.fa \
  --input ./samplesheet_18_11_2025.csv \
  --outdir ./results_17_11_2025 \
  --blacklist ./mm10-blacklist.v2.bed \
  --keep_regions_bed /mnt/storage/jonni/UCSC_mouse_genome/mm10_canonical_chromosomes2.bed \
  --comparison S2SvsS3,S2LvsS3,S2SvsS4,S2LvsS4,S4vsS3 \
  --compare_groups MUT_NPCvsWT_NPC,WT_NPCvsWT_ESC,MUT_ESCvsWT_ESC \
  --differential_solubility true \
  --compartmentalization_analysis \
  --gtf /mnt/storage/jonni/UCSC_mouse_genome/mm10.gtf


  nextflow run nf-core/sammyseq -r compartments_subworkflow --compartmentalization_analysis