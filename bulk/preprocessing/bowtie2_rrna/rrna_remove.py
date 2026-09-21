#!/usr/bin/env python3
"""Remove rRNA reads with Bowtie2, then re-pair unmapped mates.

Requires Bowtie2 and fastqcombinepairedend.py on PATH.
Override paths via environment variables or edit the defaults below.
"""
import os
import subprocess

FILE_PREFIX = os.environ.get("RRNA_FILE_PREFIX", "trimmed_sample_")
BOWTIE2_INDEX = os.environ.get("RRNA_BOWTIE2_INDEX", "./bowtie_index/mouse_rrnas")
READ1 = os.environ.get("RRNA_READ1", "./trimmed_sample_1.fq.gz")
READ2 = os.environ.get("RRNA_READ2", "./trimmed_sample_2.fq.gz")
OUTDIR = os.environ.get("RRNA_OUTDIR", "./bowtie2_out")
THREADS = os.environ.get("RRNA_THREADS", "8")

os.makedirs(OUTDIR, exist_ok=True)

def run(cmd):
    print(cmd)
    subprocess.check_call(cmd, shell=True)

# Align each mate; keep unmapped reads
run(
    f"bowtie2 -p {THREADS} -x {BOWTIE2_INDEX} -U {READ1} "
    f"-S {OUTDIR}/{FILE_PREFIX}1.sam --un {OUTDIR}/{FILE_PREFIX}1_unmapped"
)
run(
    f"bowtie2 -p {THREADS} -x {BOWTIE2_INDEX} -U {READ2} "
    f"-S {OUTDIR}/{FILE_PREFIX}2.sam --un {OUTDIR}/{FILE_PREFIX}2_unmapped"
)

# Sort unmapped FASTQ by read name
for mate in ("1", "2"):
    run(
        f"cat {OUTDIR}/{FILE_PREFIX}{mate}_unmapped | paste - - - - | "
        f"sort -k1,1 -t ' ' | tr '\\t' '\\n' > {OUTDIR}/{FILE_PREFIX}{mate}_unmapped_sorted"
    )

# Re-pair (adjust barcode prefix '@A00' to match your FASTQ headers)
run(
    f"python fastqcombinepairedend.py '@A00' ' ' "
    f"{OUTDIR}/{FILE_PREFIX}1_unmapped_sorted {OUTDIR}/{FILE_PREFIX}2_unmapped_sorted"
)
