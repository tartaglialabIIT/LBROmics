# Tested system and hardware

## Operating system

Analyses and packaging checks for this repository were run on Linux. A representative tested distribution is:

| Field | Value |
|-------|-------|
| OS | Ubuntu 22.04.5 LTS |
| Codename | Jammy Jellyfish |
| Kernel family | Linux |

```text
PRETTY_NAME="Ubuntu 22.04.5 LTS"
NAME="Ubuntu"
VERSION_ID="22.04"
VERSION="22.04.5 LTS (Jammy Jellyfish)"
ID=ubuntu
```

## Hardware expectations

The analysis workflows were run and tested on Linux-based computing environments. No GPU or other specialized accelerator hardware is required. Raw sequencing preprocessing and nf-core workflows may require a multi-core Linux server or workstation depending on dataset size, whereas the downstream analysis scripts and demonstration can be run on a standard workstation.

**Non-standard hardware:** none required for the downstream analysis code or demonstration. Large-scale preprocessing of sequencing data was performed on Linux servers and does not require GPUs or specialized hardware.
