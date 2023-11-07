# deep_sequencing_influenza_variants_identifying

This repository contains scripts and data for setting up and running influenza deep sequencing analysis with variants filtration.

## Setup

### Conda enviroment
The steps to installing the environments are as follows:
- Ensure that you have the Anaconda or Miniconda package manager
- Clone this repo
- Run 
```
conda env create -f requirements.yml
```
- Activate new conda enviroment

## Get reference HA gene sequence

[Download reference manually](https://www.ncbi.nlm.nih.gov/nuccore/KF848938.1?report=fasta) in `fasta` format.

## Proccess experimental and control reads  

You can repeat protocol described in report by using `Snakemake` script provided in this dir. For this purpose execute following command:
```
Snakemake --cores all -p filtered_roommate_variants.csv
```
This script was created through manually implementation of each step:

1. Download and uzip reads in `fastqs` format
2. Align your reads to the reference sequence

