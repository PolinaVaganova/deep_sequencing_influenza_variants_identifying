# deep_sequencing_influenza_variants_identifying

This repository contains scripts and data for setting up and running influenza deep sequencing analysis with rare variants filtration.

## Setup

### Conda enviroment
The steps to installing the environments are as follows:
- Ensure that you have the Anaconda or Miniconda package manager
- Clone this repo
- Run 
```
conda env create -f ds_influenza.yaml
```
- Activate new conda enviroment

## Get reference HA gene sequence

[Download reference manually](https://www.ncbi.nlm.nih.gov/nuccore/KF848938.1?report=fasta) in `fasta` format (send to -> file -> FASTA -> Create file) and place it in this repo dir as `deep_sequencing_influenza_variants_identifying/sequence.fasta`.

## Proccess experimental and control reads  

You can repeat protocol described in report by using `Snakemake` script provided in this dir. For this purpose execute following commands:
```
Snakemake --cores all -p reference.control_1.variants.csv reference.control_2.variants.csv reference.control_3.variants.csv filtered.reference.roommate.variants.csv
```
Data about filtered rare SNPs will be collected in `filtered.reference.roommate.variants.csv` file.

