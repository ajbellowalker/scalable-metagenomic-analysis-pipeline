# Scalable Metagenomic Analysis Pipeline (Nextflow + AWS)

## Overview
This repository contains a scalable, end-to-end metagenomic workflow for processing large-scale shotgun sequencing datasets. The pipeline was developed for rumen microbiome analysis and supports taxonomic profiling, assembly, genome binning, MAG reconstruction, and functional annotation.

The workflow is designed to run on cloud infrastructure (AWS EC2) using Docker and can be adapted for HPC environments.

## Key Features
- Taxonomic classification using Kraken2
- Assembly using MEGAHIT
- Read mapping using Bowtie2 + Samtools
- Genome binning using MetaBAT2 and SemiBin2
- Bin refinement using DAS Tool
- MAG quality assessment using CheckM2
- Dereplication using dRep
- Taxonomy assignment using GTDB-Tk
- Functional annotation using eggNOG
- Abundance estimation using CoverM

## Dataset Scale
- ~190 GB raw sequencing data
- ~2.3 billion reads processed
- 1,429 MAGs reconstructed
- ~250 medium/high-quality MAGs retained
- ~80 putative novel MAGs identified
- 400+ nitrogen metabolism genes detected

## Tech Stack
- Nextflow
- Docker
- AWS EC2 (Ubuntu)
- Bash
- Python
- R

## Repository Structure
├── main.nf
├── nextflow.config
├── scripts/
├── workflows/
├── docs/
├── modules/
└── assets/

## Pipeline Workflow

### 1. Environment Setup
```bash
conda create -n microbiome python=3.10 -y
conda activate microbiome

conda install -c bioconda -c conda-forge \
megahit bowtie2 samtools metabat2 gtdbtk bedtools hmmer -y

pip install semibin

### 2. Taxonomic Profiling (Kraken2)
```bash
kraken2 \
--db k2_standard \
--paired sample_R1.fastq.gz sample_R2.fastq.gz \
--threads 32 \
--report sample.report

Generates:
•	Taxonomic reports
•	Genus-level abundance matrix

### 3. Assembly (MEGAHIT)
```bash
megahit \
-1 reads_R1.fastq.gz \
-2 reads_R2.fastq.gz \
-o assembly_output \
-t 32

### 4. Read Mapping
```bash
bowtie2 -x contigs_index \
-1 reads_R1.fastq.gz \
-2 reads_R2.fastq.gz | samtools view -bS - > output.bam

### 5. Genome Binning
•	MetaBAT2
•	SemiBin2
•	DAS Tool (integration)

### 6. MAG Quality Filtering
```bash
checkm2 predict \
-i bins \
-o quality_report

Filtering thresholds:
•	High quality: ≥90% completeness, ≤5% contamination
•	Medium quality: ≥50% completeness, ≤10% contamination

### 7. Dereplication
```bash
dRep dereplicate output_dir \
-g bins/*.fa \
-p 16

### 8. Taxonomic Classification
```bash
gtdbtk classify_wf \
--genome_dir bins \
--out_dir gtdbtk_out

### 9. Functional Annotation
```bash
emapper.py \
-i genes.faa \
--data_dir eggnog_db \
-o annotations

### 10. Abundance Estimation
```bash
coverm genome \
--bam-files *.bam \
--genome-fasta-directory genomes \
-o abundance.tsv

## Key Outputs
•	Taxonomic profiles
•	Assembly contigs
•	MAG catalogue
•	MAG quality metrics
•	Functional annotation tables
•	Nitrogen metabolism gene profiles
•	Abundance matrices

## Reproducibility
•	Modular pipeline design
•	Docker-compatible
•	Cloud-ready (AWS)
•	Easily adaptable to HPC clusters

## Notes
•	Raw sequencing data and databases are not included
•	Paths and credentials have been removed for security

## Author
Ayemere J. Bellowalker
PhD Researcher – Microbial Genomics & Bioinformatics


