#!/bin/bash

sudo apt update && sudo apt upgrade -y
sudo apt install -y wget curl git build-essential unzip docker.io openjdk-17-jre

# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Setup environment
conda create -n microbiome python=3.10 -y
conda activate microbiome

conda install -c bioconda -c conda-forge \
megahit bowtie2 samtools metabat2 gtdbtk bedtools hmmer -y

pip install semibin