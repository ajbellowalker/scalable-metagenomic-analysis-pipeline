#!/bin/bash

# Run Kraken
bash scripts/kraken_processing.py

# Assembly
megahit ...

# Binning
bash scripts/mag_processing.sh

# Functional analysis
bash scripts/functional_analysis.sh
