# 🧬 Scalable Metagenomic Analysis Pipeline (Nextflow + AWS-ready)

A modular, reproducible metagenomics pipeline built using Nextflow DSL2 for taxonomic profiling, assembly, genome binning, and MAG analysis.

## 🚀 Features

- ⚙️ Fully modular Nextflow DSL2 workflow  
- 🧬 Taxonomic classification using Kraken2  
- 🧱 Assembly using MEGAHIT  
- 🧪 Genome binning with MetaBAT2 (Linux/cloud)  
- 📊 MAG quality assessment with CheckM2  
- 🌍 Taxonomic classification with GTDB-Tk  
- 📈 Automated reporting with MultiQC  
- ☁️ Cloud-ready (AWS / HPC compatible)  

## 🧠 Pipeline Overview

Reads → Kraken2 → MEGAHIT → MetaBAT2 → CheckM2 → GTDB-Tk → MultiQC

## 📦 Installation

```bash
git clone https://github.com/ajbellowalker/scalable-metagenomic-analysis-pipeline
cd scalable-metagenomic-analysis-pipeline
```

## ▶️ Run Pipeline (Local)

```bash
nextflow run main.nf -profile conda
```

## ☁️ Run Full Pipeline (Cloud / HPC)

```bash
nextflow run main.nf -profile conda \
--run_metabat2 true \
--run_checkm2 true \
--run_gtdbtk true
```

## 📊 Output

results/		
├── kraken2/		
├── megahit/		
├── metabat2/			
├── checkm2/		
├── gtdbtk/		
├── multiqc/		

## 📈 Example Output

```markdown
![MultiQC Report](figures/multiqc_report.png)
```

See /figures for more example outputs.

## 🧑‍💻 Author

Ayemere J. Bellowalker		 Bioinformatics | Microbiome | Computational Biology

## 📄 License

MIT