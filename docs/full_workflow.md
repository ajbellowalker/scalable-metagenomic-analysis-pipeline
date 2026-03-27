### LINUX workflow
## Make Directories
mkdir -p ~/metagenome/{raw_reads,qc_reads,host_removed,kraken,humann,scripts,assemblies,bins,annotation,logs}

ls ~/metagenome
mkdir databases
cd databases

## Setup miniconda 
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

source ~/.bashrc

conda –version

conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

conda create -n microbiome python=3.10 -y
conda activate microbiome

conda install -c bioconda -c conda-forge \
megahit bowtie2 samtools metabat2 gtdbtk -y
pip install semibin
conda install -c bioconda bedtools -y
conda install -c bioconda hmmer -y
sudo apt install seqkit

which seqkit
which megahit
which bowtie2
which samtools
which SemiBin2
which metabat2
which gtdbtk
which bedtools
which hmmsearch

## Install Java
sudo apt update
sudo apt install openjdk-17-jre -y
java -version

## Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version

## Install Docker
sudo apt update
sudo apt install docker.io -y
sudo systemctl start docker
sudo systemctl enable docker
sudo usermod -aG docker ubuntu
sudo chmod 666 /var/run/docker.sock
newgrp docker
docker run hello-world
docker ps

## Pull HUMAnN container
docker pull biobakery/humann
docker run biobakery/humann humann --version
docker run -v ~/metagenome:/data biobakery/humann ls /data

## Download Sequence files
rsync -avz --progress \
-e "ssh -i ~/ .pem key path \
/raw reads path/ \
ubuntu@{IP address}:~/metagenome/raw_reads/
ls ~/metagenome/raw_reads | wc -l

## Download gtdbtk database
mkdir -p databases/gtdbtk
cd ~/metagenome/databases/gtdbtk

wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz

pv gtdbtk_r226_data.tar.gz | tar -xz

du -sh gtdbtk

export GTDBTK_DATA_PATH=~/metagenome/databases/gtdbtk/release226
echo $GTDBTK_DATA_PATH

nano ~/.bashrc 
export GTDBTK_DATA_PATH=/home/ubuntu/metagenome/databases/gtdbtk/release226

source ~/.bashrc

gtdbtk check_install

## Download Kraken database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240112.tar.gz

tar -xvzf k2_standard_20240112.tar.gz

ls -lh ~/metagenome/databases
du -sh ~/metagenome/databases/*

kraken2 --db ~/metagenome/databases/k2_standard help

## Run Kraken2
for r1 in raw_reads/*_R1_001_fastp.fastq.gz
do
base=$(basename $r1 _R1_001_fastp.fastq.gz)
r2=raw_reads/${base}_R2_001_fastp.fastq.gz

kraken2 \
--db /data/home/ubuntu/metagenome/databases/k2_standard \
--paired $r1 $r2 \
--threads 32 \
--use-names \
--report kraken_reports/${base}.report \
--output kraken_reports/${base}.kraken
done

watch -n 15 'du -sh kraken_reports && ls -lh kraken_reports'
ls kraken_reports/*.report | wc -l

awk '$4=="G"' kraken_reports/*.report > genus_abundance.txt


## Create genus matrix
python3 << 'EOF'
import glob
import pandas as pd

files = sorted(glob.glob("kraken_reports/*.report"))
matrix = {}

for f in files:
    sample = f.split("/")[-1].replace(".report","")
    taxa = {}
    
    with open(f) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) > 5 and cols[3] == "G":
                taxa[cols[5].strip()] = float(cols[0])
    
    matrix[sample] = taxa
df = pd.DataFrame(matrix).fillna(0)
df.to_csv("kraken_genus_matrix.csv")
print("Saved:", df.shape)
EOF

head kraken_genus_matrix.csv

## Run Cluster Samples
python3 << 'EOF'
import pandas as pd

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage,fcluster
df = pd.read_csv("kraken_genus_matrix.csv", index_col=0)
dist = pdist(df.T, metric="braycurtis")
link = linkage(dist, method="average")
clusters = fcluster(link,6,criterion="maxclust")
print("\nSample clusters:\n")
for sample,cluster in zip(df.columns,clusters):
    print(sample,"→ cluster",cluster)
EOF

## Run MEGAHIT
# Individual
mkdir -p assemblies/individual

megahit \
-1 raw_reads/*sample*_R1_001_fastp.fastq.gz \
-2 raw_reads/*sample*_R2_001_fastp.fastq.gz \
-o assemblies/individual/sample \
--tmp-dir /data/tmp \
-t 32 \
-m 0.85 \
--k-min 41 \
--k-max 141 \
--k-step 12 \
--min-contig-len 1000
done

watch -n 15 'du -sh assemblies/individual && ls -lh assemblies/individual'

# Coassembly
ls raw_reads/*sample*_R1_001_fastp.fastq.gz \
   raw_reads/*sample *_R1_001_fastp.fastq.gz \
   raw_reads/*sample *_R1_001_fastp.fastq.gz \
   raw_reads/*sample *_R1_001_fastp.fastq.gz \
   raw_reads/*sample *_R1_001_fastp.fastq.gz > coassembly_r1.txt

sed 's/_R1_/_R2_/' coassembly_r1.txt > coassembly_r2.txt
wc -l coassembly_r1.txt coassembly_r2.txt

megahit \
-1 $(paste -sd, coassembly_r1.txt) \
-2 $(paste -sd, coassembly_r2.txt) \
-o assemblies/coassembly_5samples \
--tmp-dir /data/tmp \
-t 32 \
-m 0.85 \
--k-min 41 \
--k-max 141 \
--k-step 12 \
--min-contig-len 1500

watch -n 15 'du -sh assemblies/coassembly_5samples && ls -lh assemblies/coassembly_5samples/intermediate_contigs’

## Download results backup
# Check file size & info
[ -f assemblies/coassembly_5samples/final.contigs.fa ] && echo "ASSEMBLY FINISHED" || echo "STILL RUNNING"

ls -lh assemblies/coassembly_5samples/intermediate_contigs/*final.contigs.fa

grep -c ">" assemblies/coassembly_5samples/intermediate_contigs/k101.contigs.fa

awk '/^>/ {if (seqlen){print seqlen};seqlen=0;next}{seqlen+=length($0)}END{print seqlen}' assemblies/coassembly_5samples/intermediate_contigs/k101.contigs.fa | awk '{sum+=$1} END {print sum/1000000000 " Gb"}'

seqkit stats assemblies/coassembly_5samples/intermediate_contigs/k101.contigs.fa

# Compress & download files to local
gzip assemblies/coassembly_5samples/intermediate_contigs/k101.contigs.fa

tar -czvf results_backup.tar.gz \
assemblies/individual/Ag_009/final.contigs.fa \
assemblies/coassembly_5samples/intermediate_contigs/k101.final.contigs.fa \
assemblies/hybrid_coassembly/k65.final.contigs.fa \
kraken_genus_matrix.csv \
genus_abundance.txt

rm -r assemblies/coassembly_5samples/intermediate_contigs

## Process Assembly 
# Filter the assembly
seqkit seq -m 1500 k101.contigs.fa > coassembly_1.5kb.fa
seqkit seq -m 1500 final.contigs.fa > individual_1.5kb.fa

# Merge Contigs
cat coassembly_1.5kb.fa individual_1.5kb.fa > merged_contigs.fa

# Remove duplicates
seqkit rmdup -s merged_contigs.fa > merged_contigs_dedup.fa

# Rename & check file
seqkit rename merged_contigs_dedup.fa > contigs_final.fa
seqkit stats contigs_final.fa

# Build Bowtie2 index
bowtie2-build contigs_final.fa contigs_index

# Filter the final contigs
seqkit seq -m 2000 contigs_final.fa > contigs_2kb.fa
seqkit stats contigs_2kb.fa

# Map reads to contigs
for r1 in raw_reads/*_R1_001_fastp.fastq.gz
do
sample=$(basename $r1 _R1_001_fastp.fastq.gz)

bowtie2 -x assemblies/contigs_index \
-1 $r1 \
-2 raw_reads/${sample}_R2_001_fastp.fastq.gz \
-p 12 | samtools view -bS - > mapping/${sample}.bam

done

watch -n 15 ls -lh mapping

# Sort BAM files
for bam in *.bam
do
samtools sort -@ 8 -o ${bam%.bam}_sorted.bam $bam
done

# Index sorted BAM files
for bam in *_sorted.bam
do
samtools index $bam
done

# Generate depth file
jjgi_summarize_bam_contig_depths \
--outputDepth depth.txt \
mapping/*_sorted.bam

## Run MetaBAT2
mkdir bins/metabat2

metabat2 \
-i contigs_2kb.fa \
-a depth.txt \
-o bins/metabat2/bin \
-t 16

cd bins/metabat2

for f in *.fa; do
bin=$(basename $f .fa)
grep ">" $f | sed "s/>//g" | awk -v b=$bin '{print $1"\t"b}'
done > metabat_scaffolds2bin.tsv

## Run SemiBin2
mkdir bins/semibin

SemiBin2 single_easy_bin \
-i contigs_2kb.fa \
-b mapping/*_sorted.bam \
-o bins/semibin \
--threads 16

## Run DAS Tool
docker pull quay.io/biocontainers/das_tool:1.1.7--r44hdfd78af_1
docker run --rm quay.io/biocontainers/das_tool:1.1.7--r44hdfd78af_1 DAS_Tool --help

mkdir bins/dastool

docker run --rm \
-v $PWD:$PWD \
-w $PWD \
quay.io/biocontainers/das_tool:1.1.7--r44hdfd78af_1 \
DAS_Tool \
-i bins/metabat2/metabat_final.tsv,bins/semibin/semibin_final.tsv \
-l metabat2,semibin \
-c contigs_2kb.fa \
-o bins/dastool/dastool \
-t 16 \
--score_threshold 0

ls bins/dastool/dastool_DASTool_bins | wc -l

cut -f1 bins/dastool/dastool_DASTool_summary.tsv | tail -n +2 > dastool_selected_bins.txt

wc -l dastool_selected_bins.txt

mkdir bins/dastool_bins

while read bin; do
    grep -w "$bin" bins/dastool/dastool_DASTool_contig2bin.tsv | cut -f1 \
    | seqkit grep -f - contigs_2kb.fa \
    > bins/dastool_bins/${bin}.fa
done < dastool_selected_bins.txt

## MAG Quality filtering
checkm2 predict \
-i bins/dastool_bins \
-x fa \
-o bins/checkm2_dastool \
--threads 16

awk -F'\t' '$2>=90 && $3<=5' bins/checkm2_dastool/quality_report.tsv | wc -l
awk -F'\t' '$2>=50 && $3<=10' bins/checkm2_dastool/quality_report.tsv | wc -l
awk -F"\t" '$2>=50 && $3<=10 {print $1}' bins/checkm2_dastool/quality_report.tsv

awk -F"\t" '
NR>1{
if($2>=90 && $3<=5) HQ++
else if($2>=50 && $3<=10) MQ++
else LQ++
}
END{
print "High quality:",HQ
print "Medium quality:",MQ
print "Low quality:",LQ
}' bins/checkm2_dastool/quality_report.tsv

awk -F"\t" 'NR>1 && $2>=50 && $3<=10 {print $1}' bins/checkm2_dastool/quality_report.tsv > good_MAGs.txt

wc -l good_MAGs.txt

mkdir bins/MAGs_filtered

while read bin
do
cp bins/dastool_bins/${bin}.fa bins/MAGs_filtered/
done < good_MAGs.txt

ls bins/MAGs_filtered | wc -l

tar -czvf MAGs_filtered.tar.gz bins/MAGs_filtered

## Generate MAG QC summary table
awk -F'\t' 'NR==1 || $2>=50 && $3<=10' bins/checkm2_dastool/quality_report.tsv \
> MAG_quality_filtered.tsv

for f in bins/drep_out/dereplicated_genomes/*.fa; do
echo -e "$(basename $f .fa)\t$(grep -v ">" $f | wc -c)\t$(grep -c ">" $f)"
done > MAG_basic_stats.tsv

## Run DeRep
mamba create -n drep_env -c bioconda -c conda-forge \
python=3.10 drep fastani mash pandas=1.5
conda activate drep_env
dRep -h

dRep dereplicate bins/drep_out \
-g bins/MAGs_filtered/*.fa \
-p 16 \
--ignoreGenomeQuality

ls bins/drep_out/dereplicated_genomes | wc -l

## Assign Taxonomy with GTDBTK
gtdbtk classify_wf \
--genome_dir bins/drep_out/dereplicated_genomes \
--extension fa \
--out_dir bins/gtdbtk_out \
--cpus 16

head bins/gtdbtk_out/classify/gtdbtk.bac120.summary.tsv
tar -czvf rumen_MAG_catalogue.tar.gz bins/drep_out/dereplicated_genomes bins/gtdbtk_out

## Generate MAG Taxa summary table
paste \
bins/checkm2_dastool/quality_report.tsv \
bins/gtdbtk_out/classify/gtdbtk.bac120.summary.tsv \
> MAG_summary_raw.tsv

nano mag_summary.sh

#!/bin/bash

echo -e "MAG\tCompleteness\tContamination\tGenomeSize\tTaxonomy" > MAG_summary.tsv

tail -n +2 bins/checkm2_dastool/quality_report.tsv | while read line
do
MAG=$(echo $line | cut -f1)
COMP=$(echo $line | cut -f2)
CONT=$(echo $line | cut -f3)
SIZE=$(echo $line | cut -f4)

TAX=$(grep -w $MAG bins/gtdbtk_out/classify/gtdbtk.bac120.summary.tsv | cut -f2)

echo -e "$MAG\t$COMP\t$CONT\t$SIZE\t$TAX" >> MAG_summary.tsv
done

bash mag_summary.sh

## Detect novel rumen species
awk -F'\t' '$19 < 95' bins/gtdbtk_out/classify/gtdbtk.bac120.summary.tsv > novel_MAGs.tsv
wc -l novel_MAGs.tsv

## Build  MAG phylogenetic tree
ls bins/gtdbtk_out
find bins/gtdbtk_out -name "*.tree"
Go to https://itol.embl.de

## Build MAG + Abundance 
cat bins/drep_out/dereplicated_genomes/*.fa > MAG_catalogue.fa
grep ">" MAG_catalogue.fa | wc -l
bowtie2-build MAG_catalogue.fa MAG_index
mkdir MAG_mapping

for R1 in raw_reads/*_R1_001_fastp.fastq.gz
do
BASE=$(basename $R1 _R1_001_fastp.fastq.gz)
R2=raw_reads/${BASE}_R2_001_fastp.fastq.gz

bowtie2 \
-x MAG_index \
-1 $R1 \
-2 $R2 \
-p 16 \
| samtools view -bS - \
| samtools sort -o MAG_mapping/${BASE}.bam

samtools index MAG_mapping/${BASE}.bam

done

mamba install -c bioconda coverm

coverm genome \
--bam-files MAG_mapping/*.bam \
--genome-fasta-directory bins/drep_out/dereplicated_genomes \
--methods relative_abundance \
-o MAG_abundance.tsv

join -1 1 -2 1 \
<(sort MAG_abundance.tsv) \
<(cut -f1,2 bins/gtdbtk_out/classify/gtdbtk.bac120.summary.tsv | sort) \
> MAG_abundance_taxonomy.tsv

sort -k2 -nr MAG_abundance.tsv | head -30 > top_MAGs.tsv

## Run eggNOG functional annotations
mkdir -p databases/eggnog
cd databases/eggnog

wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz

gunzip eggnog.db.gz
gunzip eggnog_proteins.dmnd.gz
tar -xzf eggnog.taxa.tar.gz

ls -lh

cd ~/metagenome

emapper.py \
-i dram_output/genes.faa \
--data_dir databases/eggnog \
-o MAG_eggnog \
--cpu 12

cut -f1,12 MAG_eggnog.emapper.annotations > MAG_KOs.tsv
cut -f1,7,12,13 MAG_eggnog.emapper.annotations > eggnog_annotations.tsv

# Extract Nitrogen Genes
grep -Ei "glycoside|cellulase|xylanase|CAZy|GH|GT|PL|CE|CBM" eggnog_annotations.tsv > CAZymes.tsv
grep -E "K02586|K02591|K02588|K00370|K00362|K04561|K00376|K01428|K01915" MAG_KOs.tsv > nitrogen_genes.tsv

cut -d "_" -f1 nitrogen_genes.tsv | sort | uniq -c > nitrogen_genes_per_MAG.tsv
cut -f3 nitrogen_genes.tsv | sort | uniq -c | sort -nr > nitrogen_gene_counts.tsv
cut -f1,2 bins/gtdbtk_out/classify/gtdbtk.bac120.summary.tsv > MAG_taxonomy.tsv
cut -d "_" -f1 nitrogen_genes.tsv | sort | uniq -c

sort -k2,2 nitrogen_genes_per_MAG.tsv > sorted_genes.tsv
sort -k1,1 MAG_taxonomy.tsv > sorted_tax.tsv


### R workflow
install.packages("pheatmap")
install.packages("tidyverse")
install.packages(c("ggplot2", "vegan", "ggpubr", "dplyr"))
install.packages("BiocManager")
BiocManager::install("pairwiseAdonis", ask = FALSE, update = FALSE)
install.packages("devtools")
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(devtools)
library(pairwiseAdonis)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(vegan)
library(pheatmap)
library(tidyverse)

data <- read.table("MAG_abundance.tsv", 
                   header=TRUE, 
                   row.names=1,
                   sep="\t")

data_log <- log10(data + 1e-6)

pheatmap(data_log,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         clustering_method="complete",
         scale="row")

pdf("MAG_abundance_heatmap.pdf", width=8, height=10)
pheatmap(data_log,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         clustering_method="complete",
         scale="row")

dev.off()

presence <- data > 0.001
prevalence <- rowSums(presence)
core_MAGs <- prevalence >= ncol(data)/2
sum(core_MAGs)
core_names <- names(core_MAGs[core_MAGs])

grep -Ff core_names bins/gtdbtk_out/classify/gtdbtk.bac120.summary.tsv

hist(prevalence,
     breaks=20,
     col="steelblue",
     main="MAG prevalence across rumen samples",
     xlab="Number of samples")

## Load data for statistical analyses
setwd("working directory path")

abund <- read.csv("MAG_Abundance.csv")
tax <- read.csv("MAG_Taxonomy.csv")

abund_rel <- abund %>%
  filter(MAG != "unmapped") %>%
  mutate(across(-MAG, as.numeric)) %>%
  mutate(across(-MAG, ~ . / sum(.) * 100))

head(abund_rel)

colSums(abund_rel[,-1])

setdiff(abund_rel$MAG, n_wide$MAG)
setdiff(n_wide$MAG, abund_rel$MAG)


n_genes <- read.csv("Nitrogen_genes.csv")

n_long <- n_genes %>%
  separate_rows(N_genes, sep = ",") %>%
  mutate(N_genes = trimws(N_genes)) %>%
  distinct(MAG, Nitrogen_gene_abundance, Gene_function, N_genes)

n_summary <- n_long %>%
  group_by(MAG,Nitrogen_gene_abundance, Gene_function, N_genes) %>%
  summarise(count = n(), .groups = "drop")

n_wide <- n_summary %>%
  pivot_wider(
    names_from = N_genes,
    values_from = count,
    values_fill = 0
  )

clean_MAG <- function(x) {
  x %>%
    gsub("SemIBin\\.", "bin.", .) %>%
    gsub("_", ".", .) %>%
    trimws()
}

abund_rel$MAG <- clean_MAG(abund_rel$MAG)
tax$MAG <- clean_MAG(tax$MAG)
n_wide$MAG <- clean_MAG(n_wide$MAG)

length(intersect(abund_rel$MAG, n_wide$MAG))

n_wide <- n_wide %>%
  filter(MAG %in% abund_rel$MAG)

merged <- abund_rel %>%
  left_join(tax, by = "MAG") %>%
  left_join(n_wide, by = "MAG")

merged[is.na(merged)] <- 0

head(merged)

write.csv(merged, "merged.csv", row.names = FALSE)

## MAG Community Composition and Diversity
# Beta Diversity
abund_mat <- abund_rel %>%
  column_to_rownames("MAG")

head(abund_mat)

# Transpose
abund_t <- t(abund_mat)

# Convert to dataframe
abund_t <- as.data.frame(abund_t)

abund_log <- log10(abund_t + 1e-6)

pca <- prcomp(abund_log, scale. = TRUE)

pca_df <- as.data.frame(pca$x)
pca_df$Sample <- rownames(pca_df)

metadata <- read.csv("Metadata.csv", check.names = FALSE)

# Remove empty columns
metadata <- metadata[, colSums(is.na(metadata)) < nrow(metadata)]
head(metadata$Sample)

pca_df <- pca_df %>%
  left_join(metadata, by = "Sample")

head(pca_df)

write.csv(pca_df, "pca_df.csv", row.names = FALSE)

pca_df$Treatment <- as.character(pca_df$Treatment)

# Fix naming inconsistencies
pca_df$Treatment[pca_df$Treatment == "Control"] <- "C"
pca_df$Treatment[pca_df$Treatment == "Agolin"] <- "O"
pca_df$Treatment[pca_df$Treatment == "AgolinYeaSacc"] <- "OY"

# Convert to factor
pca_df$Treatment <- factor(pca_df$Treatment, levels = c("C", "O", "OY"))
pca_df$Cow <- as.factor(pca_df$Cow)
pca_df$Period <- as.factor(pca_df$Period)

unique(pca_df$Treatment)
unique(pca_df$Cow)

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2)) +

  # Points (treatment only)
  geom_point(aes(color = Treatment), size = 5, alpha = 0.9) +

  # Ellipses
  stat_ellipse(
    aes(color = Treatment, group = Treatment),
    linetype = 2,
    linewidth = 1
  ) +

  # Axis lines at 0
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +

  # Clean theme
 theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black")) +

  # Nice colours
  scale_color_manual(values = c(
    "C" = "#ede109",
    "O" = "#0ec60e",
    "OY" = "#1089cf"
  )) +

  # Labels
  labs(
    x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")
  ) +

  # Make it look like a paper figure
  theme(
    text = element_text(size = 14),
    legend.position = "right"
  )

ggsave("PCA_MAG_abundance.tiff", width = 8, height = 6, dpi = 300)

# Load Abundance Data
abundance_raw <- read.csv("MAG_Abundance.csv", row.names = 1)

# Transpose Data
abundance <- t(abundance_raw)

# Convert to dataframe
abundance <- as.data.frame(abundance)

# PCA metadata
metadata <- pca_df
rownames(metadata) <- metadata$Sample

# Check matching
intersect(rownames(abundance), rownames(metadata))

# Align Datasets
metadata <- metadata[rownames(abundance), ]

# Remove NA
abundance[is.na(abundance)] <- 0

rownames(abundance)
rownames(metadata)

sum(is.na(abundance))

# BRAY-CURTIS + PERMANOVA
library(vegan)

dist_mat <- vegdist(abundance, method = "bray")

table(metadata$Treatment)

adonis_res <- adonis2(dist_mat ~ Treatment, data = metadata)
print(adonis_res)

write.csv(as.data.frame(adonis_res), "Bray-Curtis_PERMANOVA_results.csv", row.names = TRUE)

# Pairwise PERMANOVA
pairwise_permanova <- function(dist_mat, metadata, group_var) {
  
  groups <- unique(metadata[[group_var]])
  results <- data.frame()
  
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      
      g1 <- groups[i]
      g2 <- groups[j]
      
      subset_idx <- metadata[[group_var]] %in% c(g1, g2)
      
      sub_dist <- as.dist(as.matrix(dist_mat)[subset_idx, subset_idx])
      sub_meta <- metadata[subset_idx, ]
      
      res <- adonis2(sub_dist ~ sub_meta[[group_var]])
      
      results <- rbind(results, data.frame(
        Comparison = paste(g1, "vs", g2),
        F = res$F[1],
        R2 = res$R2[1],
        P = res$`Pr(>F)`[1]
      ))
    }
  }
  
  return(results)
}

pairwise_res <- pairwise_permanova(dist_mat, metadata, "Treatment")
pairwise_res$P_adj <- p.adjust(pairwise_res$P, method = "BH")
print(pairwise_res)

write.csv(pairwise_res, "Bray-Curtis_pairwise_PERMANOVA_results.csv", row.names = FALSE)

# Alpha Diversity
abund_t <- t(abund_mat)
abund_t <- as.data.frame(abund_t)

# Shannon
shannon <- diversity(abund_t, index = "shannon")

# Simpson
simpson <- diversity(abund_t, index = "simpson")

# Richness
richness <- specnumber(abund_t)

alpha_df <- data.frame(
  Sample = rownames(abund_t),
  Shannon = shannon,
  Simpson = simpson,
  Richness = richness
)

alpha_df <- alpha_df %>%
  left_join(metadata, by = "Sample")

write.csv(alpha_df, "alpha_diversity.csv", row.names = FALSE)

ggplot(alpha_df, aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic() +
  scale_fill_manual(values = c(
    "Control" = "#ede109",
    "Agolin" = "#0ec60e",
    "AgolinYeaSacc" = "#1089cf"
  )) +
  labs(
    title = "Shannon Diversity of MAG Communities",
    x = "Treatment",
    y = "Shannon Index"
  )

ggsave("Shannon_diversity.png", width = 8, height = 6, dpi = 300)

ggplot(alpha_df, aes(x = Treatment, y = Simpson, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic() +
  scale_fill_manual(values = c(
    "Control" = "#ede109",
    "Agolin" = "#0ec60e",
    "AgolinYeaSacc" = "#1089cf"
  )) +
  labs(
    title = "Simpson Diversity of MAG Communities",
    x = "Treatment",
    y = "Simpson Index"
  )

ggsave("Simpson_diversity.png", width = 8, height = 6, dpi = 300)

ggplot(alpha_df, aes(x = Treatment, y = Richness, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic() +
  scale_fill_manual(values = c(
    "Control" = "#ede109",
    "Agolin" = "#0ec60e",
    "AgolinYeaSacc" = "#1089cf"
  )) +
  labs(
    title = "MAG Richness",
    x = "Treatment",
    y = "Number of MAGs"
  )

ggsave("MAG_richness.png", width = 8, height = 6, dpi = 300)

kruskal.test(Shannon ~ Treatment, data = alpha_df)
kruskal.test(Simpson ~ Treatment, data = alpha_df)
kruskal.test(Richness ~ Treatment, data = alpha_df)

## Plots and Correlations
# Physiological parameters by treatment
ggplot(metadata, aes(x = Treatment, y = Rumen_NH4, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic() +
  scale_fill_manual(values = c(
    "Control" = "#ede109",
    "Agolin" = "#0ec60e",
    "AgolinYeaSacc" = "#1089cf"
  )) +
  labs(
    title = "Rumen Ammonia by Treatment",
    x = "Treatment",
    y = "NH₃ (mg/L)"
  )

ggsave("Rumen_NH3_by_Treatment.png", width = 8, height = 6, dpi = 300)

ggplot(metadata, aes(x = Treatment, y = UREA, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic() +
  scale_fill_manual(values = c(
    "Control" = "#ede109",
    "Agolin" = "#0ec60e",
    "AgolinYeaSacc" = "#1089cf"
  )) +
  labs(
    title = "Urea by Treatment",
    x = "Treatment",
    y = "Urea (mmol/L)"
  )

ggsave("Urea_by_Treatment.png", width = 8, height = 6, dpi = 300)

ggplot(metadata, aes(x = Treatment, y = pH, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic() +
  scale_fill_manual(values = c(
    "Control" = "#ede109",
    "Agolin" = "#0ec60e",
    "AgolinYeaSacc" = "#1089cf"
  )) +
  labs(
    title = "pH by Treatment",
    x = "Treatment",
    y = "pH"
  )

ggsave("pH_by_Treatment.png", width = 8, height = 6, dpi = 300)

# Correlation between PCA and Physiology
PC1_cor_Rumen_NH4 <- cor.test(pca_df$PC1, pca_df$Rumen_NH4)
write.csv(data.frame(
  Correlation = PC1_cor_Rumen_NH4$estimate,
  t = PC1_cor_Rumen_NH4$statistic,
  df = PC1_cor_Rumen_NH4$parameter,
  p_value = PC1_cor_Rumen_NH4$p.value,
  CI_low = PC1_cor_Rumen_NH4$conf.int[1],
  CI_high = PC1_cor_Rumen_NH4$conf.int[2]
), "PC1_Rumen_NH4_correlation.csv", row.names = FALSE)

PC2_cor_Rumen_NH4 <- cor.test(pca_df$PC2, pca_df$Rumen_NH4)
write.csv(data.frame(
  Correlation = PC2_cor_Rumen_NH4$estimate,
  t = PC2_cor_Rumen_NH4$statistic,
  df = PC2_cor_Rumen_NH4$parameter,
  p_value = PC2_cor_Rumen_NH4$p.value,
  CI_low = PC2_cor_Rumen_NH4$conf.int[1],
  CI_high = PC2_cor_Rumen_NH4$conf.int[2]
), "PC2_Rumen_NH4_correlation.csv", row.names = FALSE)

PC1_cor_UREA <- cor.test(pca_df$PC1, pca_df$UREA)
write.csv(data.frame(
  Correlation = PC1_cor_UREA$estimate,
  t = PC1_cor_UREA$statistic,
  df = PC1_cor_UREA$parameter,
  p_value = PC1_cor_UREA$p.value,
  CI_low = PC1_cor_UREA$conf.int[1],
  CI_high = PC1_cor_UREA$conf.int[2]
), "PC1_UREA_correlation.csv", row.names = FALSE)

PC2_cor_UREA <- cor.test(pca_df$PC2, pca_df$UREA)
write.csv(data.frame(
  Correlation = PC2_cor_UREA$estimate,
  t = PC2_cor_UREA$statistic,
  df = PC2_cor_UREA$parameter,
  p_value = PC2_cor_UREA$p.value,
  CI_low = PC2_cor_UREA$conf.int[1],
  CI_high = PC2_cor_UREA$conf.int[2]
), "PC2_UREA_correlation.csv", row.names = FALSE)

PC1_cor_pH <- cor.test(pca_df$PC1, pca_df$pH)
write.csv(data.frame(
  Correlation = PC1_cor_pH$estimate,
  t = PC1_cor_pH$statistic,
  df = PC1_cor_pH$parameter,
  p_value = PC1_cor_pH$p.value,
  CI_low = PC1_cor_pH$conf.int[1],
  CI_high = PC1_cor_pH$conf.int[2]
), "PC1_pH_correlation.csv", row.names = FALSE)

PC2_cor_pH <- cor.test(pca_df$PC2, pca_df$pH)
write.csv(data.frame(
  Correlation = PC2_cor_pH$estimate,
  t = PC2_cor_pH$statistic,
  df = PC2_cor_pH$parameter,
  p_value = PC2_cor_pH$p.value,
  CI_low = PC2_cor_pH$conf.int[1],
  CI_high = PC2_cor_pH$conf.int[2]
), "PC2_pH_correlation.csv", row.names = FALSE)

ggplot(pca_df, aes(x = PC1, y = Rumen_NH4, color = Treatment)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE) +
    scale_color_manual(values = c(
  "C" = "#ede109",
  "O" = "green",
  "OY" = "blue"
)) +
  theme_classic()

ggsave("Rumen_NH4_correlation.png", width = 8, height = 6, dpi = 300)

ggplot(pca_df, aes(x = PC1, y = UREA, color = Treatment)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(
  "C" = "#ede109",
  "O" = "green",
  "OY" = "blue"
)) +
  theme_classic()

ggsave("UREA_correlation.png", width = 8, height = 6, dpi = 300)

ggplot(pca_df, aes(x = PC1, y = pH, color = Treatment)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE) +
    scale_color_manual(values = c(
  "C" = "#ede109",
  "O" = "green",
  "OY" = "blue"
)) +
  theme_classic()

ggsave("Rumen_pH_correlation.png", width = 8, height = 6, dpi = 300)

## Nitrogen Metabolism Genes
# Convert abundance to long format
abund_long <- abund_rel %>%
  pivot_longer(-MAG, names_to = "Sample", values_to = "Abundance")

unique(n_long$MAG)[1:10]
unique(abund_long$MAG)[1:10]

intersect(unique(abund_long$MAG), unique(n_long$MAG))

# Merge with Data
n_sample <- abund_long %>%
  left_join(n_long, by = "MAG")

sum(is.na(n_sample$Gene_function))

n_sample_weighted <- n_sample %>%
  group_by(Sample, Gene_function) %>%
  summarise(weighted_abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

n_sample_wide <- n_sample_weighted %>%
  pivot_wider(
    names_from = N_genes,
    values_from = weighted_abundance,
    values_fill = 0
  )

n_sample_wide <- n_sample_wide %>%
  left_join(metadata, by = "Sample")

colnames(n_sample_wide) <- colnames(n_sample_wide) %>%
  gsub("-", "_", .) %>%
  gsub(" ", "_", .)

kruskal.test(Gln_synt_C ~ Treatment, data = n_sample_wide)

genes <- c("Amidohydro_1","Cytochrom_C552", "Fer2_BFD", "Fer4_NifH", "Gln_synt_C", 
"GSIII_N", "Gln_synt_N", "Molybdopterin","Molydop_binding", 
"Nitr_red_alph_N", "Nitr_red_bet_C","Oxidored_nitro", 
"Pyr_redox_2", "Urease_alpha","Urease_beta")

for (g in genes) {
  print(
    ggplot(n_sample_wide, aes_string(x = "Treatment", y = g, fill = "Treatment")) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.15) +
      theme_classic() +
      ggtitle(g)
  )
    ggsave(paste0(g, "_by_Treatment.png"), width = 8, height = 6, dpi = 300)
}

for (g in genes) {
  ggplot(n_sample_wide, aes_string(x = g, y = "Rumen_NH4", color = "Treatment")) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic()

    ggsave(paste0(g, "_NH4_correlation.png"), width = 8, height = 6, dpi = 300)
}

for (g in genes) {
  ggplot(n_sample_wide, aes_string(x = g, y = "pH", color = "Treatment")) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic()

    ggsave(paste0(g, "_pH_correlation.png"), width = 8, height = 6, dpi = 300)
}

for (g in genes) {
  ggplot(n_sample_wide, aes_string(x = g, y = "UREA", color = "Treatment")) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic()

    ggsave(paste0(g, "_UREA_correlation.png"), width = 8, height = 6, dpi = 300)
}

## Linear Discriminant Analysis
library(dplyr)
library(MASS)

# Select Only Gene Columns
gene_cols <- c("Amidohydro_1","Cytochrom_C552", "Fer2_BFD", "Fer4_NifH", "Gln_synt_C", 
"GSIII_N", "Gln_synt_N", "Molybdopterin","Molydop_binding", 
"Nitr_red_alph_N", "Nitr_red_bet_C","Oxidored_nitro", 
"Pyr_redox_2", "Urease_alpha","Urease_beta"
)

lda_data <- n_sample_wide %>%
  dplyr::select(all_of(gene_cols), Treatment) %>%
  na.omit()

lda_data$Treatment <- as.factor(lda_data$Treatment)

lda_model <- lda(Treatment ~ ., data = lda_data)
write.csv(as.data.frame(lda_model$means), "LDA_means.csv", row.names = TRUE)

lda_results <- as.data.frame(lda_model$scaling)
lda_results$Feature <- rownames(lda_results)

lda_results <- lda_results %>%
  arrange(desc(abs(LD1)))

write.csv(lda_results, "LDA_results.csv", row.names = FALSE)

lda_pred <- predict(lda_model)

lda_df <- as.data.frame(lda_pred$x)
lda_df$Treatment <- lda_data$Treatment

lda_agg <- aggregate(LD1 ~ Treatment, data = lda_df, mean)

write.csv(lda_agg, "LDA_aggregated.csv", row.names = FALSE)

lda_results <- read.csv("LDA_results.csv")

lda_results$Direction <- ifelse(lda_results$LD1 > 0, "C/OY", "O")

ggplot(lda_results, aes(x = reorder(Feature, LD1), y = LD1, fill = Direction)) +
  geom_col() +
  coord_flip() +
  theme_classic()

ggsave("LDA_feature_importance.png", width = 8, height = 6, dpi = 300)

lda_var <- lda_model$svd^2
lda_var_percent <- lda_var / sum(lda_var) * 100

ggplot(lda_df, aes(x = LD1, y = LD2, color = Treatment)) +
  geom_point(size = 4) +

  # Dashed ellipses like your PCA plot
  stat_ellipse(
    level = 0.95,
    linewidth = 1,
    linetype = "dashed"
  ) +

  theme_classic() +
  labs(
  title = "LDA of N Gene Abundance",
  x = paste0("LD1 (", round(lda_var_percent[1], 1), "%)"),
  y = paste0("LD2 (", round(lda_var_percent[2], 1), "%)")
  ) +

  scale_color_manual(values = c(
    "Control" = "#ede109",
    "Agolin" = "#0ec60e",
    "AgolinYeaSacc" = "#1089cf"
  ))
ggsave("LDA_N_Gene_plot.png", width = 8, height = 6, dpi = 300)

pvals_Kruskal <- n_sample_wide %>%
  select(all_of(gene_cols), Treatment) %>%
  pivot_longer(-Treatment, names_to = "Feature", values_to = "Value") %>%
  group_by(Feature) %>%
  summarise(
    p_value = kruskal.test(Value ~ Treatment)$p.value
  ) %>%
  arrange(p_value)

write.csv(pvals_Kruskal, "gene_pvalues_Kruskal.csv", row.names = FALSE)

final_results_Kruskal <- lda_results %>%
  left_join(pvals_Kruskal, by = "Feature") %>%
  arrange(desc(abs(LD1)))

write.csv(final_results_Kruskal, "LDA_with_pvalues_Kruskal.csv", row.names = FALSE)

final_results_Kruskal$p_adj <- p.adjust(final_results_Kruskal$p_value, method = "BH")

write.csv(final_results_Kruskal, "n_genes_stats_results.csv", row.names = FALSE)

## Save Stats Results
# Alpha Diversity
alpha_stats <- data.frame(
  Metric = c("Shannon", "Simpson", "Richness"),
  Chi_sq = c(
    kruskal.test(Shannon ~ Treatment, data = alpha_df)$statistic,
    kruskal.test(Simpson ~ Treatment, data = alpha_df)$statistic,
    kruskal.test(Richness ~ Treatment, data = alpha_df)$statistic
  ),
  df = c(
    kruskal.test(Shannon ~ Treatment, data = alpha_df)$parameter,
    kruskal.test(Simpson ~ Treatment, data = alpha_df)$parameter,
    kruskal.test(Richness ~ Treatment, data = alpha_df)$parameter
  ),
  p_value = c(
    kruskal.test(Shannon ~ Treatment, data = alpha_df)$p.value,
    kruskal.test(Simpson ~ Treatment, data = alpha_df)$p.value,
    kruskal.test(Richness ~ Treatment, data = alpha_df)$p.value
  )
)

write.csv(alpha_stats, "alpha_diversity_stats.csv", row.names = FALSE)

# Beta Diversity
permanova <- adonis2(
  abund_log ~ Treatment,
  data = metadata,
  method = "euclidean"
)

beta_stats <- as.data.frame(permanova)

write.csv(beta_stats, "beta_diversity_PERMANOVA.csv")

# Correlation Results
vars <- c("Amidohydro_1","Cytochrom_C552", "Fer2_BFD", "Fer4_NifH", "Gln_synt_C", 
"GSIII_N", "Gln_synt_N", "Molybdopterin","Molydop_binding", 
"Nitr_red_alph_N", "Nitr_red_bet_C","Oxidored_nitro", 
"Pyr_redox_2", "Urease_alpha","Urease_beta")

NH4_cor_results <- lapply(vars, function(v) {
  test <- cor.test(n_sample_wide[[v]], n_sample_wide$Rumen_NH4)
  
  data.frame(
    Variable = v,
    Correlation = test$estimate,
    p_value = test$p.value,
    CI_low = test$conf.int[1],
    CI_high = test$conf.int[2]
  )
})

NH4_cor_results_df <- do.call(rbind, NH4_cor_results)

write.csv(NH4_cor_results_df, "NH4_multi_correlation_results.csv", row.names = FALSE)

pH_cor_results <- lapply(vars, function(v) {
  test <- cor.test(n_sample_wide[[v]], n_sample_wide$pH)
  
  data.frame(
    Variable = v,
    Correlation = test$estimate,
    p_value = test$p.value,
    CI_low = test$conf.int[1],
    CI_high = test$conf.int[2]
  )
})

pH_cor_results_df <- do.call(rbind, pH_cor_results)

write.csv(pH_cor_results_df, "pH_multi_correlation_results.csv", row.names = FALSE)

urea_cor_results <- lapply(vars, function(v) {
  test <- cor.test(n_sample_wide[[v]], n_sample_wide$UREA)
  
  data.frame(
    Variable = v,
    Correlation = test$estimate,
    p_value = test$p.value,
    CI_low = test$conf.int[1],
    CI_high = test$conf.int[2]
  )
})

urea_cor_results_df <- do.call(rbind, urea_cor_results)

write.csv(urea_cor_results_df, "UREA_multi_correlation_results.csv", row.names = FALSE)

