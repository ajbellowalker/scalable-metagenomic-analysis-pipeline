# LINUX workflow
## Make Directories
```bash
mkdir -p ~/metagenome/{raw_reads,qc_reads,host_removed,kraken,humann,scripts,assemblies,bins,annotation,logs}

ls ~/metagenome
mkdir databases
cd databases
```

## Setup miniconda 
```bash
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
```

## Install Java
```bash
sudo apt update
sudo apt install openjdk-17-jre -y
java -version
```

## Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
```

## Install Docker
```bash
sudo apt update
sudo apt install docker.io -y
sudo systemctl start docker
sudo systemctl enable docker
sudo usermod -aG docker ubuntu
sudo chmod 666 /var/run/docker.sock
newgrp docker
docker run hello-world
docker ps
```

## Pull HUMAnN container
```bash
docker pull biobakery/humann
docker run biobakery/humann humann --version
docker run -v ~/metagenome:/data biobakery/humann ls /data
```

## Download gtdbtk database
```bash
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
```

## Download Kraken database
```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240112.tar.gz

tar -xvzf k2_standard_20240112.tar.gz

ls -lh ~/metagenome/databases
du -sh ~/metagenome/databases/*

kraken2 --db ~/metagenome/databases/k2_standard help
```

## Run Kraken2
```bash
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
```

## Create genus matrix
```py
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
```

## Run Cluster Samples
```py
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
```

## Run MEGAHIT
### Individual
```bash
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
```

### Coassembly
```bash
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
```

## Download results backup
### Check file size & info
```bash
[ -f assemblies/coassembly_5samples/final.contigs.fa ] && echo "ASSEMBLY FINISHED" || echo "STILL RUNNING"

ls -lh assemblies/coassembly_5samples/intermediate_contigs/*final.contigs.fa

grep -c ">" assemblies/coassembly_5samples/intermediate_contigs/k101.contigs.fa

awk '/^>/ {if (seqlen){print seqlen};seqlen=0;next}{seqlen+=length($0)}END{print seqlen}' assemblies/coassembly_5samples/intermediate_contigs/k101.contigs.fa | awk '{sum+=$1} END {print sum/1000000000 " Gb"}'

seqkit stats assemblies/coassembly_5samples/intermediate_contigs/k101.contigs.fa
```

### Compress & download files to local
```bash
gzip assemblies/coassembly_5samples/intermediate_contigs/k101.contigs.fa

tar -czvf results_backup.tar.gz \
assemblies/individual/Ag_009/final.contigs.fa \
assemblies/coassembly_5samples/intermediate_contigs/k101.final.contigs.fa \
assemblies/hybrid_coassembly/k65.final.contigs.fa \
kraken_genus_matrix.csv \
genus_abundance.txt

rm -r assemblies/coassembly_5samples/intermediate_contigs
```

## Process Assembly 
### Filter the assembly
```bash
seqkit seq -m 1500 k101.contigs.fa > coassembly_1.5kb.fa
seqkit seq -m 1500 final.contigs.fa > individual_1.5kb.fa
```

### Merge Contigs
```bash
cat coassembly_1.5kb.fa individual_1.5kb.fa > merged_contigs.fa
```

### Remove duplicates
```bash
seqkit rmdup -s merged_contigs.fa > merged_contigs_dedup.fa
```

### Rename & check file
```bash
seqkit rename merged_contigs_dedup.fa > contigs_final.fa
seqkit stats contigs_final.fa
```

### Build Bowtie2 index
```bash
bowtie2-build contigs_final.fa contigs_index
```

### Filter the final contigs
```bash
seqkit seq -m 2000 contigs_final.fa > contigs_2kb.fa
seqkit stats contigs_2kb.fa
```

### Map reads to contigs
```bash
for r1 in raw_reads/*_R1_001_fastp.fastq.gz
do
sample=$(basename $r1 _R1_001_fastp.fastq.gz)

bowtie2 -x assemblies/contigs_index \
-1 $r1 \
-2 raw_reads/${sample}_R2_001_fastp.fastq.gz \
-p 12 | samtools view -bS - > mapping/${sample}.bam

done

watch -n 15 ls -lh mapping
```

### Sort BAM files
```bash
for bam in *.bam
do
samtools sort -@ 8 -o ${bam%.bam}_sorted.bam $bam
done
```

### Index sorted BAM files
```bash
for bam in *_sorted.bam
do
samtools index $bam
done
```

### Generate depth file
```bash
jjgi_summarize_bam_contig_depths \
--outputDepth depth.txt \
mapping/*_sorted.bam
```

## Run MetaBAT2
```bash
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
```

## Run SemiBin2
```bash
mkdir bins/semibin

SemiBin2 single_easy_bin \
-i contigs_2kb.fa \
-b mapping/*_sorted.bam \
-o bins/semibin \
--threads 16
```

## Run DAS Tool
```bash
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
```

## MAG Quality filtering
```bash
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
```

## Generate MAG QC summary table
```bash
awk -F'\t' 'NR==1 || $2>=50 && $3<=10' bins/checkm2_dastool/quality_report.tsv \
> MAG_quality_filtered.tsv

for f in bins/drep_out/dereplicated_genomes/*.fa; do
echo -e "$(basename $f .fa)\t$(grep -v ">" $f | wc -c)\t$(grep -c ">" $f)"
done > MAG_basic_stats.tsv
```

## Run DeRep
```bash
mamba create -n drep_env -c bioconda -c conda-forge \
python=3.10 drep fastani mash pandas=1.5
conda activate drep_env
dRep -h

dRep dereplicate bins/drep_out \
-g bins/MAGs_filtered/*.fa \
-p 16 \
--ignoreGenomeQuality

ls bins/drep_out/dereplicated_genomes | wc -l
```

## Assign Taxonomy with GTDBTK
```bash
gtdbtk classify_wf \
--genome_dir bins/drep_out/dereplicated_genomes \
--extension fa \
--out_dir bins/gtdbtk_out \
--cpus 16

head bins/gtdbtk_out/classify/gtdbtk.bac120.summary.tsv
tar -czvf rumen_MAG_catalogue.tar.gz bins/drep_out/dereplicated_genomes bins/gtdbtk_out
```

## Generate MAG Taxa summary table
```bash
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
```

## Detect novel rumen species
```bash
awk -F'\t' '$19 < 95' bins/gtdbtk_out/classify/gtdbtk.bac120.summary.tsv > novel_MAGs.tsv
wc -l novel_MAGs.tsv
```

## Build  MAG phylogenetic tree
```bash
ls bins/gtdbtk_out
find bins/gtdbtk_out -name "*.tree"
Go to https://itol.embl.de
```

## Build MAG + Abundance 
```bash
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
```

## Run eggNOG functional annotations
```bash
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
```

## Extract Nitrogen Genes
```bash
grep -Ei "glycoside|cellulase|xylanase|CAZy|GH|GT|PL|CE|CBM" eggnog_annotations.tsv > CAZymes.tsv
grep -E "K02586|K02591|K02588|K00370|K00362|K04561|K00376|K01428|K01915" MAG_KOs.tsv > nitrogen_genes.tsv

cut -d "_" -f1 nitrogen_genes.tsv | sort | uniq -c > nitrogen_genes_per_MAG.tsv
cut -f3 nitrogen_genes.tsv | sort | uniq -c | sort -nr > nitrogen_gene_counts.tsv
cut -f1,2 bins/gtdbtk_out/classify/gtdbtk.bac120.summary.tsv > MAG_taxonomy.tsv
cut -d "_" -f1 nitrogen_genes.tsv | sort | uniq -c

sort -k2,2 nitrogen_genes_per_MAG.tsv > sorted_genes.tsv
sort -k1,1 MAG_taxonomy.tsv > sorted_tax.tsv
```
