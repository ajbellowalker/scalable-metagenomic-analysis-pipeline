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
