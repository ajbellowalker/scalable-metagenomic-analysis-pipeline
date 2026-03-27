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

