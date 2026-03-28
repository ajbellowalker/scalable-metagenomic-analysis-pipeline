from Bio import SeqIO
import matplotlib.pyplot as plt

lengths = [len(record.seq) for record in SeqIO.parse("results/megahit/sample_1_assembly/final.contigs.fa", "fasta")]

plt.figure(figsize=(6,4))
plt.hist(lengths, bins=50)
plt.xlabel("Contig length (bp)")
plt.ylabel("Frequency")
plt.title("Assembly Contig Length Distribution", fontsize=12)

plt.tight_layout()
plt.savefig("figures/assembly_stats.png", dpi=300)
