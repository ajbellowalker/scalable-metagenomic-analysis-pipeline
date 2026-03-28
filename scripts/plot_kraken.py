import pandas as pd
import matplotlib.pyplot as plt
import glob

# Automatically find report
file = glob.glob("results/kraken2/*.report")[0]

df = pd.read_csv(file, sep="\t", header=None)
df.columns = ["percent", "reads_clade", "reads_taxon", "rank", "taxid", "name"]

df["name"] = df["name"].str.strip()

top = df.sort_values("percent", ascending=False).head(10)

plt.figure(figsize=(8,5))
plt.barh(top["name"], top["percent"])
plt.xlabel("Relative abundance (%)")
plt.title("Top 10 Taxa (Kraken2)", fontsize=12)
plt.gca().invert_yaxis()

plt.tight_layout()
plt.savefig("figures/kraken_taxa_barplot.png", dpi=300)