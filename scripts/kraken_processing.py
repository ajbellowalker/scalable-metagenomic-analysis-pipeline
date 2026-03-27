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
