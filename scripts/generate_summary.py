import pandas as pd

df = pd.read_csv("results/kraken2/sample_1.report", sep="\t", header=None)

summary = {
    "Total taxa": len(df),
    "Top taxon": df.iloc[0, 5]
}

pd.DataFrame([summary]).to_csv("results/pipeline_summary.csv", index=False)
