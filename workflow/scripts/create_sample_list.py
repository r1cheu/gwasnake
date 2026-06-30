import pandas as pd

phenotype = pd.read_csv(snakemake.config["phenotype"], sep="\t")[
    ["FID", "IID", snakemake.wildcards.phenotype]
].dropna()

phenotype.to_csv(snakemake.output.phenotype, sep="\t", index=False)
phenotype.iloc[:, :2].to_csv(snakemake.output.sample_list, sep="\t", index=False)
