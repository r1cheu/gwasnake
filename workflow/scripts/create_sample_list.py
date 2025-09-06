import pandas as pd

phenotype = pd.read_csv(snakemake.config["phenotype"], sep="\t")
halfsib = snakemake.params.halfsib[snakemake.wildcards.group]
filtered = phenotype[phenotype["IID"].apply(lambda x: any(sub in x for sub in halfsib))]
filtered.dropna().to_csv(snakemake.output.phenotype, sep="\t", index=False)
filtered.iloc[:, :2].to_csv(snakemake.output.sample_list, sep="\t", index=False)
