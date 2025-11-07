import pandas as pd

phenotype = pd.read_csv(snakemake.config["phenotype"], sep="\t")[
    snakemake.config["phenotype_name"]
]
halfsib = snakemake.params.halfsib[snakemake.wildcards.group]
filtered = phenotype[phenotype["IID"].apply(lambda x: any(sub in x for sub in halfsib))]

# male = phenotype[
#     (~phenotype["IID"].str.contains("GS1")) & (~phenotype["IID"].str.contains("GS2"))
# ]
# result = pd.concat([filtered, male], axis=0, ignore_index=True)
result = filtered.copy()
result.dropna().to_csv(snakemake.output.phenotype, sep="\t", index=False)
result.iloc[:, :2].to_csv(snakemake.output.sample_list, sep="\t", index=False)
