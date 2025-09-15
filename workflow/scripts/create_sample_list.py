from functools import partial

import pandas as pd

phenotype = pd.read_csv(snakemake.config["phenotype"], sep="\t")
halfsib = snakemake.params.halfsib[snakemake.wildcards.group]
filtered = phenotype[phenotype["IID"].apply(lambda x: any(sub in x for sub in halfsib))]
filtered.dropna().to_csv(snakemake.output.phenotype, sep="\t", index=False)
filtered.iloc[:, :2].to_csv(snakemake.output.sample_list, sep="\t", index=False)


def give_sub(x, halfsib):
    for sub in halfsib:
        if sub in x:
            return sub


filtered["mother"] = filtered["IID"].apply(partial(give_sub, halfsib=halfsib))


def get_father(x):
    if "~" in x:
        id_list = x.split("~")
        for id in id_list:
            if "GS3" in id:
                print(id)
                return id
    return "NA"


father = list(set(filtered["IID"].apply(get_father).to_list()))

filtered["father"] = filtered["IID"].apply(partial(give_sub, halfsib=father))

filtered.loc[:, ["FID", "IID", "mother", "father"]].to_csv(
    snakemake.output.covar, sep="\t", index=False, na_rep="NA"
)
