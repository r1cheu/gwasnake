import pandas as pd
from plot_utils import manhattan_plot


data = pd.read_csv(snakemake.input.summary, sep="\s+", names=['ID','b','se','p'])
data['Chr'] = data['ID'].apply(lambda x: int(x.split("_")[0][-2:])).astype(str)
data['bp'] = data['ID'].str.split("_").str[1].astype(int)
manhattan_plot(
    data,
    chrom_col="Chr",
    pos_col="bp",
    p_col="p",
    threshold=6,
    log=True,
    title=f"{snakemake.wildcards.group}-{snakemake.wildcards.phenotype} (emmax)",
    output_file=snakemake.output.png,
    figsize=(18, 6),
)