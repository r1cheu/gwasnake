import pandas as pd
from plot_utils import manhattan_plot


data = pd.read_csv(snakemake.input.summary, sep="\s+")

manhattan_plot(
    data,
    chrom_col="Chr",
    pos_col="bp",
    p_col="p",
    threshold=6,
    log=True,
    title=f"{snakemake.wildcards.group}-{snakemake.wildcards.phenotype} (GCTA)",
    output_file=snakemake.output.png,
    figsize=(18, 6),
)