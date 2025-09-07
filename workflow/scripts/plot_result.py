import pandas as pd
from plot_utils import manhattan_plot


data = pd.read_csv(snakemake.input.summary, sep=" ")

manhattan_plot(
    data,
    chrom_col="CHROM",
    pos_col="GENPOS",
    p_col="LOG10P",
    threshold=6,
    title=snakemake.params.title,
    output_file=snakemake.output.png,
    figsize=(18, 6),
)
