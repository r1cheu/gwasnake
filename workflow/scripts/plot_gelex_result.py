import pandas as pd
from plot_utils import manhattan_plot


data = pd.read_csv(snakemake.input.summary, sep="\t")

manhattan_plot(
    data,
    chrom_col="CHR",
    pos_col="BP",
    p_col="P",
    threshold=6,
    log=True,
    title=f"{snakemake.wildcards.group}-{snakemake.wildcards.phenotype} ({snakemake.wildcards.model})",
    output_file=snakemake.output.png,
    figsize=(18, 6),
)
