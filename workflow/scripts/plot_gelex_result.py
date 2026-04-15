import pandas as pd
from plot_utils import manhattan_plot


data = pd.read_csv(snakemake.input.summary, sep="\t")
group = snakemake.wildcards.group
phenotype = snakemake.wildcards.phenotype
mode = snakemake.params.mode

if mode == "joint":
    plots = [("P_A", "joint-A", snakemake.output.png_a),
             ("P_D", "joint-D", snakemake.output.png_d)]
else:
    plots = [("P", mode, snakemake.output.png)]

for p_col, label, out in plots:
    manhattan_plot(
        data,
        chrom_col="CHR",
        pos_col="BP",
        p_col=p_col,
        threshold=6,
        log=True,
        title=f"{group}-{phenotype} ({label})",
        output_file=out,
        figsize=(18, 6),
    )
