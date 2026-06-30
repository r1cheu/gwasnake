import pandas as pd

qtl = pd.read_csv(snakemake.config["qtl_list"], sep="\t")
qtl = qtl[qtl["PHENOTYPE"] == snakemake.wildcards.phenotype]

qtl.to_csv(snakemake.output.qtl, sep="\t", index=False)

# plink2 range format: CHR START END [set ID]
qtl[["CHR", "START", "END", "SNP"]].to_csv(
    snakemake.output.regions, sep="\t", header=False, index=False
)
