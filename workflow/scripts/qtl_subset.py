import numpy as np
import pandas as pd
from bed_reader import open_bed

qtl = pd.read_csv(snakemake.config["qtl_list"], sep="\t")
qtl = qtl[qtl["PHENOTYPE"] == snakemake.wildcards.phenotype].reset_index(drop=True)

# regions cover ALL QTL: a locus's region is removed from the background GRM
# regardless of whether the SNP itself is separately estimable in this subsample
# plink2 range format: CHR START END [set ID]
qtl[["CHR", "START", "END", "SNP"]].to_csv(
    snakemake.output.regions, sep="\t", header=False, index=False
)

# step2 bed is already the phenotype subsample; a SNP missing a genotype class
# here makes its additive and dominance covariates collinear (not identifiable)
bed = open_bed(snakemake.params.bfile + ".bed")
sid_pos = {s: i for i, s in enumerate(bed.sid)}
col_idx = qtl["SNP"].map(sid_pos)
genotype = bed.read(
    index=np.s_[:, col_idx.dropna().astype(int).to_numpy()], dtype="float64"
)

n_classes = pd.Series(0, index=qtl.index)
for k, idx in enumerate(qtl.index[col_idx.notna()]):
    col = genotype[:, k]
    n_classes[idx] = len(np.unique(col[~np.isnan(col)]))

keep = n_classes == 3
qtl[keep].to_csv(snakemake.output.qtl, sep="\t", index=False)
