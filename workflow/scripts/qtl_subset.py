import numpy as np
import pandas as pd
from bed_reader import open_bed

# qtl_scope selects the candidate QTL pool for THIS phenotype's subsample:
#   per_phenotype -> only this phenotype's rows
#   all           -> every QTL in the list (deduped by SNP), one shared core-loci panel
scope = snakemake.config["qtl_scope"]
qtl = pd.read_csv(snakemake.config["qtl_list"], sep="\t")
if scope == "per_phenotype":
    qtl = qtl[qtl["PHENOTYPE"] == snakemake.wildcards.phenotype]
elif scope != "all":
    raise ValueError(f"qtl_scope must be 'all' or 'per_phenotype', got {scope!r}")
qtl = qtl.drop_duplicates(subset="SNP").reset_index(drop=True)

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

qtl = qtl[n_classes == 3].reset_index(drop=True)
qtl.to_csv(snakemake.output.qtl, sep="\t", index=False)

# regions derive from the SAME gated set, so background removal, the QTL GRM, and the
# fixed effects can never diverge. plink2 range format: CHR START END [set ID]
qtl[["CHR", "START", "END", "SNP"]].to_csv(
    snakemake.output.regions, sep="\t", header=False, index=False
)
