import gelex
import numpy as np
import pandas as pd
from bed_reader import open_bed
from bed_utils import align_to_bed

qtl = pd.read_csv(snakemake.input.qtl, sep="\t")
samples = pd.read_csv(snakemake.input.samples, sep="\t")
samples.index = samples["IID"].astype(str)

bed = open_bed(snakemake.params.bfile + ".bed")
snp_ids = qtl["SNP"].to_numpy()
sid_pos = {s: i for i, s in enumerate(bed.sid)}
sid_idx = [sid_pos[s] for s in snp_ids]

iids, row_idx = align_to_bed(bed, set(samples.index))

genotype = bed.read(index=np.s_[row_idx, sid_idx], dtype="float64")
additive = np.copy(genotype, order="F")
dominant = np.copy(genotype, order="F")
gelex.encode_inplace(
    additive, effect=gelex.GeneticMode.A, method=gelex.GenotypeMethod.Center
)
gelex.encode_inplace(
    dominant, effect=gelex.GeneticMode.D, method=gelex.GenotypeMethod.Center
)

qcovar = pd.DataFrame({"FID": samples.loc[iids, "FID"].to_numpy(), "IID": iids})
for j, snp in enumerate(snp_ids):
    qcovar[f"{snp}_a"] = additive[:, j]
    qcovar[f"{snp}_d"] = dominant[:, j]

qcovar.to_csv(snakemake.output.qcovar, sep="\t", index=False)
