import gelex
import numpy as np
import pandas as pd
from bed_reader import open_bed

qtl = pd.read_csv(snakemake.input.qtl, sep="\t")
samples = pd.read_csv(snakemake.input.samples, sep="\t")

bed = open_bed(snakemake.params.bfile + ".bed")
snp_ids = qtl["SNP"].to_numpy()
sid_idx = [int(np.where(bed.sid == s)[0][0]) for s in snp_ids]

row_for_iid = {iid: i for i, iid in enumerate(bed.iid)}
samples = samples[samples["IID"].astype(str).isin(row_for_iid)]
iids = samples["IID"].astype(str).to_numpy()
row_idx = [row_for_iid[iid] for iid in iids]

genotype = bed.read(index=np.s_[row_idx, sid_idx], dtype="float64")
additive = np.copy(genotype, order="F")
dominant = np.copy(genotype, order="F")
gelex.encode_inplace(
    additive, effect=gelex.GeneticMode.A, method=gelex.GenotypeMethod.NOIACenter
)
gelex.encode_inplace(
    dominant, effect=gelex.GeneticMode.D, method=gelex.GenotypeMethod.NOIACenter
)

qcovar = pd.DataFrame({"FID": samples["FID"].to_numpy(), "IID": iids})
for j, snp in enumerate(snp_ids):
    qcovar[f"{snp}_a"] = additive[:, j]
    qcovar[f"{snp}_d"] = dominant[:, j]

qcovar.to_csv(snakemake.output.qcovar, sep="\t", index=False)
