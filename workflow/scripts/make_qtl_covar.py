import sys

import gelex
import numpy as np
import pandas as pd
from bed_reader import open_bed
from bed_utils import align_to_bed

# drop a locus whose a/d columns are ~fully explained (>= 1 - RANK_TOL) by kept columns
RANK_TOL = 1e-6


def full_rank_loci(additive, dominant):
    """Filter loci to a full-rank design via greedy modified Gram-Schmidt.

    A locus is kept only if BOTH its centered additive and dominance columns
    add rank over the intercept and all previously kept columns. This ensures
    the assembled QTL covariate matrix is full column rank.
    """
    n = additive.shape[0]
    basis = [np.full(n, 1.0 / np.sqrt(n))]  # normalized intercept
    keep = np.zeros(additive.shape[1], dtype=bool)
    for j in range(additive.shape[1]):
        trial = list(basis)
        for col in (additive[:, j], dominant[:, j]):
            v = col.astype(np.float64).copy()
            norm0 = np.linalg.norm(v)
            for q in trial:
                v -= (q @ v) * q
            if norm0 == 0 or np.linalg.norm(v) <= RANK_TOL * norm0:
                break
            trial.append(v / np.linalg.norm(v))
        else:
            # executed only if neither column broke out of the loop -> both
            # a and d are linearly independent of the current basis
            basis = trial
            keep[j] = True
    return keep


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
    additive, effect=gelex.GeneticMode.A, method=gelex.GenotypeMethod.NOIACenter
)
gelex.encode_inplace(
    dominant, effect=gelex.GeneticMode.D, method=gelex.GenotypeMethod.NOIACenter
)

keep = full_rank_loci(additive, dominant)
dropped = snp_ids[~keep]
if len(dropped):
    print(
        f"make_qtl_covar: dropped {len(dropped)} collinear loci from conditional "
        f"design: {', '.join(dropped)}",
        file=sys.stderr,
    )
snp_ids = snp_ids[keep]
additive = additive[:, keep]
dominant = dominant[:, keep]

covars = {"FID": samples.loc[iids, "FID"].to_numpy(), "IID": iids}
for j, snp in enumerate(snp_ids):
    covars[f"{snp}_a"] = additive[:, j]
    covars[f"{snp}_d"] = dominant[:, j]

qcovar = pd.DataFrame(covars)
qcovar.to_csv(snakemake.output.qcovar, sep="\t", index=False)
