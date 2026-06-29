from pathlib import Path

import gelex
import numpy as np
import pandas as pd
from bed_reader import open_bed


def create_random_slope_matrix(bed_file: str, family_file: str, snp_id):
    family = pd.read_csv(family_file, sep="\t", index_col="IID")
    bed = open_bed(bed_file)

    bed_iid = bed.iid
    common_iid = np.intersect1d(family.index.to_numpy(dtype=str), bed_iid)

    family = family.loc[common_iid]
    row_for_iid = {iid: i for i, iid in enumerate(bed_iid)}
    row_idx = [row_for_iid[iid] for iid in common_iid]
    sid_idx = np.where(bed.sid == snp_id)[0]
    genotype = bed.read(index=np.s_[row_idx, sid_idx], dtype="float64")

    additive = np.copy(genotype, order="F")
    dominant = np.copy(genotype, order="F")

    gelex.encode_inplace(
        additive, effect=gelex.GeneticMode.A, method=gelex.GenotypeMethod.NOIACenter
    )
    gelex.encode_inplace(
        dominant, effect=gelex.GeneticMode.D, method=gelex.GenotypeMethod.NOIACenter
    )

    dummies = pd.get_dummies(family["family"])
    family_names = dummies.columns
    incident_mat = dummies.to_numpy(dtype=np.float64)

    Za = incident_mat * additive
    Zd = incident_mat * dominant

    for mat, path in ((Za, "Za.tsv"), (Zd, "Zd.tsv")):
        df = pd.DataFrame(mat, columns=family_names)
        df.insert(0, "IID", common_iid)
        df.insert(0, "FID", common_iid)
        df.to_csv(path, sep="\t", index=False)
