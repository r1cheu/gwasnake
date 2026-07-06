from pathlib import Path

import gelex
import numpy as np
import pandas as pd
from bed_reader import open_bed
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

from bed_utils import align_to_bed

n_strata = snakemake.params.n_strata
het_min = snakemake.params.het_min
hom_min = snakemake.params.hom_min
min_strata = snakemake.params.min_strata
loci = pd.read_csv(snakemake.input.qtl, sep="\t")["SNP"].tolist()

# background strata replace half-sib families: a stratum is ~orthogonal to the
# focal genotype, so each carries all three genotype classes -> within-stratum
# additive and dominance slopes are separable (a half-sib family cannot do this)
pcs = pd.read_csv(snakemake.input.pca, sep="\t")
pcs.columns = [c.lstrip("#") for c in pcs.columns]
pcs["IID"] = pcs["IID"].astype(str)
pc_cols = [c for c in pcs.columns if c.startswith("PC")]

bed = open_bed(snakemake.params.bfile + ".bed")
sid_pos = {s: i for i, s in enumerate(bed.sid)}
common, row_idx = align_to_bed(bed, set(pcs["IID"]))
pcs = pcs.set_index("IID").loc[common]
fid = pcs["FID"].to_numpy()

# k-means on standardized background PCs; fixed seed -> reproducible strata
X = StandardScaler().fit_transform(pcs[pc_cols].to_numpy(dtype=np.float64))
label = KMeans(n_clusters=n_strata, random_state=0, n_init=10).fit_predict(X)
stratum = pd.Series(label)

pd.DataFrame({"FID": fid, "IID": common, "stratum": label}).to_csv(
    snakemake.output.strata, sep="\t", index=False
)

loci_dir = Path(snakemake.output.loci_dir)
loci_dir.mkdir(parents=True, exist_ok=True)

for snp_id in loci:
    if snp_id not in sid_pos:
        continue

    genotype = bed.read(index=np.s_[row_idx, [sid_pos[snp_id]]], dtype="float64")
    g = genotype[:, 0]
    counts = (
        pd.DataFrame(
            {"stratum": stratum, "het": g == 1, "hom_ref": g == 0, "hom_alt": g == 2}
        )
        .groupby("stratum")
        .agg(het=("het", "sum"), hom_ref=("hom_ref", "sum"), hom_alt=("hom_alt", "sum"))
    )
    informative = counts.index[
        (counts["het"] >= het_min)
        & (counts["hom_ref"] >= hom_min)
        & (counts["hom_alt"] >= hom_min)
    ].to_numpy()
    if len(informative) < min_strata:
        continue

    additive = np.copy(genotype, order="F")
    dominant = np.copy(genotype, order="F")
    gelex.encode_inplace(
        additive, effect=gelex.GeneticMode.A, method=gelex.GenotypeMethod.NOIACenter
    )
    gelex.encode_inplace(
        dominant, effect=gelex.GeneticMode.D, method=gelex.GenotypeMethod.NOIACenter
    )
    additive = additive[:, 0]
    dominant = dominant[:, 0]

    # rows = all individuals; columns = three-class strata only
    incident = (
        pd.get_dummies(stratum)
        .reindex(columns=informative, fill_value=False)
        .to_numpy(dtype=np.float64)
    )
    za = incident * additive[:, None]
    zd = incident * dominant[:, None]

    out = loci_dir / snp_id
    out.mkdir(parents=True, exist_ok=True)
    for mat, name in ((za, "za.tsv"), (zd, "zd.tsv")):
        df = pd.DataFrame(mat, columns=[f"s{c}" for c in informative])
        df.insert(0, "IID", common)
        df.insert(0, "FID", fid)
        df.to_csv(out / name, sep="\t", index=False)

    qcovar = pd.DataFrame({"FID": fid, "IID": common})
    qcovar[f"{snp_id}_a"] = additive
    qcovar[f"{snp_id}_d"] = dominant
    qcovar.to_csv(out / "qcovar", sep="\t", index=False)
