from pathlib import Path

import gelex
import numpy as np
import pandas as pd
from bed_reader import open_bed

het_min = snakemake.params.het_min
hom_min = snakemake.params.hom_min
min_families = snakemake.params.min_families
loci = snakemake.params.loci

family = pd.read_csv(
    snakemake.config["heterogeneity"]["family_file"], sep="\t"
).set_index("IID")

bed = open_bed(snakemake.params.bfile + ".bed")
row_for_iid = {iid: i for i, iid in enumerate(bed.iid)}
common = [iid for iid in bed.iid if iid in family.index]
row_idx = [row_for_iid[iid] for iid in common]
fid = family.loc[common, "FID"].to_numpy()
fam_label = pd.Series(family.loc[common, "family"].to_numpy())

loci_dir = Path(snakemake.output.loci_dir)
loci_dir.mkdir(parents=True, exist_ok=True)

diag_rows = []
for snp_id in loci:
    sid_idx = np.where(bed.sid == snp_id)[0]
    if len(sid_idx) == 0:
        diag_rows.append(
            {"locus": snp_id, "n_families": 0, "n_informative": 0,
             "n_het_total": 0, "status": "inconclusive"}
        )
        continue

    genotype = bed.read(index=np.s_[row_idx, sid_idx], dtype="float64")
    g = genotype[:, 0]
    counts = pd.DataFrame(
        {"family": fam_label, "het": g == 1, "hom": (g == 0) | (g == 2)}
    ).groupby("family").agg(het=("het", "sum"), hom=("hom", "sum"))
    informative = counts.index[
        (counts["het"] >= het_min) & (counts["hom"] >= hom_min)
    ].to_numpy()
    n_inf = len(informative)
    status = "conclusive" if n_inf >= min_families else "inconclusive"
    diag_rows.append(
        {"locus": snp_id, "n_families": int(len(counts)), "n_informative": int(n_inf),
         "n_het_total": int((g == 1).sum()), "status": status}
    )
    if status != "conclusive":
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

    # rows = all individuals; columns = informative families only
    incident = (
        pd.get_dummies(fam_label)
        .reindex(columns=informative, fill_value=False)
        .to_numpy(dtype=np.float64)
    )
    za = incident * additive[:, None]
    zd = incident * dominant[:, None]

    out = loci_dir / snp_id
    out.mkdir(parents=True, exist_ok=True)
    for mat, name in ((za, "za.tsv"), (zd, "zd.tsv")):
        df = pd.DataFrame(mat, columns=informative.astype(str))
        df.insert(0, "IID", common)
        df.insert(0, "FID", fid)
        df.to_csv(out / name, sep="\t", index=False)

    qcovar = pd.DataFrame({"FID": fid, "IID": common})
    qcovar[f"{snp_id}_a"] = additive
    qcovar[f"{snp_id}_d"] = dominant
    qcovar.to_csv(out / "qcovar", sep="\t", index=False)

pd.DataFrame(diag_rows).to_csv(snakemake.output.diagnostic, sep="\t", index=False)
