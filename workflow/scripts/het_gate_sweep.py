import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

RUN = Path("/home/rlchen/project/gwasnake/results/6983_peak")
OUT = Path(
    "/home/rlchen/tmp/claude-1000/-home-rlchen-project-gwasnake/"
    "aaeb5a3d-4ccb-40b2-a9e1-1479e5478d8a/scratchpad"
)
# current production params (config.yaml): n_strata=10 het_min=5 hom_min=2 min_strata=6
GRID_N_STRATA = [8, 10, 12, 15, 20]
GRID_HET_MIN = [3, 5, 8]
GRID_HOM_MIN = [2, 3]
GRID_MIN_STRATA = [4, 6, 8, 10]
CUR = dict(n_strata=10, het_min=5, hom_min=2, min_strata=6)

phenos = sorted(p.parent.name for p in RUN.glob("*/heterogeneity"))

by_pheno = []
# per (pheno, combo) -> informative-strata count for each candidate locus that
# passed the gate; pooled later into coverage/power summaries
calib = []  # locus-level: informative strata at CUR params vs observed support

for pheno in phenos:
    base = RUN / pheno
    bed_path = base / "common" / "step2.bed"
    pca_path = base / "common" / "pca.eigenvec"
    qtl_path = base / "qtl" / "qtl.tsv"
    if not (bed_path.exists() and pca_path.exists() and qtl_path.exists()):
        print(f"skip {pheno}: missing inputs")
        continue

    snps = pd.read_csv(qtl_path, sep="\t")["SNP"].tolist()
    pcs = pd.read_csv(pca_path, sep="\t")
    pcs.columns = [c.lstrip("#") for c in pcs.columns]
    pcs["IID"] = pcs["IID"].astype(str)
    pc_cols = [c for c in pcs.columns if c.startswith("PC")]

    bed = open_bed(str(bed_path))
    sid_pos = {s: i for i, s in enumerate(bed.sid)}
    pos = {iid: i for i, iid in enumerate(bed.iid)}
    common = [iid for iid in bed.iid if iid in set(pcs["IID"])]
    row_idx = [pos[iid] for iid in common]
    pcs = pcs.set_index("IID").loc[common]
    X = StandardScaler().fit_transform(pcs[pc_cols].to_numpy(dtype=np.float64))

    present = [s for s in snps if s in sid_pos]
    geno = {
        s: bed.read(index=np.s_[row_idx, [sid_pos[s]]], dtype="float64")[:, 0]
        for s in present
    }

    # observed support labels from the regenerated heterogeneity table
    het_tbl = base / "heterogeneity" / "heterogeneity.tsv"
    support = {}
    if het_tbl.exists():
        t = pd.read_csv(het_tbl, sep="\t")
        support = dict(zip(t["locus"], t.get("stable_support", pd.Series(dtype=str))))

    for n_strata in GRID_N_STRATA:
        label = KMeans(n_clusters=n_strata, random_state=0, n_init=10).fit_predict(X)
        stratum = pd.Series(label)
        # per locus: per-stratum (het, hom_ref, hom_alt) counts
        per_locus = {}
        for s in present:
            g = geno[s]
            c = (
                pd.DataFrame(
                    {
                        "stratum": stratum,
                        "het": g == 1,
                        "hom_ref": g == 0,
                        "hom_alt": g == 2,
                    }
                )
                .groupby("stratum")
                .agg(
                    het=("het", "sum"),
                    hom_ref=("hom_ref", "sum"),
                    hom_alt=("hom_alt", "sum"),
                )
            )
            per_locus[s] = c

        for het_min, hom_min, min_strata in itertools.product(
            GRID_HET_MIN, GRID_HOM_MIN, GRID_MIN_STRATA
        ):
            if min_strata > n_strata:
                continue
            inf_counts = []
            for s in present:
                c = per_locus[s]
                n_inf = int(
                    (
                        (c["het"] >= het_min)
                        & (c["hom_ref"] >= hom_min)
                        & (c["hom_alt"] >= hom_min)
                    ).sum()
                )
                if n_inf >= min_strata:
                    inf_counts.append(n_inf)
                # calibration snapshot at current params
                if (
                    n_strata == CUR["n_strata"]
                    and het_min == CUR["het_min"]
                    and hom_min == CUR["hom_min"]
                    and min_strata == CUR["min_strata"]
                    and n_inf >= min_strata
                ):
                    calib.append(
                        {
                            "pheno": pheno,
                            "locus": s,
                            "n_inf_strata": n_inf,
                            "support": support.get(s, "NA"),
                        }
                    )
            by_pheno.append(
                {
                    "pheno": pheno,
                    "n_strata": n_strata,
                    "het_min": het_min,
                    "hom_min": hom_min,
                    "min_strata": min_strata,
                    "n_candidates": len(present),
                    "n_pass": len(inf_counts),
                    "med_inf_strata": float(np.median(inf_counts)) if inf_counts else 0.0,
                    "min_inf_strata": int(np.min(inf_counts)) if inf_counts else 0,
                }
            )
    print(f"done {pheno}: {len(present)} candidate loci")

bp = pd.DataFrame(by_pheno)
bp.to_csv(OUT / "sweep_by_pheno.tsv", sep="\t", index=False)

# pooled across phenotypes
pooled = (
    bp.groupby(["n_strata", "het_min", "hom_min", "min_strata"])
    .apply(
        lambda d: pd.Series(
            {
                "n_pass": d["n_pass"].sum(),
                "n_candidates": d["n_candidates"].sum(),
                "frac_pass": d["n_pass"].sum() / d["n_candidates"].sum(),
                # coverage-weighted median informative strata across passing loci
                "med_inf_strata": np.average(
                    d["med_inf_strata"], weights=d["n_pass"].clip(lower=1)
                ),
                "min_inf_strata": d["min_inf_strata"][d["n_pass"] > 0].min()
                if (d["n_pass"] > 0).any()
                else 0,
            }
        ),
        include_groups=False,
    )
    .reset_index()
)
pooled.to_csv(OUT / "sweep_pooled.tsv", sep="\t", index=False)

cal = pd.DataFrame(calib)
cal.to_csv(OUT / "sweep_calibration.tsv", sep="\t", index=False)

print("\n=== pooled sweep (sorted by n_pass) ===")
print(pooled.sort_values("n_pass", ascending=False).to_string(index=False))

print("\n=== calibration @ current params: n_inf_strata vs observed support ===")
if len(cal):
    print(
        cal.groupby("support")["n_inf_strata"]
        .agg(["count", "min", "median", "max"])
        .to_string()
    )
