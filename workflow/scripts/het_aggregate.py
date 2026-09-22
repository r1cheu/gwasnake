import math
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

# classify on BH-FDR, not the raw boundary-LRT p-value
ALPHA = 0.05
RESID_FLOOR_FRAC = 1e-4

rows = []
for full_summary, null_summary, locus in zip(
    snakemake.input.full_summaries,
    snakemake.input.null_summaries,
    snakemake.params.conclusive_loci,
):
    full = pd.read_csv(full_summary, sep="\t")
    null = pd.read_csv(null_summary, sep="\t")
    if (
        "CONVERGENCE_FAILED" in full["term"].values
        or "CONVERGENCE_FAILED" in null["term"].values
    ):
        continue
    sd = full[full["term"] == "zd"].iloc[0]
    sa = full[full["term"] == "za"].iloc[0]
    za_ratio = float(sa["ratio"])
    za_constrained = str(sa["pvalue"]).strip() == "-"
    zd_constrained = str(sd["pvalue"]).strip() == "-"
    za_cols = pd.read_csv(
        Path(full_summary).parent / "za.tsv", sep="\t", nrows=0
    ).columns
    n_strata = len(za_cols) - 2  # drop FID, IID

    resid_full = float(full.loc[full["term"] == "Residual", "estimate"].iloc[0])
    resid_null = float(null.loc[null["term"] == "Residual", "estimate"].iloc[0])
    full_var = full.loc[full["type"] == "variance", "estimate"].astype(float).sum()
    null_var = null.loc[null["type"] == "variance", "estimate"].astype(float).sum()
    degenerate = (
        resid_full / full_var < RESID_FLOOR_FRAC
        or resid_null / null_var < RESID_FLOOR_FRAC
    )

    logl_full = float(full.loc[full["term"] == "logL", "estimate"].iloc[0])
    logl_null = float(null.loc[null["term"] == "logL", "estimate"].iloc[0])
    lrt = max(0.0, 2.0 * (logl_full - logl_null))
    # boundary test: 0.5*chi2_0 + 0.5*chi2_1, so p = 0.5 * P(chi2_1 >= lrt)
    p_lrt = 0.5 * math.erfc(math.sqrt(lrt / 2.0)) if lrt > 0 else 1.0

    rows.append(
        {
            "locus": locus,
            "sigma2_sd": float(sd["estimate"]),
            "se_sd": float(sd["se"]),
            "ratio_sd": float(sd["ratio"]),
            "residual": resid_full,
            "logL_full": logl_full,
            "logL_null": logl_null,
            "LRT": lrt,
            "p_lrt": p_lrt,
            "za_ratio": za_ratio,
            "za_constrained": za_constrained,
            "zd_constrained": zd_constrained,
            "n_strata": n_strata,
            "degenerate": degenerate,
        }
    )

# degenerate fits have no valid p-value -> drop them entirely, do not report
table = pd.DataFrame(rows)
table = table[~table["degenerate"]].drop(columns="degenerate").reset_index(drop=True)

# BH-FDR over the retained loci (statsmodels)
table["fdr"] = np.nan
table["classification"] = pd.Series(dtype=str)
if len(table):
    p = table["p_lrt"].to_numpy()
    _, fdr, _, _ = multipletests(p, alpha=ALPHA, method="fdr_bh")
    table["fdr"] = fdr
    table["classification"] = np.where(
        table["fdr"] < ALPHA, "background_dependent", "stable"
    )
    # split "stable": powered (zd or companion za estimable) vs low_power (both
    # slopes clamped at the floor -> the design could not test stability here)
    both_floored = table["zd_constrained"] & table["za_constrained"]
    table["stable_support"] = np.where(
        table["classification"] == "background_dependent",
        "background_dependent",
        np.where(both_floored, "low_power", "stable"),
    )

table.to_csv(snakemake.output.table, sep="\t", index=False)
