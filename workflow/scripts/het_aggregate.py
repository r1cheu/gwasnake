import math
from pathlib import Path

import numpy as np
import pandas as pd

# classify on BH-FDR, not the raw boundary-LRT p-value
ALPHA = 0.05
# residual below this fraction of total variance = AI-REML collapsed it to the
# ~1e-6*Vp floor -> logL surface degenerate, the LRT is void (not a real signal)
RESID_FLOOR_FRAC = 1e-4
# a slope-variance ratio below this = AI-REML pinned it at the ~1e-6*Vp floor
SLOPE_FLOOR_FRAC = 1e-4

rows = []
for full_summary, null_summary, locus in zip(
    snakemake.input.full_summaries,
    snakemake.input.null_summaries,
    snakemake.params.conclusive_loci,
):
    full = pd.read_csv(full_summary, sep="\t")
    null = pd.read_csv(null_summary, sep="\t")
    sd = full[full["term"].str.endswith("/zd.tsv")].iloc[0]
    # companion additive slope on the SAME strata: if it too is pinned at the
    # floor the strata carry no estimable slope signal -> a floored zd is "no
    # power", not "genuinely stable". n_strata is the fit-independent info count.
    sa = full[full["term"].str.endswith("/za.tsv")].iloc[0]
    za_ratio = float(sa["ratio"])
    za_cols = pd.read_csv(Path(full_summary).parent / "za.tsv", sep="\t", nrows=0).columns
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
            "n_strata": n_strata,
            "degenerate": degenerate,
        }
    )

# degenerate fits have no valid p-value -> drop them entirely, do not report
table = pd.DataFrame(rows)
table = table[~table["degenerate"]].drop(columns="degenerate").reset_index(drop=True)

# BH-FDR over the retained loci
table["fdr"] = np.nan
table["classification"] = pd.Series(dtype=str)
if len(table):
    p = table["p_lrt"].to_numpy()
    n = len(p)
    order = p.argsort()
    q = np.clip(np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1], 0, 1)
    fdr = np.empty(n)
    fdr[order] = q
    table["fdr"] = fdr
    table["classification"] = np.where(
        table["fdr"] < ALPHA, "background_dependent", "stable"
    )
    # split "stable": powered (zd or companion za estimable) vs low_power (both
    # slopes pinned at the floor -> the design could not test stability here)
    both_floored = (table["ratio_sd"] < SLOPE_FLOOR_FRAC) & (
        table["za_ratio"] < SLOPE_FLOOR_FRAC
    )
    table["stable_support"] = np.where(
        table["classification"] == "background_dependent",
        "background_dependent",
        np.where(both_floored, "low_power", "stable"),
    )

table.to_csv(snakemake.output.table, sep="\t", index=False)
