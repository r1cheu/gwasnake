import math

import numpy as np
import pandas as pd

# sigma2_sd LRT p-value below this -> dominance effect is background-dependent
ALPHA = 0.05

diag = pd.read_csv(snakemake.input.diagnostic, sep="\t").set_index("locus")

conclusive = {}
for full_summary, null_summary, locus in zip(
    snakemake.input.full_summaries,
    snakemake.input.null_summaries,
    snakemake.params.conclusive_loci,
):
    full = pd.read_csv(full_summary, sep="\t").set_index("term")
    null = pd.read_csv(null_summary, sep="\t").set_index("term")
    sd = full[full.index.str.endswith("/zd.tsv")].iloc[0]

    logl_full = float(full.loc["logL", "estimate"])
    logl_null = float(null.loc["logL", "estimate"])
    lrt = max(0.0, 2.0 * (logl_full - logl_null))
    # boundary test: 0.5*chi2_0 + 0.5*chi2_1, so p = 0.5 * P(chi2_1 >= lrt)
    p_lrt = 0.5 * math.erfc(math.sqrt(lrt / 2.0)) if lrt > 0 else 1.0

    conclusive[locus] = {
        "sigma2_sd": sd["estimate"],
        "se_sd": sd["se"],
        "ratio_sd": sd["ratio"],
        "logL_full": logl_full,
        "logL_null": logl_null,
        "LRT": lrt,
        "p_lrt": p_lrt,
        "classification": "background_dependent" if p_lrt < ALPHA else "stable",
    }

inconclusive = {
    "sigma2_sd": np.nan,
    "se_sd": np.nan,
    "ratio_sd": np.nan,
    "logL_full": np.nan,
    "logL_null": np.nan,
    "LRT": np.nan,
    "p_lrt": np.nan,
    "classification": "inconclusive",
}

rows = [
    {"locus": locus, "n_informative": int(d["n_informative"]),
     **conclusive.get(locus, inconclusive)}
    for locus, d in diag.iterrows()
]
pd.DataFrame(rows).to_csv(snakemake.output.table, sep="\t", index=False)
