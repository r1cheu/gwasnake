import numpy as np
import pandas as pd

summary = pd.read_csv(snakemake.input.summary, sep="\t")
qtl = pd.read_csv(snakemake.input.qtl, sep="\t")
qcovar = pd.read_csv(snakemake.input.qcovar, sep="\t")
phenotype = pd.read_csv(snakemake.input.phenotype, sep="\t")

var_y = phenotype[snakemake.wildcards.phenotype].var(ddof=1)
fixed = summary[summary["type"] == "fixed"].set_index("term")

rows = []
for _, q in qtl.iterrows():
    snp = q["SNP"]
    term_a, term_d = f"{snp}_a", f"{snp}_d"
    if term_a not in fixed.index or term_d not in fixed.index:
        continue
    a = float(fixed.loc[term_a, "estimate"])
    d = float(fixed.loc[term_d, "estimate"])
    xa = qcovar[term_a].to_numpy()
    xd = qcovar[term_d].to_numpy()
    var_a = a**2 * xa.var(ddof=1)
    var_d = d**2 * xd.var(ddof=1)
    cov_ad = 2 * a * d * np.cov(xa, xd, ddof=1)[0, 1]
    rows.append(
        {
            "CHR": q["CHR"],
            "SNP": snp,
            "BP": q["BP"],
            "A1": q["A1"],
            "A2": q["A2"],
            "A": a,
            "A_SE": float(fixed.loc[term_a, "se"]),
            "D": d,
            "D_SE": float(fixed.loc[term_d, "se"]),
            "PVE_A": var_a / var_y,
            "PVE_D": var_d / var_y,
            "PVE_AD": (var_a + var_d + cov_ad) / var_y,
        }
    )

pd.DataFrame(
    rows,
    columns=[
        "CHR",
        "SNP",
        "BP",
        "A1",
        "A2",
        "A",
        "A_SE",
        "D",
        "D_SE",
        "PVE_A",
        "PVE_D",
        "PVE_AD",
    ],
).to_csv(snakemake.output.table, sep="\t", index=False)
