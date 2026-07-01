import pandas as pd

summary = pd.read_csv(snakemake.input.summary, sep="\t")
qtl = pd.read_csv(snakemake.input.qtl, sep="\t")
effects = pd.read_csv(snakemake.input.effects, sep="\t")

fixed = summary[summary["type"] == "fixed"].set_index("term")

# RINT detaches the raw phenotype variance from the model scale, so the PVE denominator is
# rebuilt from the model: variance of the fitted QTL fixed effects (each effects column is
# already X*beta per individual) plus every REML variance component (background GRMs + residual).
qtl_terms = [c for c in effects.columns if c in fixed.index and c != "Intercept"]
var_fixed = effects[qtl_terms].sum(axis=1).var(ddof=1)
var_components = (
    summary.loc[summary["type"] == "variance", "estimate"].astype(float).sum()
)
var_p = var_fixed + var_components

rows = []
for _, q in qtl.iterrows():
    snp = q["SNP"]
    term_a, term_d = f"{snp}_a", f"{snp}_d"
    if term_a not in fixed.index or term_d not in fixed.index:
        continue
    rows.append(
        {
            "CHR": q["CHR"],
            "SNP": snp,
            "BP": q["BP"],
            "A1": q["A1"],
            "A2": q["A2"],
            "A": float(fixed.loc[term_a, "estimate"]),
            "A_SE": float(fixed.loc[term_a, "se"]),
            "D": float(fixed.loc[term_d, "estimate"]),
            "D_SE": float(fixed.loc[term_d, "se"]),
            "PVE_A": effects[term_a].var(ddof=1) / var_p,
            "PVE_D": effects[term_d].var(ddof=1) / var_p,
            "PVE_AD": (effects[term_a] + effects[term_d]).var(ddof=1) / var_p,
        }
    )

pd.DataFrame(rows).to_csv(snakemake.output.table, sep="\t", index=False)
