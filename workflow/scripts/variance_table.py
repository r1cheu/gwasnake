import pandas as pd

summary = pd.read_csv(snakemake.input.summary, sep="\t")
var = summary[summary["type"] == "variance"].copy()
var["estimate"] = var["estimate"].astype(float)
var["se"] = var["se"].astype(float)
total = var["estimate"].sum()


def component_name(term):
    if term == "Residual":
        return "residual"
    base = term.rsplit("/", 1)[-1]  # background.add | qtl_region.dom
    src, mode = base.rsplit(".", 1)
    return f"{'qtl' if src == 'qtl_region' else src}_{mode}"


var["component"] = var["term"].map(component_name)
# gelex normalizes each GRM by trace/n, so ratio (share of total variance) is the
# comparable quantity across components; qtl PVE = qtl_add ratio + qtl_dom ratio
var["ratio"] = var["estimate"] / total
var["ratio_se"] = pd.to_numeric(var["ratio_se"], errors="coerce")

var[["component", "estimate", "se", "ratio", "ratio_se"]].to_csv(
    snakemake.output.table, sep="\t", index=False
)
