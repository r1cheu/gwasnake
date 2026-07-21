# Common helper functions shared across rule files.

from pathlib import Path


def analysis_targets(analysis):
    """Return the target output files for a given analysis type."""
    if analysis == "gwas":
        return expand(
            "results/{run_id}/{phenotype}/joint_manhattan_{suffix}.png",
            run_id=config["run_id"],
            phenotype=PHENOTYPE,
            suffix=["A", "D", "AD"],
        )
    if analysis == "conditional":
        return expand(
            "results/{run_id}/{phenotype}/conditional/qtl_effects.tsv",
            run_id=config["run_id"],
            phenotype=QTL_PHENOTYPE,
        )
    if analysis == "variance":
        return expand(
            "results/{run_id}/{phenotype}/variance/variance.tsv",
            run_id=config["run_id"],
            phenotype=QTL_PHENOTYPE,
        )
    if analysis == "heterogeneity":
        return expand(
            "results/{run_id}/{phenotype}/heterogeneity/heterogeneity.tsv",
            run_id=config["run_id"],
            phenotype=QTL_PHENOTYPE,
        )
    raise ValueError(f"Unknown analysis: {analysis}")


def qtl_regions(wildcards):
    """Return the QTL regions file path based on config.qtl_scope."""
    scope = config["qtl_scope"]
    if scope == "all":
        return f"results/{wildcards.run_id}/all_qtl_regions.txt"
    if scope == "per_phenotype":
        return f"results/{wildcards.run_id}/{wildcards.phenotype}/qtl/regions.txt"
    raise ValueError(f"qtl_scope must be 'all' or 'per_phenotype', got {scope!r}")


def het_locus_file(name):
    """Return a function that resolves a locus file within the checkpoint output."""
    def inner(wildcards):
        checkpoints.het_gate.get(run_id=wildcards.run_id, phenotype=wildcards.phenotype)
        return f"results/{wildcards.run_id}/{wildcards.phenotype}/heterogeneity/loci/{wildcards.locus}/{name}"
    return inner


def het_conclusive_loci(wildcards):
    """Return sorted list of locus directory names from the het_gate checkpoint."""
    ck = checkpoints.het_gate.get(run_id=wildcards.run_id, phenotype=wildcards.phenotype)
    return sorted(p.name for p in Path(ck.output.loci_dir).iterdir() if p.is_dir())


def het_conclusive_summaries(model):
    """Return a function that expands to all conclusive locus summary paths for a given model."""
    def inner(wildcards):
        return expand(
            "results/{run_id}/{phenotype}/heterogeneity/loci/{locus}/reml.{model}.summary",
            run_id=wildcards.run_id,
            phenotype=wildcards.phenotype,
            locus=het_conclusive_loci(wildcards),
            model=model,
        )
    return inner
