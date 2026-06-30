# Task B: is each QTL's dominance effect stable across families?
# Family-specific slopes are passed as --qrand (gelex forms Za Za' / Zd Zd' kernels).
# A GATE (checkpoint) keeps only families informative for the dominance slope and
# skips underpowered loci; the LRT on sigma2_sd (full vs null) tests stability.


wildcard_constraints:
    locus="[^/]+",


checkpoint het_gate:
    input:
        bfile=rules.extract_bed_step2.output.bfile,
    output:
        diagnostic="results/{run_id}/{phenotype}/heterogeneity/het_diagnostic.tsv",
        loci_dir=directory("results/{run_id}/{phenotype}/heterogeneity/loci"),
    params:
        bfile=rules.extract_bed_step2.params.output_prefix,
        loci=lambda wildcards: QTL_LOCI[wildcards.phenotype],
        het_min=config["heterogeneity"]["het_min"],
        hom_min=config["heterogeneity"]["hom_min"],
        min_families=config["heterogeneity"]["min_families"],
    script:
        "../scripts/het_gate.py"


def _locus_file(name):
    def inner(wildcards):
        checkpoints.het_gate.get(run_id=wildcards.run_id, phenotype=wildcards.phenotype)
        return f"results/{wildcards.run_id}/{wildcards.phenotype}/heterogeneity/loci/{wildcards.locus}/{name}"

    return inner


# Null model: additive family slope only (za as nuisance), no dominance slope.
rule reml_het_null:
    input:
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=_locus_file("qcovar"),
        za=_locus_file("za.tsv"),
        grm=rules.background_grm.output.grm,
    output:
        summary="results/{run_id}/{phenotype}/heterogeneity/loci/{locus}/reml.null.summary",
        effects="results/{run_id}/{phenotype}/heterogeneity/loci/{locus}/reml.null.effects",
    threads: config["gelex"]["reml_threads"]
    resources:
        cpus_per_task=threads,
    params:
        grm_prefix=rules.background_grm.params.prefix,
        out_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/heterogeneity/loci/{wildcards.locus}/reml.null",
        transform=f"--transform {config['transform']['heterogeneity']}" if config["transform"]["heterogeneity"] != "none" else "",
    shell:
        """
        gelex reml -p {input.phenotype} --grm {params.grm_prefix}.add {params.grm_prefix}.dom --qcovar {input.qcovar} --qrand {input.za} {params.transform} -o {params.out_prefix} -t {threads}
        """


# Full model: additive + dominance family slopes. LRT vs null tests sigma2_sd = 0.
rule reml_het_full:
    input:
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=_locus_file("qcovar"),
        za=_locus_file("za.tsv"),
        zd=_locus_file("zd.tsv"),
        grm=rules.background_grm.output.grm,
    output:
        summary="results/{run_id}/{phenotype}/heterogeneity/loci/{locus}/reml.full.summary",
        effects="results/{run_id}/{phenotype}/heterogeneity/loci/{locus}/reml.full.effects",
    threads: config["gelex"]["reml_threads"]
    resources:
        cpus_per_task=threads,
    params:
        grm_prefix=rules.background_grm.params.prefix,
        out_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/heterogeneity/loci/{wildcards.locus}/reml.full",
        transform=f"--transform {config['transform']['heterogeneity']}" if config["transform"]["heterogeneity"] != "none" else "",
    shell:
        """
        gelex reml -p {input.phenotype} --grm {params.grm_prefix}.add {params.grm_prefix}.dom --qcovar {input.qcovar} --qrand {input.za} {input.zd} {params.transform} -o {params.out_prefix} -t {threads}
        """


def _conclusive_loci(wildcards):
    ck = checkpoints.het_gate.get(run_id=wildcards.run_id, phenotype=wildcards.phenotype)
    diag = pd.read_csv(ck.output.diagnostic, sep="\t")
    return diag.loc[diag["status"] == "conclusive", "locus"].tolist()


def _conclusive_summaries(model):
    def inner(wildcards):
        return expand(
            "results/{run_id}/{phenotype}/heterogeneity/loci/{locus}/reml.{model}.summary",
            run_id=wildcards.run_id,
            phenotype=wildcards.phenotype,
            locus=_conclusive_loci(wildcards),
            model=model,
        )

    return inner


rule het_aggregate:
    input:
        diagnostic=rules.het_gate.output.diagnostic,
        full_summaries=_conclusive_summaries("full"),
        null_summaries=_conclusive_summaries("null"),
    output:
        table="results/{run_id}/{phenotype}/heterogeneity/heterogeneity.tsv",
    conda:
        "../envs/base.yml"
    params:
        conclusive_loci=_conclusive_loci,
    script:
        "../scripts/het_aggregate.py"
