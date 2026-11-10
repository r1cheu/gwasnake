from pathlib import Path


def het_locus_file(name):
    def inner(wildcards):
        checkpoints.het_gate.get(run_id=wildcards.run_id, phenotype=wildcards.phenotype)
        return f"results/{wildcards.run_id}/{wildcards.phenotype}/heterogeneity/loci/{wildcards.locus}/{name}"

    return inner


def het_conclusive_loci(wildcards):
    ck = checkpoints.het_gate.get(
        run_id=wildcards.run_id, phenotype=wildcards.phenotype
    )
    return sorted(p.name for p in Path(ck.output.loci_dir).iterdir() if p.is_dir())


def het_conclusive_summaries(model):
    def inner(wildcards):
        return expand(
            "results/{run_id}/{phenotype}/heterogeneity/loci/{locus}/reml.{model}.summary",
            run_id=wildcards.run_id,
            phenotype=wildcards.phenotype,
            locus=het_conclusive_loci(wildcards),
            model=model,
        )

    return inner


wildcard_constraints:
    locus="[^/]+",


checkpoint het_gate:
    input:
        bfile=rules.extract_bed_step2.output.bfile,
        qtl=rules.qtl_subset.output.qtl,
        pca=rules.pca.output.pca,
    output:
        loci_dir=directory("results/{run_id}/{phenotype}/heterogeneity/loci"),
        strata="results/{run_id}/{phenotype}/heterogeneity/strata.tsv",
    log:
        "logs/{run_id}/{phenotype}/het_gate.log",
    params:
        bfile=rules.extract_bed_step2.params.output_prefix,
        n_strata=config["heterogeneity"]["n_strata"],
        het_min=config["heterogeneity"]["het_min"],
        hom_min=config["heterogeneity"]["hom_min"],
        min_strata=config["heterogeneity"]["min_strata"],
    script:
        "../scripts/het_gate.py"


# Null model: additive stratum slope only (za as nuisance), no dominance slope.
rule reml_het_null:
    input:
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=het_locus_file("qcovar"),
        za=het_locus_file("za.tsv"),
        grm=rules.background_grm.output.grm,
    output:
        multiext(
            "results/{run_id}/{phenotype}/heterogeneity/loci/{locus}/reml.null",
            summary=".summary",
            effects=".effects",
        ),
    log:
        "logs/{run_id}/{phenotype}/heterogeneity/{locus}/reml_null.log",
    threads: config["gelex"]["reml_threads"]
    resources:
        cpus_per_task=threads,
    params:
        grm_prefix=rules.background_grm.params.prefix,
        out_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/heterogeneity/loci/{wildcards.locus}/reml.null",
        transform=config["transform"],
    shell:
        "gelex reml -p {input.phenotype} "
        "--random {params.grm_prefix}.A:g {params.grm_prefix}.D:g {input.za}:q "
        "--fixed {input.qcovar}:q "
        "--transform {params.transform} "
        "-o {params.out_prefix} -t {threads} &>{log} || {{ "
        "printf 'term\\ttype\\testimate\\tse\\tratio\\tratio_se\\tpvalue\\n"
        "CONVERGENCE_FAILED\\tstatus\\tnan\\tnan\\tnan\\tnan\\tnan\\n' "
        ">{params.out_prefix}.summary; : >{params.out_prefix}.effects; }}"


# Full model: additive + dominance stratum slopes. LRT vs null tests sigma2_sd = 0.
rule reml_het_full:
    input:
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=het_locus_file("qcovar"),
        za=het_locus_file("za.tsv"),
        zd=het_locus_file("zd.tsv"),
        grm=rules.background_grm.output.grm,
    output:
        multiext(
            "results/{run_id}/{phenotype}/heterogeneity/loci/{locus}/reml.full",
            summary=".summary",
            effects=".effects",
        ),
    log:
        "logs/{run_id}/{phenotype}/heterogeneity/{locus}/reml_full.log",
    threads: config["gelex"]["reml_threads"]
    resources:
        cpus_per_task=threads,
    params:
        grm_prefix=rules.background_grm.params.prefix,
        out_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/heterogeneity/loci/{wildcards.locus}/reml.full",
        transform=config["transform"],
    shell:
        "gelex reml -p {input.phenotype} "
        "--random {params.grm_prefix}.A:g {params.grm_prefix}.D:g {input.za}:q {input.zd}:q "
        "--fixed {input.qcovar}:q "
        "--transform {params.transform} "
        "-o {params.out_prefix} -t {threads} &>{log} || {{ "
        "printf 'term\\ttype\\testimate\\tse\\tratio\\tratio_se\\tpvalue\\n"
        "CONVERGENCE_FAILED\\tstatus\\tnan\\tnan\\tnan\\tnan\\tnan\\n' "
        ">{params.out_prefix}.summary; : >{params.out_prefix}.effects; }}"


rule het_aggregate:
    input:
        full_summaries=het_conclusive_summaries("full"),
        null_summaries=het_conclusive_summaries("null"),
    output:
        table="results/{run_id}/{phenotype}/heterogeneity/heterogeneity.tsv",
    log:
        "logs/{run_id}/{phenotype}/het_aggregate.log",
    conda:
        "../envs/base.yml"
    params:
        conclusive_loci=het_conclusive_loci,
    script:
        "../scripts/het_aggregate.py"
