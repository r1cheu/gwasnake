# Task A: QTL as fixed effects, conditional on background + each other.


rule make_qtl_covar:
    input:
        qtl=rules.qtl_subset.output.qtl,
        samples=rules.create_sample_list.output.sample_list,
        bfile=rules.extract_bed_step2.output.bfile,
    output:
        qcovar="results/{run_id}/{phenotype}/conditional/qcovar",
    log:
        "logs/{run_id}/{phenotype}/make_qtl_covar.log",
    params:
        bfile=rules.extract_bed_step2.params.output_prefix,
    script:
        "../scripts/make_qtl_covar.py"


rule reml_conditional:
    input:
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=rules.make_qtl_covar.output.qcovar,
        grm=rules.background_grm.output.grm,
    output:
        multiext(
            "results/{run_id}/{phenotype}/conditional/reml",
            summary=".summary",
            effects=".effects",
        ),
    log:
        "logs/{run_id}/{phenotype}/reml_conditional.log",
    threads: config["gelex"]["reml_threads"]
    resources:
        cpus_per_task=threads,
    params:
        grm_prefix=rules.background_grm.params.prefix,
        out_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/conditional/reml",
        transform=config["transform"]["conditional"],
    shell:
        "gelex reml -p {input.phenotype} "
        "--grm {params.grm_prefix}.A {params.grm_prefix}.D "
        "--qcovar {input.qcovar} "
        "--transform {params.transform} "
        "-o {params.out_prefix} -t {threads} &>{log}"


rule extract_qtl_effects:
    input:
        summary=rules.reml_conditional.output.summary,
        effects=rules.reml_conditional.output.effects,
        qtl=rules.qtl_subset.output.qtl,
    output:
        table="results/{run_id}/{phenotype}/conditional/qtl_effects.tsv",
    log:
        "logs/{run_id}/{phenotype}/extract_qtl_effects.log",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/extract_qtl_effects.py"
