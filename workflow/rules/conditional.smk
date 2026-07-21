# Task A: QTL as fixed effects, conditional on background + each other.

rule make_qtl_covar:
    input:
        qtl=rules.qtl_subset.output.qtl,
        samples=rules.create_sample_list.output.sample_list,
        bfile=rules.extract_bed_step2.output.bfile,
    output:
        qcovar="results/{run_id}/{phenotype}/conditional/qcovar",
    params:
        bfile=rules.extract_bed_step2.params.output_prefix,
    log:
        "logs/{run_id}/{phenotype}/make_qtl_covar.log",
    script:
        "../scripts/make_qtl_covar.py"


rule reml_conditional:
    input:
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=rules.make_qtl_covar.output.qcovar,
        grm=rules.background_grm.output.grm,
    output:
        multiext("results/{run_id}/{phenotype}/conditional/reml", summary=".summary", effects=".effects"),
    threads: config["gelex"]["reml_threads"]
    resources:
        cpus_per_task=threads,
    params:
        grm_prefix=rules.background_grm.params.prefix,
        out_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/conditional/reml",
        transform=f"--transform {config['transform']['conditional']}" if config["transform"]["conditional"] != "none" else "",
    log:
        "logs/{run_id}/{phenotype}/reml_conditional.log",
    shell:
        """
        gelex reml -p {input.phenotype} --grm {params.grm_prefix}.add {params.grm_prefix}.dom --qcovar {input.qcovar} {params.transform} -o {params.out_prefix} -t {threads} &> {log}
        """


rule extract_qtl_effects:
    input:
        summary=rules.reml_conditional.output.summary,
        effects=rules.reml_conditional.output.effects,
        qtl=rules.qtl_subset.output.qtl,
    output:
        table="results/{run_id}/{phenotype}/conditional/qtl_effects.tsv",
    conda:
        "../envs/base.yml"
    log:
        "logs/{run_id}/{phenotype}/extract_qtl_effects.log",
    script:
        "../scripts/extract_qtl_effects.py"
