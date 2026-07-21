# Variance partitioning: the QTL regions are captured by their own GRM (the
# conditional GRM) instead of QTL fixed effects. REML with two GRM sources
# (shared background + QTL region) splits genetic variance into background vs
# QTL-region additive and dominance components.
#
# The QTL GRM is built from the SAME gated regions (qtl_subset.output.regions) that
# qtl.smk removes from the background, so the two GRMs are a clean partition of step2
# (no SNP shared) under either scope. The background is the shared background_grm.


rule qtl_region_bfile:
    input:
        bfile=rules.extract_bed_step2.output.bfile,
        regions=rules.qtl_subset.output.regions,
    output:
        bfile=temp(
            multiext(
                "results/{run_id}/{phenotype}/variance/qtl_region",
                ".bed",
                ".bim",
                ".fam",
            )
        ),
    conda:
        "../envs/plink2.yml"
    log:
        "logs/{run_id}/{phenotype}/qtl_region_bfile.log",
    params:
        step2=rules.extract_bed_step2.params.output_prefix,
        prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/variance/qtl_region",
    shell:
        """
        plink2 --bfile {params.step2} --extract range {input.regions} --make-bed --out {params.prefix} &> {log}
        """


rule qtl_region_grm:
    input:
        bfile=rules.qtl_region_bfile.output.bfile,
    output:
        grm=multiext(
            "results/{run_id}/{phenotype}/variance/qtl_region",
            ".add.bin",
            ".add.id",
            ".dom.bin",
            ".dom.id",
        ),
    threads: config["gelex"]["grm_threads"]
    resources:
        cpus_per_task=threads,
    params:
        prefix=rules.qtl_region_bfile.params.prefix,
    log:
        "logs/{run_id}/{phenotype}/qtl_region_grm.log",
    shell:
        "gelex grm -b {params.prefix} --mode AD -o {params.prefix} -t {threads} --gm NS &> {log}"


rule reml_variance:
    input:
        phenotype=rules.create_sample_list.output.phenotype,
        bg_grm=rules.background_grm.output.grm,
        qtl_grm=rules.qtl_region_grm.output.grm,
    output:
        multiext("results/{run_id}/{phenotype}/variance/reml", summary=".summary", effects=".effects"),
    threads: config["gelex"]["reml_threads"]
    resources:
        cpus_per_task=threads,
    params:
        bg_prefix=rules.background_grm.params.prefix,
        qtl_prefix=rules.qtl_region_grm.params.prefix,
        out_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/variance/reml",
        transform=f"--transform {config['transform']['variance']}" if config["transform"]["variance"] != "none" else "",
    log:
        "logs/{run_id}/{phenotype}/reml_variance.log",
    shell:
        """
        gelex reml -p {input.phenotype} --grm {params.bg_prefix}.add {params.qtl_prefix}.add {params.bg_prefix}.dom {params.qtl_prefix}.dom {params.transform} -o {params.out_prefix} -t {threads} &> {log}
        """


rule variance_table:
    input:
        summary=rules.reml_variance.output.summary,
    output:
        table="results/{run_id}/{phenotype}/variance/variance.tsv",
    conda:
        "../envs/base.yml"
    log:
        "logs/{run_id}/{phenotype}/variance_table.log",
    script:
        "../scripts/variance_table.py"
