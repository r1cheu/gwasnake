# Task B: is each QTL's dominance effect stable across genetic-background strata?
# Strata come from k-means on background PCs; being ~orthogonal to the focal
# genotype, each stratum carries all three genotype classes, so within-stratum
# additive and dominance slopes are separable (half-sib families cannot do this).
# Stratum-specific slopes are passed as --qrand (gelex forms Za Za' / Zd Zd'
# kernels). A GATE (checkpoint) keeps only three-class strata and skips loci with
# too few; the LRT on sigma2_sd (full vs null) tests stability.

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
    params:
        bfile=rules.extract_bed_step2.params.output_prefix,
        n_strata=config["heterogeneity"]["n_strata"],
        het_min=config["heterogeneity"]["het_min"],
        hom_min=config["heterogeneity"]["hom_min"],
        min_strata=config["heterogeneity"]["min_strata"],
    log:
        "logs/{run_id}/{phenotype}/het_gate.log",
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
        multiext("results/{run_id}/{phenotype}/heterogeneity/loci/{locus}/reml.null", summary=".summary", effects=".effects"),
    threads: config["gelex"]["reml_threads"]
    resources:
        cpus_per_task=threads,
    params:
        grm_prefix=rules.background_grm.params.prefix,
        out_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/heterogeneity/loci/{wildcards.locus}/reml.null",
        transform=f"--transform {config['transform']['heterogeneity']}" if config["transform"]["heterogeneity"] != "none" else "",
    log:
        "logs/{run_id}/{phenotype}/heterogeneity/{locus}/reml_null.log",
    shell:
        """
        gelex reml -p {input.phenotype} --grm {params.grm_prefix}.add {params.grm_prefix}.dom --qcovar {input.qcovar} --qrand {input.za} {params.transform} -o {params.out_prefix} -t {threads} &> {log}
        """


# Full model: additive + dominance stratum slopes. LRT vs null tests sigma2_sd = 0.
rule reml_het_full:
    input:
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=het_locus_file("qcovar"),
        za=het_locus_file("za.tsv"),
        zd=het_locus_file("zd.tsv"),
        grm=rules.background_grm.output.grm,
    output:
        multiext("results/{run_id}/{phenotype}/heterogeneity/loci/{locus}/reml.full", summary=".summary", effects=".effects"),
    threads: config["gelex"]["reml_threads"]
    resources:
        cpus_per_task=threads,
    params:
        grm_prefix=rules.background_grm.params.prefix,
        out_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/heterogeneity/loci/{wildcards.locus}/reml.full",
        transform=f"--transform {config['transform']['heterogeneity']}" if config["transform"]["heterogeneity"] != "none" else "",
    log:
        "logs/{run_id}/{phenotype}/heterogeneity/{locus}/reml_full.log",
    shell:
        """
        gelex reml -p {input.phenotype} --grm {params.grm_prefix}.add {params.grm_prefix}.dom --qcovar {input.qcovar} --qrand {input.za} {input.zd} {params.transform} -o {params.out_prefix} -t {threads} &> {log}
        """


rule het_aggregate:
    input:
        full_summaries=het_conclusive_summaries("full"),
        null_summaries=het_conclusive_summaries("null"),
    output:
        table="results/{run_id}/{phenotype}/heterogeneity/heterogeneity.tsv",
    conda:
        "../envs/base.yml"
    log:
        "logs/{run_id}/{phenotype}/het_aggregate.log",
    params:
        conclusive_loci=het_conclusive_loci,
    script:
        "../scripts/het_aggregate.py"
