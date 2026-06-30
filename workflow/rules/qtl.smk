# Shared QTL skeleton: subset QTL, remove their regions, LD-prune, build background GRM.
# Consumed by both conditional.smk (task A) and heterogeneity.smk (task B).


rule qtl_subset:
    output:
        qtl="results/{run_id}/{phenotype}/qtl/qtl.tsv",
        regions="results/{run_id}/{phenotype}/qtl/regions.txt",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/qtl_subset.py"


# Remove QTL regions from step2, two-step LD-prune the remainder.
rule background_bfile:
    input:
        bfile=rules.extract_bed_step2.output.bfile,
        regions=rules.qtl_subset.output.regions,
    output:
        bfile=temp(
            multiext(
                "results/{run_id}/{phenotype}/qtl/background",
                ".bed",
                ".bim",
                ".fam",
            )
        ),
    conda:
        "../envs/plink2.yml"
    log:
        "logs/{run_id}/{phenotype}/background_bfile.log",
    params:
        step2=rules.extract_bed_step2.params.output_prefix,
        prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/qtl/background",
        ld0=config["ld_prune"][0],
        ld1=config["ld_prune"][1],
    shell:
        """
        plink2 --bfile {params.step2} --exclude range {input.regions} --indep-pairwise {params.ld0} --out {params.prefix}.p1 &> {log}
        plink2 --bfile {params.step2} --extract {params.prefix}.p1.prune.in --indep-pairwise {params.ld1} --out {params.prefix}.p2 &>> {log}
        plink2 --bfile {params.step2} --extract {params.prefix}.p2.prune.in --make-bed --out {params.prefix} &>> {log}
        """


# Single genome-wide background GRM (QTL regions removed -> no LOCO needed).
rule background_grm:
    input:
        bfile=rules.background_bfile.output.bfile,
    output:
        grm=multiext(
            "results/{run_id}/{phenotype}/qtl/background",
            ".add.bin",
            ".add.id",
            ".dom.bin",
            ".dom.id",
        ),
    threads: config["gelex"]["grm_threads"]
    resources:
        cpus_per_task=threads,
    params:
        prefix=rules.background_bfile.params.prefix,
    shell:
        "gelex grm -b {params.prefix} --add --dom -o {params.prefix} -t {threads} --gm NS"
