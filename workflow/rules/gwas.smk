wildcard_constraints:
    phenotype="[^_]+",


rule gelex_grm:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        grm=temp(
            multiext(
                "results/{run_id}/{phenotype}/grm",
                ".A.bin",
                ".A.id",
                ".D.bin",
                ".D.id",
            )
        ),
        loco=temp(
            expand(
                "results/{{run_id}}/{{phenotype}}/grm.chr{chr:02d}.{grm_type}.{ext}",
                grm_type=["A", "D"],
                chr=range(1, 13),
                ext=["bin", "id"],
            )
        ),
    log:
        "logs/{run_id}/{phenotype}/gelex_grm.log",
    threads: config["gelex"]["grm_threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile_prefix=rules.extract_bed_step1.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/grm",
    shell:
        """
        gelex grm -b {params.bfile_prefix} --mode AD -o {params.output_prefix} -t {threads} --gm NS &>{log}
        gelex grm -b {params.bfile_prefix} --mode AD -o {params.output_prefix} -t {threads} --loco --gm NS &>>{log}
        """


rule gelex_assoc:
    input:
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=rules.clean_pca_eigenvec.output.covar,
        bfile=rules.extract_bed_step2.output.bfile,
        grm=rules.gelex_grm.output,
    output:
        assoc="results/{run_id}/{phenotype}/assoc.gwas.tsv",
    log:
        "logs/{run_id}/{phenotype}/gelex_assoc.log",
    threads: config["gelex"]["assoc_threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile_prefix=rules.extract_bed_step2.params.output_prefix,
        grm_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/grm",
        transform=config["transform"]["gwas"],
    shell:
        "gelex assoc -b {params.bfile_prefix} -p {input.phenotype} "
        "--mode AD --grm {params.grm_prefix}.A {params.grm_prefix}.D "
        "--transform {params.transform} --qcovar {input.qcovar} --gm NC "
        "-o results/{wildcards.run_id}/{wildcards.phenotype}/assoc "
        "-t {threads} --loco &>{log}"


rule plot_gelex_joint:
    input:
        summary="results/{run_id}/{phenotype}/assoc.gwas.tsv",
    output:
        multiext(
            "results/{run_id}/{phenotype}/manhattan",
            png_a="_A.png",
            png_d="_D.png",
            png_ad="_AD.png",
        ),
    log:
        "logs/{run_id}/{phenotype}/plot_gelex.log",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/plot_gelex_result.py"
