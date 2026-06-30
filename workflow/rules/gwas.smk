wildcard_constraints:
    phenotype="[^_]+",


rule gelex_grm:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        grm=temp(
            multiext(
                "results/{run_id}/{phenotype}/grm",
                ".add.bin",
                ".add.id",
                ".dom.bin",
                ".dom.id",
            )
        ),
        loco=temp(
            expand(
                "results/{{run_id}}/{{phenotype}}/grm.{grm_type}.chr{chr}.{ext}",
                grm_type=["add", "dom"],
                chr=range(1, 13),
                ext=["bin", "id"],
            )
        ),
    threads: config["gelex"]["grm_threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile_prefix=rules.extract_bed_step1.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/grm",
    shell:
        """
        gelex grm -b {params.bfile_prefix} --add --dom -o {params.output_prefix} -t {threads} --gm NS
        gelex grm -b {params.bfile_prefix} --add --dom -o {params.output_prefix} -t {threads} --loco --gm NS
        """


rule gelex_assoc:
    input:
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=rules.clean_pca_eigenvec.output.covar,
        bfile=rules.extract_bed_step2.output.bfile,
        grm=rules.gelex_grm.output,
    output:
        assoc="results/{run_id}/{phenotype}/joint_assoc.gwas.tsv",
    threads: config["gelex"]["assoc_threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile_prefix=rules.extract_bed_step2.params.output_prefix,
        grm_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/grm",
        transform=f"--transform {config['transform']['gwas']}" if config["transform"]["gwas"] != "none" else "",
    shell:
        """
        gelex assoc -b {params.bfile_prefix} -p {input.phenotype} --grm {params.grm_prefix}.add {params.grm_prefix}.dom --test joint {params.transform} --qcovar {input.qcovar} --gm NC -o results/{wildcards.run_id}/{wildcards.phenotype}/joint_assoc -t {threads} --loco
        """


rule plot_gelex_joint:
    input:
        summary="results/{run_id}/{phenotype}/joint_assoc.gwas.tsv",
    output:
        png_a="results/{run_id}/{phenotype}/joint_manhattan_A.png",
        png_d="results/{run_id}/{phenotype}/joint_manhattan_D.png",
        png_ad="results/{run_id}/{phenotype}/joint_manhattan_AD.png",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/plot_gelex_result.py"
