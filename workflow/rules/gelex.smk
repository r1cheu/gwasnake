wildcard_constraints:
    phenotype="[^_]+",
    model="add|dom",


rule gelex_grm_add:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        grm=temp(
            multiext("results/{run_id}/{group}/{phenotype}/grm", ".add.bin", ".add.id")
        ),
        loco=temp(
            expand(
                "results/{{run_id}}/{{group}}/{{phenotype}}/grm.add.chr{chr}.{ext}",
                chr=range(1, 13),
                ext=["bin", "id"],
            )
        ),
    threads: config["gelex"]["grm_threads"]
    resources:
        cpus_per_task=threads,
    conda:
        "../envs/gelex.yml"
    params:
        bfile_prefix=rules.extract_bed_step1.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/grm.add",
        output_prefix_loco=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/grm",
    shell:
        """
        gelex grm -b {params.bfile_prefix} --add -o {params.output_prefix} -t {threads}
        gelex grm -b {params.bfile_prefix} --add -o {params.output_prefix_loco} -t {threads} --loco
        """


rule gelex_grm_dom:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        grm=temp(
            multiext("results/{run_id}/{group}/{phenotype}/grm", ".dom.bin", ".dom.id")
        ),
        loco=temp(
            expand(
                "results/{{run_id}}/{{group}}/{{phenotype}}/grm.dom.chr{chr}.{ext}",
                chr=range(1, 13),
                ext=["bin", "id"],
            )
        ),
    threads: config["gelex"]["grm_threads"]
    resources:
        cpus_per_task=threads,
    conda:
        "../envs/gelex.yml"
    params:
        bfile_prefix=rules.extract_bed_step1.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/grm.dom",
        output_prefix_loco=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/grm",
    shell:
        """
        gelex grm -b {params.bfile_prefix} --dom -o {params.output_prefix} -t {threads}
        gelex grm -b {params.bfile_prefix} --dom -o {params.output_prefix_loco} -t {threads} --loco
        """


rule gelex_assoc_add:
    input:
        grm_add=rules.gelex_grm_add.output.grm,
        grm_loco=rules.gelex_grm_add.output.loco,
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=rules.clean_pca_eigenvec.output.covar,
        bfile=rules.extract_bed_step2.output.bfile,
    output:
        assoc="results/{run_id}/{group}/{phenotype}/add_assoc.gwas.tsv",
    threads: config["gelex"]["assoc_threads"]
    resources:
        cpus_per_task=threads,
    conda:
        "../envs/gelex.yml"
    params:
        bfile_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/common/step2",
        grm_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/grm",
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/add_assoc",
        transform=lambda wildcards: (
            f"--transform {config['transform'][wildcards.phenotype]}"
            if wildcards.phenotype in config.get("transform", {})
            else ""
        ),
    shell:
        """
        gelex assoc -b {params.bfile_prefix} -p {input.phenotype} --grm {params.grm_prefix}.add {params.transform} --qcovar {input.qcovar} -o {params.output_prefix} -t {threads} --loco
        """


rule gelex_assoc_dom:
    input:
        grm_add=rules.gelex_grm_add.output.grm,
        grm_dom=rules.gelex_grm_dom.output.grm,
        grm_add_loco=rules.gelex_grm_add.output.loco,
        grm_dom_loco=rules.gelex_grm_dom.output.loco,
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=rules.clean_pca_eigenvec.output.covar,
        bfile=rules.extract_bed_step2.output.bfile,
    output:
        assoc="results/{run_id}/{group}/{phenotype}/dom_assoc.gwas.tsv",
    threads: config["gelex"]["assoc_threads"]
    resources:
        cpus_per_task=threads,
    conda:
        "../envs/gelex.yml"
    params:
        bfile_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/common/step2",
        grm_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/grm",
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/dom_assoc",
        transform=lambda wildcards: (
            f"--transform {config['transform'][wildcards.phenotype]}"
            if wildcards.phenotype in config.get("transform", {})
            else ""
        ),
    shell:
        """
        gelex assoc -b {params.bfile_prefix} -p {input.phenotype} --grm {params.grm_prefix}.add {params.grm_prefix}.dom --model d {params.transform} --qcovar {input.qcovar} -o {params.output_prefix} -t {threads} --loco
        """


rule plot_gelex:
    input:
        summary="results/{run_id}/{group}/{phenotype}/{model}_assoc.gwas.tsv",
    output:
        png="results/{run_id}/{group}/{phenotype}/{model}_manhattan.png",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/plot_gelex_result.py"
