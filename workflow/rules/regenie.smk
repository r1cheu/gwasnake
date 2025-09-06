rule rg_step1:
    input:
        covar=rules.clean_pca_eigenvec.output.covar,
        bfile=rules.extract_bed.output.bfile,
        phenotype=rules.create_sample_list.output.phenotype,
    output:
        output=multiext(
            "results/{run_id}/{group}/step1", "_1.loco", ".log", "_pred.list"
        ),
    conda:
        "../envs/regenie.yml"
    threads: 8
    resources:
        cpus_per_task=threads,
    params:
        bfile=rules.extract_bed.params.output_prefix,
        bsize=config["regenie"]["bsize"],
        threads=config["regenie"]["step1_threads"],
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/step1",
    shell:
        """
        regenie --step 1 --bed {params.bfile} --covarFile {input.covar} --bsize {params.bsize} --phenoFile {input.phenotype} --out {params.output_prefix} --threads {params.threads}
        """


rule rg_step2:
    input:
        covar=rules.clean_pca_eigenvec.output.covar,
        bfile=rules.extract_bed.output.bfile,
        phenotype=rules.create_sample_list.output.phenotype,
        step1=rules.rg_step1.output.output,
    output:
        output="results/{run_id}/{group}/step2_{phenotype}.regenie",
    threads: 8
    resources:
        cpus_per_task=threads,
    conda:
        "../envs/regenie.yml"
    params:
        bfile=config["bfile"]["step2"],
        bsize=config["regenie"]["bsize"],
        threads=config["regenie"]["step1_threads"],
        prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}",
    shell:
        """
            regenie --step 2 --bed {params.bfile} --covarFile {input.covar} --bsize {params.bsize} --pred {params.prefix}/step1_pred.list --out {params.prefix}/step2 --phenoFile {input.phenotype}  --threads 8
        """


rule plot_rg:
    input:
        summary="results/{run_id}/{group}/step2_{phenotype}.regenie",
    output:
        png="results/{run_id}/{group}/{phenotype}.png",
    conda:
        "../envs/base.yml"
    params:
        title=lambda wildcards: f"{wildcards.group}-{wildcards.phenotype}",
    script:
        "../scripts/plot_result.py"
