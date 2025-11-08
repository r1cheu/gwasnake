# REGENIE Step 1: Fit null model and compute LOCO predictions
rule rg_step1:
    input:
        covar=rules.clean_pca_eigenvec.output.covar,
        bfile=rules.extract_bed_step1.output.bfile,
        phenotype=rules.create_sample_list.output.phenotype,
    output:
        output=multiext(
            "results/{run_id}/{group}/{phenotype}/regenie/step1",
            "_1.loco",
            ".log",
            "_pred.list",
        ),
    conda:
        "../envs/regenie.yml"
    threads: config["regenie"]["step1_threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile=rules.extract_bed_step1.params.output_prefix,
        bsize=config["regenie"]["bsize"],
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/regenie/step1",
    shell:
        """
        regenie --step 1 --bed {params.bfile} --covarFile {input.covar} --bsize {params.bsize} --phenoFile {input.phenotype} --out {params.output_prefix} --threads {threads}
        """


# REGENIE Step 2: Perform association testing using LOCO predictions
rule rg_step2:
    input:
        covar=rules.clean_pca_eigenvec.output.covar,
        bfile=rules.extract_bed_step1.output.bfile,
        phenotype=rules.create_sample_list.output.phenotype,
        step1=rules.rg_step1.output.output,
    output:
        output="results/{run_id}/{group}/{phenotype}/regenie/step2_{phenotype}.regenie",
    threads: config["regenie"]["step2_threads"]
    resources:
        cpus_per_task=threads,
    conda:
        "../envs/regenie.yml"
    params:
        bfile=config["bfile"]["step2"],
        bsize=config["regenie"]["bsize"],
        prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/regenie",
    shell:
        """
        regenie --step 2 --bed {params.bfile} --covarFile {input.covar} --bsize {params.bsize} --pred {params.prefix}/step1_pred.list --out {params.prefix}/step2 --phenoFile {input.phenotype} --threads {threads}
        """


# Generate Manhattan plot from REGENIE results
rule plot_rg:
    input:
        summary="results/{run_id}/{group}/{phenotype}/regenie/step2_{phenotype}.regenie",
    output:
        png="results/{run_id}/{group}/{phenotype}/regenie/manhattan.png",
    conda:
        "../envs/base.yml"
    params:
        title=lambda wildcards: f"{wildcards.group}-{wildcards.phenotype}",
    script:
        "../scripts/plot_result.py"
