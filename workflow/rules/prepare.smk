# Create sample list and phenotype file for each group
rule create_sample_list:
    output:
        sample_list="results/{run_id}/{group}/{phenotype}/common/sample.list",
        phenotype="results/{run_id}/{group}/{phenotype}/common/phenotype",
    conda:
        "../envs/base.yml"
    params:
        halfsib=HALFSIB,
    log:
        "logs/{run_id}/{group}/{phenotype}/create_sample_list.log",
    script:
        "../scripts/create_sample_list.py"


# Extract samples for each group from main bfile
rule extract_bed_step1:
    input:
        sample_list=rules.create_sample_list.output.sample_list,
    output:
        bfile=temp(
            multiext(
                "results/{run_id}/{group}/{phenotype}/common/step1",
                ".bed",
                ".bim",
                ".fam",
            )
        ),
    conda:
        "../envs/plink.yml"
    log:
        "logs/{run_id}/{group}/{phenotype}/step1_plink.log",
    params:
        step1=config["bfile"]["step1"],
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/common/step1",
    shell:
        """
        plink --bfile {params.step1} --keep {input.sample_list} --maf 0.01 --geno 0.1 --out {params.output_prefix} --make-bed --threads 1 &> {log}
        """


rule extract_bed_step2:
    input:
        sample_list=rules.create_sample_list.output.sample_list,
    output:
        bfile=temp(
            multiext(
                "results/{run_id}/{group}/{phenotype}/common/step2",
                ".bed",
                ".bim",
                ".fam",
            )
        ),
    conda:
        "../envs/plink.yml"
    log:
        "logs/{run_id}/{group}/{phenotype}/step2_plink.log",
    params:
        step2=config["bfile"]["step2"],
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/common/step2",
    shell:
        """
        plink --bfile {params.step2} --keep {input.sample_list} --maf 0.001 --out {params.output_prefix} --make-bed --threads 1 &> {log}
        """


# Perform PCA analysis for covariates
rule pca:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        pca=temp("results/{run_id}/{group}/{phenotype}/common/pca.eigenvec"),
        eval=temp("results/{run_id}/{group}/{phenotype}/common/pca.eigenval"),
    conda:
        "../envs/plink2.yml"
    threads: 1
    resources:
        cpu_per_task=threads,
    params:
        bfile=rules.extract_bed_step1.params.output_prefix,
        comp=config["pca"],
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/common/pca",
    log:
        "logs/{run_id}/{group}/{phenotype}/pca.log",
    shell:
        """
        plink2 --bfile {params.bfile} --pca {params.comp} --out {params.output_prefix} --threads {threads} &> {log}
        """


# Clean PCA eigenvec file for use as covariates
rule clean_pca_eigenvec:
    input:
        pca=rules.pca.output.pca,
    output:
        covar=temp("results/{run_id}/{group}/{phenotype}/common/qcovar")),
    conda:
        "../envs/base.yml"
    shell:
        "sed '1s/^#//' {input} > {output}"
