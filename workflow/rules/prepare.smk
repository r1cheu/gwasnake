# Create sample list and phenotype file for each group
rule create_sample_list:
    output:
        sample_list="results/{run_id}/{group}/common/sample.list",
        phenotype="results/{run_id}/{group}/common/phenotype",
        covar="results/{run_id}/{group}/common/covar",
    conda:
        "../envs/base.yml"
    params:
        halfsib=HALFSIB,
    log:
        "logs/{run_id}/{group}/create_sample_list.log",
    script:
        "../scripts/create_sample_list.py"


# Extract samples for each group from main bfile
rule extract_bed_step1:
    input:
        sample_list=rules.create_sample_list.output.sample_list,
    output:
        bfile=multiext("results/{run_id}/{group}/common/step1", ".bed", ".bim", ".fam"),
    conda:
        "../envs/plink.yml"
    params:
        step1=config["bfile"]["step1"],
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/common/step1",
    shell:
        """
        plink --bfile {params.step1} --keep {input.sample_list} --maf 0.01 --geno 0.1 --out {params.output_prefix} --make-bed --threads 1
        """

rule extract_bed_step2:
    input:
        sample_list=rules.create_sample_list.output.sample_list,
    output:
        bfile=multiext("results/{run_id}/{group}/common/step2", ".bed", ".bim", ".fam"),
    conda:
        "../envs/plink.yml"
    params:
        step2=config["bfile"]["step2"],
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/common/step2",
    shell:
        """
        plink --bfile {params.step2} --keep {input.sample_list} --maf 0.01 --geno 0.1 --out {params.output_prefix} --make-bed --threads 1
        """


# Perform PCA analysis for covariates
rule pca:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        pca=temp("results/{run_id}/{group}/common/pca.eigenvec"),
        eval=temp("results/{run_id}/{group}/common/pca.eigenval"),
    conda:
        "../envs/plink2.yml"
    threads: 2
    resources:
        cpu_per_task=threads,
    params:
        bfile=rules.extract_bed_step1.params.output_prefix,
        comp=config["pca"],

    log:
        "logs/{run_id}/{group}/pca.log",
    shell:
        """
        plink2 --bfile {params.bfile} --pca {params.comp} --out results/{wildcards.run_id}/{wildcards.group}/common/pca --threads {threads}
        """


# Clean PCA eigenvec file for use as covariates
rule clean_pca_eigenvec:
    input:
        pca=rules.pca.output.pca,
    output:
        covar="results/{run_id}/{group}/common/qcovar",
    conda:
        "../envs/base.yml"
    shell:
        "sed '1s/^#//' {input} > {output}"
