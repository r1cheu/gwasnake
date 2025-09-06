rule create_sample_list:
    output:
        sample_list="results/{run_id}/{group}/sample.list",
        phenotype="results/{run_id}/{group}/phenotype",
    conda:
        "../envs/base.yml"
    params:
        halfsib=HALFSIB,
    log:
        "logs/{run_id}/{group}/create_sample_list.log",
    script:
        "../scripts/create_sample_list.py"


rule extract_bed:
    input:
        sample_list=rules.create_sample_list.output.sample_list,
    output:
        bfile=multiext("results/{run_id}/{group}/step1", ".bed", ".bim", ".fam"),
    conda:
        "../envs/plink.yml"
    params:
        step1=config["bfile"]["step1"],
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/step1",
    shell:
        """
        plink --bfile {params.step1} --keep {input.sample_list} --maf 0.01 --out {params.output_prefix} --make-bed --threads 1
        """


rule pca:
    input:
        bfile=rules.extract_bed.output.bfile,
    output:
        pca=temp("results/{run_id}/{group}/pca.eigenvec"),
        eval=temp("results/{run_id}/{group}/pca.eigenval"),
    conda:
        "../envs/plink2.yml"
    threads: 2
    resources:
        cpu_per_task=threads,
    params:
        bfile=rules.extract_bed.params.output_prefix,
    log:
        "logs/{run_id}/{group}/pca.log",
    shell:
        """
        plink2 --bfile {params.bfile} --pca 20 --out results/{wildcards.run_id}/{wildcards.group}/pca --threads 2
        """


rule clean_pca_eigenvec:
    input:
        pca=rules.pca.output.pca,
    output:
        covar="results/{run_id}/{group}/covar",
    conda:
        "../envs/base.yml"
    shell:
        "sed '1s/^#//' {input} > {output}"
