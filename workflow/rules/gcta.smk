# Create full genome GRM (Genetic Relationship Matrix)
wildcard_constraints:
    phenotype="[^_]+",


rule gcta_grm_full:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        grm=temp(
            multiext(
                "results/{run_id}/{group}/gcta/grm_full",
                ".grm.bin",
                ".grm.id",
                ".grm.N.bin",
            )
        ),
    threads: config["gcta"]["grm_threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile_prefix=rules.extract_bed_step1.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/gcta/grm_full",
    shell:
        """
        gcta64 --bfile {params.bfile_prefix} --make-grm --out {params.output_prefix} --threads {threads}
        """


# Create per-chromosome GRM for leave-one-chromosome-out analysis
rule gcta_grm_chr:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        grm=temp(
            multiext(
                "results/{run_id}/{group}/gcta/grm_chr{chr}",
                ".grm.bin",
                ".grm.id",
                ".grm.N.bin",
            )
        ),
    threads: config["gcta"]["grm_threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile_prefix=rules.extract_bed_step1.params.output_prefix,
        chr_num=lambda wildcards: int(wildcards.chr),
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/gcta/grm_chr{wildcards.chr}",
    shell:
        """
        gcta64 --bfile {params.bfile_prefix} --make-grm --chr {params.chr_num} --out {params.output_prefix} --threads {threads}
        """


# Perform MLMA (Mixed Linear Model Association) with GRM subtraction
rule gcta_mlma:
    input:
        grm_full=rules.gcta_grm_full.output.grm,
        grm_chr=rules.gcta_grm_chr.output.grm,
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=rules.clean_pca_eigenvec.output.covar,
        covar=rules.create_sample_list.output.covar,
    output:
        assoc=temp("results/{run_id}/{group}/gcta/{phenotype}_{chr}.mlma"),
    threads: config["gcta"]["mlma_threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile_prefix=config["bfile"]["step2"],
        grm_full_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/gcta/grm_full",
        grm_chr_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/gcta/grm_chr{wildcards.chr}",
        chr_num=lambda wildcards: int(wildcards.chr),
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/gcta/{wildcards.phenotype}_{wildcards.chr}",
    shell:
        """
        gcta64 --mlma --grm {params.grm_full_prefix} --mlma-subtract-grm {params.grm_chr_prefix} --bfile {params.bfile_prefix} --chr {params.chr_num} --pheno {input.phenotype} --out {params.output_prefix} --thread-num {threads} --qcovar {input.qcovar} --covar {input.covar}
        """


# Merge chromosome-specific MLMA results into a single file
rule merge_gcta_results:
    input:
        lambda wildcards: expand(
            "results/{run_id}/{group}/gcta/{phenotype}_{chr}.mlma",
            run_id=wildcards.run_id,
            phenotype=wildcards.phenotype,
            group=wildcards.group,
            chr=range(1, config["gcta"]["chromosomes"] + 1),
        ),
    output:
        merged="results/{run_id}/{group}/gcta/{phenotype}.mlma",
    conda:
        "../envs/base.yml"
    shell:
        """
        head -n 1 {input[0]} > {output.merged}
        for file in {input}; do
            tail -n +2 "$file" >> {output.merged}
        done
        """


# Generate Manhattan plot from GCTA MLMA results
rule plot_gcta:
    input:
        summary="results/{run_id}/{group}/gcta/{phenotype}.mlma",
    output:
        png="results/{run_id}/{group}/gcta/{phenotype}.png",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/plot_gcta_result.py"
