# Create full genome GRM (Genetic Relationship Matrix)
wildcard_constraints:
    phenotype="[^_]+",


rule gcta_grm_add:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        grm=temp(
            multiext(
                "results/{run_id}/{group}/{phenotype}/gcta/grm",
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
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/gcta/grm",
    shell:
        """
        gcta64 --bfile {params.bfile_prefix} --make-grm --out {params.output_prefix} --threads {threads}
        """


# Create per-chromosome GRM for leave-one-chromosome-out analysis
rule gcta_grm_dom:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        grm=temp(
            multiext(
                "results/{run_id}/{group}/{phenotype}/gcta/grm",
                ".d.grm.bin",
                ".d.grm.id",
                ".d.grm.N.bin",
            )
        ),
    threads: config["gcta"]["grm_threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile_prefix=rules.extract_bed_step1.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/gcta/grm",
    shell:
        """
        gcta64 --bfile {params.bfile_prefix} --make-grm-d --out {params.output_prefix} --threads {threads}
        """


rule create_grm_txt:
    input:
        grm_add=rules.gcta_grm_add.output.grm,
        grm_d=lambda wildcards: rules.gcta_grm_dom.output.grm if config["phenotype_setting"][wildcards.phenotype] == "dom" else [],
    output:
        txt="results/{run_id}/{group}/{phenotype}/gcta/grm.txt",
    params:
        prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/gcta/grm",
        phenotype_setting=lambda wildcards: config["phenotype_setting"][wildcards.phenotype],
    shell:
        """
        if [ "{params.phenotype_setting}" = "dom" ]; then
            echo '{params.prefix}\n{params.prefix}.d' > {output.txt}
        else
            echo '{params.prefix}' > {output.txt}
        fi
        """


# Perform MLMA (Mixed Linear Model Association) with GRM subtraction
rule gcta_mlma:
    input:
        grm_a=rules.gcta_grm_add.output.grm,
        grm_d=lambda wildcards: rules.gcta_grm_dom.output.grm if config["phenotype_setting"][wildcards.phenotype] == "dom" else [],
        grm=rules.create_grm_txt.output.txt,
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=rules.clean_pca_eigenvec.output.covar,
        bfile=rules.extract_bed_step2.output.bfile,
    output:
        assoc="results/{run_id}/{group}/{phenotype}/gcta/assoc.mlma",
    threads: config["gcta"]["mlma_threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/common/step2",
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/gcta/assoc",
    shell:
        """
        gcta64 --mlma --mgrm {input.grm} --bfile {params.bfile_prefix} --pheno {input.phenotype} --out {params.output_prefix} --thread-num {threads} --qcovar {input.qcovar}
        """


# Merge chromosome-specific MLMA results into a single file
# Generate Manhattan plot from GCTA MLMA results
rule plot_gcta:
    input:
        summary="results/{run_id}/{group}/{phenotype}/gcta/assoc.mlma",
    output:
        png="results/{run_id}/{group}/{phenotype}/gcta/manhattan.png",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/plot_gcta_result.py"
