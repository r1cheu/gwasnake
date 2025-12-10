# Create full genome GRM (Genetic Relationship Matrix)
wildcard_constraints:
    phenotype="[^_]+",


rule create_step1_tfam:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        tfile=temp(
            multiext(
                "results/{run_id}/{group}/{phenotype}/emmax/step1",
                ".tfam",
                ".tped",
            )
        ),
    conda:
        "../envs/plink.yml"
    threads: config["plink"]["threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile_prefix=rules.extract_bed_step1.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/emmax/step1",
    shell:
        """
        plink --bfile {params.bfile_prefix} --recode 12 transpose --output-missing-genotype 0 --out {params.output_prefix}
        """


rule create_step2_tfam:
    input:
        bfile=rules.extract_bed_step2.output.bfile,
    output:
        tfile=temp(
            multiext(
                "results/{run_id}/{group}/{phenotype}/emmax/step2",
                ".tfam",
                ".tped",
            )
        ),
    conda:
        "../envs/plink.yml"
    threads: config["plink"]["threads"]
    resources:
        cpus_per_task=threads,
    params:
        bfile_prefix=rules.extract_bed_step2.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/emmax/step2",
    shell:
        """
        plink --bfile {params.bfile_prefix} --recode 12 transpose --output-missing-genotype 0 --out {params.output_prefix}
        """


rule emmax_create_pca:
    input:
        pca=rules.pca.output.pca,
    output:
        add_col3=temp("results/{run_id}/{group}/{phenotype}/emmax/emmax_fit_pca.tmp"),
        pca_emma="results/{run_id}/{group}/{phenotype}/emmax/emmax_fit_pca.eigenvec",
    conda:
        "../envs/base.yml"
    shell:
        """        
        awk 'BEGIN{{OFS="\t"}}{{for(i=NF+1;i>3;i--)$i=$(i-1);$3="1";print}}' {input.pca} > {output.add_col3}
        grep -v "#FID" {output.add_col3} >{output.pca_emma}    
        """


# Create per-chromosome GRM for leave-one-chromosome-out analysis
rule emmax_kin:
    input:
        tfile=rules.create_step1_tfam.output.tfile,
    output:
        IBS_kinf=temp("results/{run_id}/{group}/{phenotype}/emmax/step1.aIBS.kinf"),
    threads: config["emmax"]["threads"]
    resources:
        cpus_per_task=threads,
    params:
        tfile_prefix=rules.create_step1_tfam.params.output_prefix,
    shell:
        """
        emmax-kin-intel64 -v -s -d 10 {params.tfile_prefix}
        """


rule emmax_pheno:
    input:
        phenofile=rules.create_sample_list.output.phenotype,
        fam_file="results/{run_id}/{group}/{phenotype}/common/step1.fam",
    output:
        emmax_pheno=temp(
            "results/{run_id}/{group}/{phenotype}/emmax/{phenotype}_emmax_fit.txt"
        ),
    conda:
        "../envs/base.yml"
    params:
        awk_cmd1=r"'{print $1}'",
        awk_cmd2=r"-v value=$id '{if($1==value)print $0}'",
    shell:
        """
        awk {params.awk_cmd1} {input.fam_file} |while read id;do awk {params.awk_cmd2} {input.phenofile};done > {output.emmax_pheno}
        """


# Perform MLMA (Mixed Linear Model Association) with GRM subtraction
rule emmax_mlm:
    input:
        tfile=rules.create_step2_tfam.output.tfile,
        phenofile=rules.emmax_pheno.output.emmax_pheno,
        qcovar=rules.emmax_create_pca.output.pca_emma,
        kin=rules.emmax_kin.output.IBS_kinf,
    output:
        assoc="results/{run_id}/{group}/{phenotype}/emmax/assoc.ps",
    threads: config["emmax"]["threads"]
    resources:
        cpus_per_task=threads,
    params:
        tfile_prefix=rules.create_step2_tfam.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/emmax/assoc",
    shell:
        """
        emmax-intel64 -v -d 10 -t {params.tfile_prefix} -p {input.phenofile} -k {input.kin} -c {input.qcovar} -o {params.output_prefix}
        """


# Generate Manhattan plot from GCTA MLMA results
rule plot_emmax:
    input:
        summary="results/{run_id}/{group}/{phenotype}/emmax/assoc.ps",
    output:
        png="results/{run_id}/{group}/{phenotype}/emmax/manhattan.png",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/plot_emmax_result.py"
