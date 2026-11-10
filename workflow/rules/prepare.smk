# Create sample list and phenotype file for each phenotype
rule create_sample_list:
    output:
        sample_list="results/{run_id}/{phenotype}/common/sample.list",
        phenotype="results/{run_id}/{phenotype}/common/phenotype",
    log:
        "logs/{run_id}/{phenotype}/create_sample_list.log",
    conda:
        "../envs/base.yml"
    run:
        import pandas as pd

        phenotype = pd.read_csv(config["phenotype"], sep="\t")[
            ["FID", "IID", wildcards.phenotype]
        ].dropna()
        phenotype.to_csv(output.phenotype, sep="\t", index=False)
        phenotype.iloc[:, :2].to_csv(output.sample_list, sep="\t", index=False)


# Extract phenotyped samples from main bfile
rule extract_bed_step1:
    input:
        sample_list=rules.create_sample_list.output.sample_list,
    output:
        bfile=temp(
            multiext(
                "results/{run_id}/{phenotype}/common/step1",
                ".bed",
                ".bim",
                ".fam",
            )
        ),
    log:
        "logs/{run_id}/{phenotype}/step1_plink.log",
    conda:
        "../envs/plink2.yml"
    params:
        step1=config["bfile"]["step1"],
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/common/step1",
    shell:
        """
        plink2 --bfile {params.step1} --keep {input.sample_list} --maf 0.01 --geno 0.1 --out {params.output_prefix} --make-bed --threads 1 &>{log}
        """


rule extract_bed_step2:
    input:
        sample_list=rules.create_sample_list.output.sample_list,
    output:
        bfile=temp(
            multiext(
                "results/{run_id}/{phenotype}/common/step2",
                ".bed",
                ".bim",
                ".fam",
            )
        ),
    log:
        "logs/{run_id}/{phenotype}/step2_plink.log",
    conda:
        "../envs/plink2.yml"
    params:
        step2=config["bfile"]["step2"],
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/common/step2",
    shell:
        """
        plink2 --bfile {params.step2} --keep {input.sample_list} --maf 0.001 --out {params.output_prefix} --make-bed --threads 1 &>{log}
        """


# Observed genome-wide heterozygosity, not the PLINK inbreeding coefficient F.
rule heterozygosity:
    input:
        bfile=rules.extract_bed_step2.output.bfile,
    output:
        report=temp("results/{run_id}/{phenotype}/common/het.het"),
        covar="results/{run_id}/{phenotype}/common/het.tsv",
    log:
        "logs/{run_id}/{phenotype}/heterozygosity.log",
    conda:
        "../envs/plink2.yml"
    params:
        bfile=rules.extract_bed_step2.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/common/het",
    shell:
        """
        plink2 --bfile {params.bfile} --het cols=fid,het,nobs \
            --out {params.output_prefix} --threads 1 &>{log}
        # Columns: FID, IID, O(HET), E(HET), OBS_CT.
        awk 'BEGIN {{ OFS="\\t"; print "FID", "IID", "HET" }}
             NR > 1 {{
                 if ($5 <= 0) {{
                     print "No observed genotypes for " $1 " " $2 > "/dev/stderr"
                     exit 1
                 }}
                 print $1, $2, sprintf("%.10g", $3 / $5)
             }}' {output.report} >{output.covar} 2>>{log}
        """


# Perform PCA analysis for covariates
rule pca:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        temp(
            multiext(
                "results/{run_id}/{phenotype}/common/pca",
                pca=".eigenvec",
                eval=".eigenval",
            )
        ),
    log:
        "logs/{run_id}/{phenotype}/pca.log",
    conda:
        "../envs/plink2.yml"
    threads: config["plink2"]["pca_threads"]
    resources:
        cpu_per_task=threads,
    params:
        bfile=rules.extract_bed_step1.params.output_prefix,
        comp=config["pca"],
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.phenotype}/common/pca",
    shell:
        """
        plink2 --bfile {params.bfile} --pca {params.comp} --out {params.output_prefix} --threads {threads} &>{log}
        """


# Clean PCA eigenvec file for use as covariates
rule clean_pca_eigenvec:
    input:
        pca=rules.pca.output.pca,
    output:
        covar=temp("results/{run_id}/{phenotype}/common/qcovar"),
    log:
        "logs/{run_id}/{phenotype}/clean_pca.log",
    conda:
        "../envs/base.yml"
    shell:
        "sed '1s/^#//' {input} > {output} 2> {log}"
