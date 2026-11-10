# Select phenotyped samples and extract their genotypes.
rule create_sample_list:
    input:
        phenotype=config["phenotype"],
    output:
        sample_list="results/{run_id}/{phenotype}/common/sample.list",
        phenotype="results/{run_id}/{phenotype}/common/phenotype",
    log:
        "logs/{run_id}/{phenotype}/create_sample_list.log",
    conda:
        "../envs/base.yml"
    run:
        import pandas as pd

        phenotype = pd.read_csv(input.phenotype, sep="\t")[
            ["FID", "IID", wildcards.phenotype]
        ].dropna()
        phenotype.to_csv(output.phenotype, sep="\t", index=False)
        phenotype.iloc[:, :2].to_csv(output.sample_list, sep="\t", index=False)


rule extract_bed_step2:
    input:
        sample_list=rules.create_sample_list.output.sample_list,
        bfile=multiext(config["bfile"]["step2"], ".bed", ".bim", ".fam"),
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
