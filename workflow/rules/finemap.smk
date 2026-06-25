wildcard_constraints:
    grm_type="add|dom",
    chr=r"\d+",
    region=r"[^/]+",
    phenotype=r"[^/]+",
    run_id=r"[^/]+",
    group=r"[^/]+",


def _transform(wildcards):
    return config.get("transform", {}).get(wildcards.phenotype, "none")


def _grm_files(run_id, group, phenotype, grm_type):
    prefix = f"results/{run_id}/{group}/{phenotype}/grm"
    files = [f"{prefix}.{grm_type}.{ext}" for ext in ("bin", "id")]
    files += [
        f"{prefix}.{grm_type}.chr{c}.{ext}"
        for c in range(1, 13)
        for ext in ("bin", "id")
    ]
    return files


def _reml_qcovar_flag(wildcards):
    base = f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}"
    # iint already absorbed the PCs into the residual; dint/none keep them as fixed effects
    if _transform(wildcards) == "iint":
        return ""
    return f"--qcovar {base}/common/qcovar"


def _susie_input(wildcards):
    chrom = REGION_INFO[(wildcards.phenotype, wildcards.region)][0]
    base = f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}"
    return {
        "geno": f"{base}/finemap/{wildcards.region}/geno.raw",
        "linv": f"{base}/finemap/Linv.chr{chrom}.rds",
        "phenotype": f"{base}/common/transformed.phen",
        "qcovar": f"{base}/common/qcovar",
        "bim": f"{base}/common/step2.bim",
    }


# Build additive/dominance GRMs (whole-genome + per-chromosome LOCO)
rule gelex_grm:
    input:
        bfile=rules.extract_bed_step1.output.bfile,
    output:
        grm=temp(
            multiext(
                "results/{run_id}/{group}/{phenotype}/grm",
                ".{grm_type}.bin",
                ".{grm_type}.id",
            )
        ),
        loco=temp(
            expand(
                "results/{{run_id}}/{{group}}/{{phenotype}}/grm.{{grm_type}}.chr{chr}.{ext}",
                chr=range(1, 13),
                ext=["bin", "id"],
            )
        ),
    threads: config["gelex"]["grm_threads"]
    resources:
        cpus_per_task=threads,
    log:
        "logs/{run_id}/{group}/{phenotype}/grm_{grm_type}.log",
    params:
        bfile_prefix=rules.extract_bed_step1.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/grm",
    shell:
        """
        gelex grm -b {params.bfile_prefix} --{wildcards.grm_type} -o {params.output_prefix} -t {threads} &> {log}
        gelex grm -b {params.bfile_prefix} --{wildcards.grm_type} -o {params.output_prefix} -t {threads} --loco &>> {log}
        """


# Single source of truth for the RINT-transformed phenotype: shared by reml and susie
rule transform_phenotype:
    input:
        phenotype=rules.create_sample_list.output.phenotype,
        qcovar=rules.clean_pca_eigenvec.output.covar,
    output:
        phenotype="results/{run_id}/{group}/{phenotype}/common/transformed.phen",
    params:
        transform=_transform,
    conda:
        "../envs/susier.yml"
    log:
        "logs/{run_id}/{group}/{phenotype}/transform_phenotype.log",
    script:
        "../scripts/transform_phenotype.R"


# Estimate per-chromosome variance components (LOCO) on the transformed phenotype
rule gelex_reml:
    input:
        phenotype=rules.transform_phenotype.output.phenotype,
        qcovar=rules.clean_pca_eigenvec.output.covar,
        grm_add=lambda w: _grm_files(w.run_id, w.group, w.phenotype, "add"),
        grm_dom=lambda w: _grm_files(w.run_id, w.group, w.phenotype, "dom"),
    output:
        summary="results/{run_id}/{group}/{phenotype}/reml.loco.summary",
    threads: config["gelex"]["reml_threads"]
    resources:
        cpus_per_task=threads,
    log:
        "logs/{run_id}/{group}/{phenotype}/reml.log",
    params:
        grm_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/grm",
        out_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/reml",
        qcovar_flag=_reml_qcovar_flag,
    shell:
        """
        gelex reml -p {input.phenotype} --grm {params.grm_prefix}.add {params.grm_prefix}.dom {params.qcovar_flag} --transform none --loco -o {params.out_prefix} -t {threads} &> {log}
        """


# Build V = sa2*G_add_loco + sd2*G_dom_loco + se2*I for one chromosome, save whitening L^-1
rule build_whitening:
    input:
        grm_add=multiext(
            "results/{run_id}/{group}/{phenotype}/grm.add", ".bin", ".id"
        ),
        grm_add_chr=multiext(
            "results/{run_id}/{group}/{phenotype}/grm.add.chr{chr}", ".bin", ".id"
        ),
        grm_dom=multiext(
            "results/{run_id}/{group}/{phenotype}/grm.dom", ".bin", ".id"
        ),
        grm_dom_chr=multiext(
            "results/{run_id}/{group}/{phenotype}/grm.dom.chr{chr}", ".bin", ".id"
        ),
        summary=rules.gelex_reml.output.summary,
    output:
        linv=temp("results/{run_id}/{group}/{phenotype}/finemap/Linv.chr{chr}.rds"),
    conda:
        "../envs/susier.yml"
    log:
        "logs/{run_id}/{group}/{phenotype}/build_whitening_chr{chr}.log",
    params:
        chr="{chr}",
    script:
        "../scripts/build_whitening.R"


# Extract a region's genotypes as additive dosages
rule extract_region:
    input:
        bfile=rules.extract_bed_step2.output.bfile,
    output:
        geno=temp("results/{run_id}/{group}/{phenotype}/finemap/{region}/geno.raw"),
    conda:
        "../envs/plink.yml"
    log:
        "logs/{run_id}/{group}/{phenotype}/extract_region_{region}.log",
    params:
        bfile_prefix=rules.extract_bed_step2.params.output_prefix,
        out_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/finemap/{wildcards.region}/geno",
        chr=lambda wildcards: REGION_INFO[(wildcards.phenotype, wildcards.region)][0],
        start=lambda wildcards: REGION_INFO[(wildcards.phenotype, wildcards.region)][1],
        end=lambda wildcards: REGION_INFO[(wildcards.phenotype, wildcards.region)][2],
    shell:
        """
        plink --bfile {params.bfile_prefix} --chr {params.chr} --from-bp {params.start} --to-bp {params.end} --recode A --keep-allele-order --out {params.out_prefix} --threads 1 &> {log}
        """


# Whiten X (NC-encoded add+dom) and y, run susieR fine-mapping
rule run_susie:
    input:
        unpack(_susie_input),
    output:
        pip="results/{run_id}/{group}/{phenotype}/finemap/{region}/pip.tsv",
        cs="results/{run_id}/{group}/{phenotype}/finemap/{region}/credible_sets.tsv",
        rds="results/{run_id}/{group}/{phenotype}/finemap/{region}/susie.rds",
    conda:
        "../envs/susier.yml"
    log:
        "logs/{run_id}/{group}/{phenotype}/run_susie_{region}.log",
    params:
        transform=_transform,
        susie=config["susie"],
    script:
        "../scripts/run_susie.R"
