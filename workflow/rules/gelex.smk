wildcard_constraints:
    phenotype="[^_]+",
    model="add|dom|joint",
    grm_type="add|dom",


def _get_transform(wildcards):
    if wildcards.phenotype in config.get("transform", {}):
        return f"--transform {config['transform'][wildcards.phenotype]}"
    return ""


def _grm_files(wildcards, grm_type):
    prefix = f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/grm"
    grm = [f"{prefix}.{grm_type}.{ext}" for ext in ("bin", "id")]
    loco = [
        f"{prefix}.{grm_type}.chr{c}.{ext}"
        for c in range(1, 13)
        for ext in ("bin", "id")
    ]
    return grm + loco


def _assoc_input(wildcards):
    base = {
        "phenotype": rules.create_sample_list.output.phenotype,
        "qcovar": rules.clean_pca_eigenvec.output.covar,
        "bfile": rules.extract_bed_step2.output.bfile,
        "grm_add": _grm_files(wildcards, "add"),
    }
    if wildcards.model in ("dom", "joint"):
        base["grm_dom"] = _grm_files(wildcards, "dom")
    return base


def _assoc_flags(wildcards):
    grm_prefix = f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/grm"
    if wildcards.model == "add":
        return f"--grm {grm_prefix}.add"
    if wildcards.model == "dom":
        return f"--grm {grm_prefix}.add {grm_prefix}.dom --model d"
    return f"--grm {grm_prefix}.add {grm_prefix}.dom --test joint"


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
    conda:
        "../envs/gelex.yml"
    params:
        bfile_prefix=rules.extract_bed_step1.params.output_prefix,
        output_prefix=lambda wildcards: f"results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/grm",
    shell:
        """
        gelex grm -b {params.bfile_prefix} --{wildcards.grm_type} -o {params.output_prefix} -t {threads}
        gelex grm -b {params.bfile_prefix} --{wildcards.grm_type} -o {params.output_prefix} -t {threads} --loco
        """


rule gelex_assoc:
    input:
        unpack(_assoc_input),
    output:
        assoc="results/{run_id}/{group}/{phenotype}/{model}_assoc.gwas.tsv",
    threads: config["gelex"]["assoc_threads"]
    resources:
        cpus_per_task=threads,
    conda:
        "../envs/gelex.yml"
    params:
        bfile_prefix=rules.extract_bed_step2.params.output_prefix,
        flags=_assoc_flags,
        transform=_get_transform,
    shell:
        """
        gelex assoc -b {params.bfile_prefix} -p {input.phenotype} {params.flags} {params.transform} --qcovar {input.qcovar} -o results/{wildcards.run_id}/{wildcards.group}/{wildcards.phenotype}/{wildcards.model}_assoc -t {threads} --loco
        """


rule plot_gelex:
    input:
        summary="results/{run_id}/{group}/{phenotype}/{model}_assoc.gwas.tsv",
    output:
        png="results/{run_id}/{group}/{phenotype}/{model}_manhattan.png",
    conda:
        "../envs/base.yml"
    script:
        "../scripts/plot_gelex_result.py"
