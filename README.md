# Snakemake workflow: `gwasnake`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

A Snakemake workflow for `GWAS in NCII population`

- [Snakemake workflow: `gwasnake`](#snakemake-workflow-gwasnake)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## Usage

Detailed information about input data and workflow configuration can also be found in the [`config/README.md`](config/README.md).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Deployment options

To run the workflow from command line, change the working directory. This workflow need to run with **conda**.

```bash
cd gwasnake
```

Place your plink bfile in the `bfile/` directory and your phenotype file in the `phenotype/` directory. The exact paths are flexible, as long as they match the entries in `config/config.yaml`. The default configuration file is located at `config/config.yml`.

Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

You can prepare the environments with:

```bash

snakemake --sdm conda --conda-create-envs-only

```

Run the workflow with 20 cores:

```bash
snakemake --cores 20 --sdm conda --keep-going
```

`--keep-going` will continue running independent jobs even after a job fails, which is useful if jobs for small populations may fail in `regenie`.

It's recommended to run the workflow with slurm, and do not forget to change the account in the `slurm/config.yaml` file if you are using slurm.

```bash
snakemake --sdm conda --profile slurm --keep-going

```

Use the `--ri` to rerun the uncompleted jobs:

```bash
snakemake --sdm conda --profile slurm --ri

```

## Authors

- RuLei Chen
  - Develop the snakemake workflow
  - Center for Excellence in Molecular Plant Sciences
  - [ORCID](https://orcid.org/0009-0000-9645-0951)
  - [github](https://github.com/r1cheu)

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
