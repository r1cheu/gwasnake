# GWASnake — QTL heterogeneity

Test whether QTL dominance effects vary across genetic-background strata using
PCA clustering and null/full REML models.

1. Set phenotype, genotype, and QTL inputs in [`config/config.yaml`](config/config.yaml).
2. Install and run:

```bash
pixi install
pixi run snakemake --dry-run
pixi run snakemake --sdm conda --cores 8
```

Output: `results/{run_id}/{phenotype}/heterogeneity/heterogeneity.tsv`.
See [`config/README.md`](config/README.md) for input formats and gate settings.

Author: RuLei Chen. Please cite this repository when using the workflow.
