# GWASnake — QTL fixed effects

Estimate additive and dominance QTL effects, conditional on other QTL and
background genetic relatedness.

1. Set phenotype, genotype, and QTL inputs in [`config/config.yaml`](config/config.yaml).
2. Install and run:

```bash
pixi install
pixi run snakemake --dry-run
pixi run snakemake --sdm conda --cores 8
```

Output: `results/{run_id}/{phenotype}/conditional/qtl_effects.tsv`.
See [`config/README.md`](config/README.md) for input formats.

Author: RuLei Chen. Please cite this repository when using the workflow.
