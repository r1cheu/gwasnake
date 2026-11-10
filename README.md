# GWASnake — QTL-region variance

Estimate additive and dominance variance components using separate GRMs for
QTL regions and the remaining genetic background.

1. Set phenotype, genotype, and QTL inputs in [`config/config.yaml`](config/config.yaml).
2. Install and run:

```bash
pixi install
pixi run snakemake --dry-run
pixi run snakemake --sdm conda --cores 8
```

Output: `results/{run_id}/{phenotype}/variance/variance.tsv`.
See [`config/README.md`](config/README.md) for input formats and plotting.

Author: RuLei Chen. Please cite this repository when using the workflow.
