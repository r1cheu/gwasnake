# QTL fixed effects

Configure `config.yaml` with:

- `phenotype`: TSV headed by `FID`, `IID`, and trait names; `phenotypes` selects traits.
- `bfile.step2`: PLINK prefix. No step1 dataset or PCA is needed.
- `qtl_list`: TSV columns `CHR SNP BP A1 A2 START END PHENOTYPE`, in that order.
- `qtl_scope`: `per_phenotype` uses the current trait's QTL; `all` uses all loci
  and excludes all input regions from the background.
- `ld_prune`: two background LD-pruning settings.
- `transform`: `none` or `dint`; `iint` removes the QTL fixed effects.
- `run_id` and `gelex`: output identifier and thread counts.

Only listed traits with QTL are targeted. Loci need all three genotype classes;
collinear loci are dropped. QTL alleles must match the BIM orientation.

```bash
pixi install
pixi run snakemake --sdm conda --cores 8
```

Output: `results/{run_id}/{phenotype}/conditional/qtl_effects.tsv`.
