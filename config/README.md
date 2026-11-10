# QTL-region variance

Configure `config.yaml` with:

- `phenotype`: TSV headed by `FID`, `IID`, and trait names; `phenotypes` selects traits.
- `bfile.step2`: PLINK prefix. No step1 dataset or PCA is needed.
- `qtl_list`: TSV columns `CHR SNP BP A1 A2 START END PHENOTYPE`, in that order.
- `qtl_scope`: `per_phenotype` uses the current trait's QTL; `all` uses all loci
  and excludes all input regions from the background.
- `ld_prune`: two background LD-pruning settings.
- `transform`: `none`, `dint`, or `iint`.
- `run_id` and `gelex`: output identifier and thread counts.

Only listed traits with QTL are targeted. Loci need all three genotype classes.
QTL GRMs use SNPs in retained loci's regions; background GRMs exclude those
regions (all input regions for `qtl_scope: all`). Both use NOIA standardized encoding.

Output: `results/{run_id}/{phenotype}/variance/variance.tsv`.

Plot one or more result tables:

```bash
pixi run python workflow/scripts/plot_grm_variance.py \
    results/RUN/*/variance/variance.tsv -o results/RUN/variance.svg
```
