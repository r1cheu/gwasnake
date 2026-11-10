# QTL heterogeneity

Configure `config.yaml` with:

- `phenotype`: TSV headed by `FID`, `IID`, and trait names; `phenotypes` selects traits.
- `bfile.step1`: LD-filtered PLINK prefix for PCA; `bfile.step2`: QTL and background markers.
- `qtl_list`: TSV columns `CHR SNP BP A1 A2 START END PHENOTYPE`, in that order.
- `qtl_scope`: `per_phenotype` uses the current trait's QTL; `all` uses all loci
  and excludes all input regions from the background.
- `ld_prune`: two background LD-pruning settings; `pca`: number of PCs for clustering.
- `transform`: `none` or `dint`, retaining QTL fixed effects in both models.
- `heterogeneity`: stratum count, minimum heterozygotes (`het_min`), minimum of
  each homozygote class (`hom_min`), and minimum informative strata (`min_strata`).
- `run_id`, `plink2`, and `gelex`: output identifier and thread counts.

Only listed traits with QTL are targeted. A checkpoint selects informative loci;
its downstream jobs appear after the gate runs. The null model includes additive
stratum slopes; the full model also includes dominance slopes. Results use a
boundary likelihood-ratio test and BH-FDR correction.

Output: `results/{run_id}/{phenotype}/heterogeneity/heterogeneity.tsv`.
Strata and per-locus results are kept alongside it. If no loci qualify, the final
table contains only its header.
