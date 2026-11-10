# GWAS configuration

`config/config.yaml` configures the GWAS workflow on `main`.

| Setting | Purpose |
| --- | --- |
| `run_id` | Output directory identifier; use a new value for a separate analysis. |
| `phenotype` | Tab-separated phenotype file, with a header starting with `FID` and `IID`. |
| `phenotypes` | Trait column names to analyze. |
| `bfile.step1` | PLINK prefix for LD-filtered markers used for PCA and GRMs. |
| `bfile.step2` | PLINK prefix for markers tested in association. |
| `pca` | Number of principal components used as covariates. |
| `transform` | Association phenotype transformation: `iint`, `dint`, or `none`. |
| `plink2.pca_threads` | Threads for PCA. |
| `gelex.grm_threads` | Threads for GRM construction. |
| `gelex.assoc_threads` | Threads for association. |

Both PLINK prefixes must have `.bed`, `.bim`, and `.fam` files. For each
phenotype, samples with missing phenotype values are excluded before genotype
extraction.

The workflow builds additive and dominance GRMs with NOIA standardized encoding
(`NS`), then runs joint additive/dominance LOCO association with NOIA centered
encoding (`NC`). Both the global and chromosome-specific GRMs use `step1`.
The current LOCO output declaration expects numeric chromosomes 1–12, grouped
by chromosome, with markers remaining on each chromosome after filtering.

The default target produces three Manhattan plots per phenotype:

```text
results/{run_id}/{phenotype}/
├── common/              # Sample selection, extracted genotypes, and PCA
├── assoc.gwas.tsv       # Joint additive/dominance association results
├── manhattan_A.png
├── manhattan_D.png
└── manhattan_AD.png
```

Intermediate genotype, PCA, and GRM files are marked temporary and may be removed
after downstream jobs complete. Logs are under `logs/{run_id}/{phenotype}/`.

No QTL list, sample-group file, or analysis selector is required.
