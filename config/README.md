# GWASnake Configuration

This directory contains configuration files for the GWASnake workflow - a Snakemake workflow for genome-wide association study (GWAS) analysis.

## Configuration Files

### config/config.yaml

Main workflow configuration file containing:

- **run_id**: Identifier for the current analysis run (e.g., "oh_ho")
- **halfsib**: Path to sample grouping file
- **phenotype**: Path to phenotype data file
- **bfile**: Path to PLINK binary files directory

### config/sample.csv

Sample grouping file with the following structure:

|                         |
| ----------------------- |
| sample1,sample2,sample3 |
| sample4,sample5,sample6 |

Each row defines a sample group for parallel processing.

## Input Data Requirements

### PLINK Binary Files

Required files in the bfile/ directory:

- .bed: Binary genotype data
- .bim: Variant information
- .fam: Sample information

### Phenotype File

Tab-separated file with columns, header is needed:

- FID: Family ID
- IID: Individual ID
- phenotype: Phenotype values (e.g., heading_date)

## Output Structure

Results are organized by run_id and group:

```
results/
  {run_id}/
    {group}/
      common/
        sample.list          # Filtered sample IDs
        phenotype            # Filtered phenotype data
      regenie/               # regenie output files
      gcta/                  # gcta output files
```

## Customization

### Adding New Groups

1. Add group definition to config/sample.csv
2. Update sample lists as comma-separated values

### Changing Phenotype Analysis

1. Update phenotype file path in config.yaml
2. Ensure phenotype column name matches your data

## Running the Workflow

See the main project README for execution instructions. The workflow automatically uses configurations from this directory.
