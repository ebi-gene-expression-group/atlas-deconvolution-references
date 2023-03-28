# atlas-deconvolution-references
Snakemake workflow to generate single-cell reference dataset for deconvolution in Atlas

It contains scripts rules for:
- SeuratObject creation
- Spliting SeuratObjects into tissues
- Reducing cell type label granularity
- Creation of references for deconvolution tools
- UMAP plots for quality control
- Creating summary file of references 

## Prerequisites

 * Snakemake
 * LSF batch scheduler
 * Set up configuration variables at `run_deconvolution_reference_creation.sh`.


## Run pipeline

```
run_deconvolution_reference_creation.sh
```
