# atlas-deconvolution-references
Snakemake workflow to generate single-cell reference dataset for bulk RNA-seq deconvolution in Atlas

It contains rules for:
- SeuratObject creation
- Spliting SeuratObjects into tissues
- Reducing cell type label granularity
- Creation of references for deconvolution tools
- UMAP plots for quality control
- Creating summary file of references 

## Prerequisites

 * Snakemake
 * LSF batch scheduler
 * Set up configuration variables at `run_deconvolution_reference_creation.sh`
 * Config.yaml with sc experiment accessions as input and tissues that are excluded for accessions (look at example.yaml how that should look like)


## Run pipeline

```
bash run_deconvolution_reference_creation.sh
```
This file contains:
* definition of variables
* runs snakemake
```
snakemake \
  -s $SNAKEFILE  \
  --profile [cluster profile] \
  --rerun-incomplete \
  --retries 3 \ # if job fails try with more memory
  --use-conda --conda-frontend mamba --keep-going \
  --config atlas_prod=[path/to/atlas_prot_dir]\
  anndata_prod=[path/to/anndata_prot_dir]\
  deconv_ref=[path/to/desination/of/reference_library] \
```
* checks if references for organims parts are generated more than one and returns duplicated organism parts which can then be included in conig.yaml

## Output 

The output of this pipeline are :

* deconvolution reference lib with three files per scxa accession + organism part (CO_scaled, C1 and )
* one summary file of alle the references created
* one png file per scxa accession + organism part with two UMAP plots (old + new cell type labels) for quality control

## Contributions
Folowing scripts are taken from the [CATD_snakemake pipeline] (https://github.com/Functional-Genomics/CATD_snakemake)
