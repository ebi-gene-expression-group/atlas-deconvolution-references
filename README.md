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
 * A batch scheduler, LSF or SLURM
 * Set up configuration variables at `run_deconvolution_reference_creation.sh`
 * Config.yaml with sc experiment accessions as input and tissues that are excluded for accessions (look at example.yaml how that should look like)
 * Download ontology files for [UBERON](http://purl.obolibrary.org/obo/uberon/basic.obo) and [CL](http://purl.obolibrary.org/obo/cl/c-basic.obo) ontology and put in folder `files`
 
## Run pipeline

```
bash run_deconvolution_reference_creation.sh
```
This file contains:
* definition of variables
* runs snakemake like this:
```
snakemake \
  -s $SNAKEFILE  \
  --profile [cluster profile 'lsf' or 'slurm'] \
  --rerun-incomplete \
  --retries 3 \ # if job fails try with more memory
  --use-conda --conda-frontend mamba --keep-going \
  --config atlas_prod=[path/to/atlas_prot_dir]\
  anndata_prod=[path/to/anndata_prot_dir]\
  deconv_ref=[path/to/desination/of/reference_library] \
```
For SLURM profile, add `--slurm`.

* checks if references for organims parts are generated more than once and returns duplicated organism parts which can then be included in config.yaml.

example:
For these UBERON ids multiple references exist:
-------------------------------------------- 
UBERON_0001159

If this message is shown it is neccessary to go look at /UMAP/UBERON_0001159_*.png and decide 
Then either:
 - delete in $DECONV_REF the three references that belong to the reference we do not want the reference from and include UBERON_id_scxa_accession in config.yaml file under "exclude_tissues_from_accessions:" OR
 - decide which reference to delete, add them to config.yaml and reurn whole generation of deconvolution reference library
 
## Output 

The output of this pipeline are :

* deconvolution reference lib with three files per scxa accession + organism part (CO_scaled, C1 and phenData)
* one summary file of all the references created
* one png file per scxa accession + organism part with two UMAP plots (old + new cell type labels) for quality control

## Contributions
The following scripts are taken from the [CATD_snakemake pipeline](https://github.com/Functional-Genomics/CATD_snakemake) which benchmarks deconvolution methods:
* scaleTable.R
* genRef_1.R
