## Reference library preperation for deconvolution in 
## Expression Atlas
##
## Description:
##
## Gets scxa accessions from config.yaml file and creates tissue
## specifc references per accession that are used for deconvolution analysis
## 
##
##

#Set config here
# config should 'accessions' to devoncolve and 'exclude_tissues_from_accessions'
configfile: 'accession_deconvolution_references.yaml'

import glob
import os
import pandas as pd

def get_tissues_per_accession(wildcards):
    """
    Returns list of files with valid combinations of accession and organism
    parts that are acquired from cell_metadata.tsv
    """
    outnames = []
    for sp in config['species']:
        if sp['name'] == wildcards['species']:
            for accession in sp['accessions']:
                # ANND experiements not yet in sc_exps
                if "E-ANND-" in accession:
                    path = os.path.join(config['anndata_prod'], "*", accession, f"{accession}.cell_metadata.tsv")
                    path = glob.glob(path)
                    if not path:
                        print(f"Error: cell metadata not found for accession {accession}")
                        continue
                    try:
                        with open(path[-1], "r") as f:
                            data = pd.read_csv(f, sep="\t", low_memory=False)
                    except Exception as e:
                        print(f"Error: Failed to read file {path[0]}: {e}")
                        continue
                #for accessions in sc_exps
                else:
                    path = os.path.join(config['sc_exps'], accession, f"{accession}.cell_metadata.tsv")
                    if not os.path.isfile(path):
                        print(f"Error: cell metadata not found for accession {accession}")
                        continue
                    try:
                        with open(path, "r") as f:
                            data = pd.read_csv(f, sep="\t", low_memory=False)
                    except Exception as e:
                        print(f"Error: Failed to read file {path}: {e}")
                        continue
                # only select tissues which have cell type ontology labels
                try:
                    uberon_paths = data[data['inferred_cell_type_._ontology_labels_ontology'].notnull()]
                except KeyError:
                    data['inferred_cell_type_._ontology_labels_ontology'] = data['authors_cell_type_ontology']
                    uberon_paths = data[data['inferred_cell_type_._ontology_labels_ontology'].notnull()]
                uberon_paths = uberon_paths.groupby('organism_part_ontology')['inferred_cell_type_._ontology_labels_ontology'].nunique()
                uberon_paths = uberon_paths[uberon_paths >= 2].index.tolist()
                filtered_paths = []
                for path in uberon_paths:
                    cell_type_counts = data[data['organism_part_ontology'] == path]['inferred_cell_type_._ontology_labels_ontology'].value_counts()
                    if len(cell_type_counts[cell_type_counts >= 50]) >= 2:
                        filtered_paths.append(path)
                uberon_paths = filtered_paths
                # uberon_paths = list(set(data[data['inferred_cell_type_-_ontology_labels_ontology'].notnull()]["organism_part_ontology"]))
                # select only items that are UBERON path
                uberon_paths = [x for x in uberon_paths if "UBERON" in str(x)]
                # extract uberons from paths
                uberons = [os.path.basename(path) for path in uberon_paths]
                outnames.append([f"UMAP/{wildcards['species']}/{uberon}_{accession}_umap.png" for uberon in uberons])
    to_remove =  [f"UMAP/{wildcards['species']}/{uberon_and_accession}_umap.png" for uberon_and_accession in sp['exclude_tissues_from_accessions']]

    outnames = [item for sublist in outnames for item in sublist]
    # remove outnames where we dont want tissue references to be generated
    outnames = [x for x in outnames if x not in to_remove]
    return outnames

def input_for_copy(wildcards):
    """
    Returns input paths for rule copyInputFiles, neccessary as ANND accessions require different
    set of input files 
    """
    if "ANND" in wildcards['experiment']:
        mtx_path = os.path.join(config['anndata_prod'], "*", wildcards['experiment'])
        mtx_path = glob.glob(mtx_path)[-1]
        mtx = mtx_path + f"/matrices/raw/matrix.mtx.gz"
        col = mtx_path + f"/matrices/raw/barcodes.tsv.gz"
        row = mtx_path + f"/matrices/raw/genes.tsv.gz"
        meta = mtx_path + f"/{wildcards['experiment']}.cell_metadata.tsv"
        return [mtx,  col, row , meta]
    else:
        return [os.path.join(config['sc_exps'], wildcards['experiment'], \ 
         f"{wildcards['experiment']}.aggregated_filtered_counts.mtx"),
                os.path.join(config['sc_exps'], wildcards['experiment'],  \ 
                f"{wildcards['experiment']}.aggregated_filtered_counts.mtx_cols"), 
                os.path.join(config['sc_exps'], wildcards['experiment'], \ 
                 f"{wildcards['experiment']}.aggregated_filtered_counts.mtx_rows"), 
                os.path.join(config['sc_exps'], wildcards['experiment'],  \ 
                f"{wildcards['experiment']}.cell_metadata.tsv")]
            
def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules 
    attemps = reiterations + 1
    Max number attemps = 4
    """
    mem_avail = [16, 32, 64, 128, 300 ]  
    if attempt > len(mem_avail):
        print(f"Attemps {attempt} exceeds the maximum number of attemps: {len(mem_avail)}")
        print(f"modify value of --restart-times or adjust mem_avail resources accordingly")
        sys.exit(1)
    else:
        return mem_avail[attempt-1] * 1000

# Rule all
rule all:
    input:
        #get_tissues_per_accession()
        expand(config['deconv_ref'] + "/{species}_summary.tsv", species=[sp['name'] for sp in config['species']])
        
rule copyInputFiles:
    """
    Rule for copying all the mtx, gene, barcode and metadata files we want to build references from.
    """
    log: "logs/copyInputFiles/{species}/{experiment}.log"
    input:
       input_for_copy
    output:
       touch("scxa_input/{species}/{experiment}.copied")
    resources: mem_mb=get_mem_mb
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        mkdir -p scxa_input/{wildcards.species}
        cp {input} scxa_input/{wildcards.species}
        """
        
rule createSeuratObject:
    """
    Rule for creating SeuratObjects for anndata and non anndata scxa experiments.
    """
    conda: "envs/createSeuratObject.yaml"
    log: "logs/createSeuratObject/{species}_{experiment}.log"
    input:
        "scxa_input/{species}/{experiment}/{experiment}.copied"
    output:
        "scxa_input/{species}/{experiment}/{experiment}.project_seurat.rds"
    resources: mem_mb=get_mem_mb
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        if [[ {wildcards.experiment} == *"ANND"* ]]; then
            Rscript {workflow.basedir}/scripts/createSeuratObjectFromANND.R {wildcards.experiment} {wildcards.species}
        else
            Rscript {workflow.basedir}/scripts/createSeuratObject.R {wildcards.experiment} {wildcards.species}
        fi
        """

rule splitByTissue:
    """
    Rule to split sc reference into single organism parts if more than one
    organism part is present.
    """
    conda: "envs/TissueSplit.yaml"
    log: "logs/splitByTissue/{species}_{tissue}_{experiment}.log"
    input:
        "scxa_input/{species}/{experiment}/{experiment}.project_seurat.rds"
    output:
        config['deconv_ref'] + "/{species}/{tissue}_{experiment}_seurat.rds"
    resources: mem_mb=get_mem_mb
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        Rscript {workflow.basedir}/scripts/splitIntoTissues.R {input} {output} {wildcards.tissue}
        """
        
rule reduce_celltype_labels:
    """
    Rule to reduce cell type labels by mapping cell types to their ancestors based on CL ontology.
    """
    conda: "envs/scONTO.yaml"
    log: "logs/reduce_celltype_labels/{species}_{tissue}_{experiment}.log"
    input:
        config['deconv_ref'] + "/{species}/{tissue}_{experiment}_seurat.rds"
    output:
        temp(config['deconv_ref'] + "/{species}/{tissue}_{experiment}_seurat_curated.rds")  
    resources: mem_mb=get_mem_mb
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        Rscript {workflow.basedir}/scripts/reduceCellTypes.R {input}
        """
        
rule generateReferences:
    """
    Rule to generate three different kind of references for deconvolution from the SeuratObject:
     - C0: sc reference [genes x cells]
     - C1: bulk reference [genes x cell types]
     - phenData: metadata for C0
    """
    conda: "envs/refgen.yaml"
    log: "logs/generateReferences/{species}_{tissue}_{experiment}.log"
    input:
        config['deconv_ref'] + "/{species}/{tissue}_{experiment}_seurat_curated.rds"
    output:
        config['deconv_ref'] + "/{species}/{tissue}_{experiment}_C1.rds",
        config['deconv_ref'] + "/{species}/{tissue}_{experiment}_phenData.rds",
        temp(config['deconv_ref'] + "/{species}/{tissue}_{experiment}_C0.rds")
    resources: mem_mb=get_mem_mb
    params:
        method = "none" #do not apply transformation to C matrix
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        Rscript {workflow.basedir}/scripts/genRef_1.R {input} {params.method}
        """

rule scale_C0_reference:
    """
    Rule to scale sc reference (C0) with 'SCTransform'.
    """
    conda: "envs/refgen.yaml"
    log: "logs/scale_C/{species}_{tissue}_{experiment}.log"
    input:
        C_table=config['deconv_ref']  + "/{species}/{tissue}_{experiment}_C0.rds",
        C_phenData=config['deconv_ref'] + "/{species}/{tissue}_{experiment}_phenData.rds"
    output:
        config['deconv_ref'] +  "/{species}/{tissue}_{experiment}_C0_scaled.rds"
    params: scaleCMethod = 'SCTransform'
    threads: 8
    resources: mem_mb=get_mem_mb
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        Rscript {workflow.basedir}/scripts/scaleTable_C.R {input.C_table} {input.C_phenData} {params.scaleCMethod} {threads}
        """

rule UMAP_plots:
    """
    Rule to plot UMAP plots for quality control of references. Created two UMAP plots
    in one png file with old and reduced cell type labels.
    """
    conda: "envs/UMAPplot.yaml"
    log: "logs/UMAP_plots/{species}_{tissue}_{experiment}.log"
    input:
        seurat=config['deconv_ref'] + "/{species}/{tissue}_{experiment}_seurat_curated.rds",
        C0=config['deconv_ref'] + "/{species}/{tissue}_{experiment}_C0_scaled.rds"
    output:
        'UMAP/{species}/{tissue}_{experiment}_umap.png'
    resources: mem_mb=get_mem_mb
    shell: 
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        mkdir -p UMAP/{wildcards.species}
        Rscript {workflow.basedir}/scripts/makeUMAPplots.R {input.seurat} {output}
        """
        
rule reference_summary:
    """
    Rule to create a summary file of the organism part references that were created.
    """
    conda: "envs/scONTO.yaml"
    log: "logs/reference_summary/{species}_summary.log"
    input:
        get_tissues_per_accession
    output:
        config['deconv_ref'] + '/{species}_summary.tsv'
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        Rscript {workflow.basedir}/scripts/makeSummary.R {input} {output}
        """
