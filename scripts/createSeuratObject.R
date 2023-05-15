#!/usr/bin/env Rscript
# create seuratObjects from mtx files for non ANND scxa experiments

suppressMessages(library(Seurat))
suppressMessages(library(plyr))
suppressMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 2) {
   stop("Not correct number of arguments. Please supply 2 arguments")
}
accession = args[1]
species = args[2]

SeuratObject_filename = paste0('scxa_input/', species, '/', accession,'/', accession, '.project_seurat.rds')
#read in counts
counts = ReadMtx(mtx = paste0('scxa_input/',accession,'/', accession, '.aggregated_filtered_counts.mtx'), 
             features= paste0('scxa_input/',accession, '/', accession,'.aggregated_filtered_counts.mtx_rows'), 
             cells=   paste0('scxa_input/',accession,'/', accession, '.aggregated_filtered_counts.mtx_cols'))
#read in metadata for cells
metadata = read.csv(paste0('scxa_input/',accession ,'/', accession, '.cell_metadata.tsv'), sep = '\t', row.names = 1)

#filter only cells that have metadata
counts = counts[,colnames(counts) %in% rownames(metadata)]

#create SeuratObject
SeuratObject = CreateSeuratObject(counts = counts, meta.data = metadata)
#rename colum names. If "inferred_cell_type_._ontology_labels_ontology"
# is not present in metadata, try SeuratObject$inferred_cell_type_._authors_labels_ontology 
#we want to extract the UBERON ids from the links --> basename
first_try <- try(SeuratObject$cellType
                 <- as.factor(basename(as.character(SeuratObject$inferred_cell_type_._ontology_labels_ontology))))
if("try-error" %in% class(first_try)){
    second_try <- try(SeuratObject$cellType
                     <- as.factor(basename(as.character(SeuratObject$authors_cell_type_ontology))))
   if("try-error" %in% class(second_try)){
         stop(paste(accession, 'does not seem to have ontology cell type labels, but this is a requirement to create references'))
   }
 }
SeuratObject$tissue = as.factor(basename(as.character(SeuratObject$organism_part_ontology)))
SeuratObject$cellID = colnames(SeuratObject)

# Remove cells without CL and UBERON IDs
SeuratObject <- SeuratObject[, grepl("CL", SeuratObject$cellType) & grepl("UBERON", SeuratObject$tissue) ]

#save SeuratObject
saveRDS(SeuratObject, SeuratObject_filename)
