#!/usr/bin/env Rscript
# create seuratObjects from mtx files for non ANND scxa experiments

suppressMessages(library(Seurat))
suppressMessages(library(plyr))
suppressMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)

sample <- args[1]
filename_O = paste0('scxa_input/',sample,'/', sample, '.project_seurat.rds')
#read in counts
counts = ReadMtx(mtx = paste0('scxa_input/',sample,'/', sample, '.aggregated_filtered_counts.mtx'), 
             features= paste0('scxa_input/',sample,'/', sample,'.aggregated_filtered_counts.mtx_rows'), 
             cells=   paste0('scxa_input/',sample,'/', sample, '.aggregated_filtered_counts.mtx_cols'))
#read in metadata for cells
metadata = read.csv(paste0('scxa_input/',sample,'/', sample, '.cell_metadata.tsv'), sep = '\t', row.names = 1)

#filter only cells that have metadata
counts = counts[,colnames(counts) %in% rownames(metadata)]

#create SeuratObject
C = CreateSeuratObject(counts = counts, meta.data = metadata)
#rename colum names. If "inferred_cell_type_._ontology_labels_ontology"
# is not present in metadata, try C$inferred_cell_type_._authors_labels_ontology 
#we want to extract the UBERON ids from the links --> basename
first_try <- try(C$cellType
                 <- as.factor(basename(as.character(C$authors_cell_type_._ontology_labels_ontology))))
if("try-error" %in% class(first_try)){
    seconf_try <- try(C$cellType
                     <- as.factor(basename(as.character(C$authors_cell_type_ontology)))
}
if("try-error" %in% class(second_try)){
    stop(paste(sample, 'does not seem to have ontology cell type labels, but this is a requirement to create references'))
}
C$tissue = as.factor(basename(as.character(C$organism_part_ontology)))
C$cellID = colnames(C)

# Remove cells without CL and UBERON IDs
C <- C[, grepl("CL", C$cellType) & grepl("UBERON", C$tissue) ]

#save SeuratObject
saveRDS(C, filename_O)
