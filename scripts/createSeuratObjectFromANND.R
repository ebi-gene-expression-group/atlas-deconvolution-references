## script to read in mtx files from ANND experiments
## and create SeuratObject
suppressMessages(library(devtools))
if (!require("sceasy")) devtools::install_github('cellgeni/sceasy')
suppressMessages(library(sceasy))
suppressMessages(library(Seurat))
suppressMessages(library(plyr))
suppressMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
   stop("Not correct number of arguments. Please supply 2 arguments")
}
accession = args[1]
species = args[2]
#filename of seuratObject
SeuratObject_filename = paste0('scxa_input/', species, '/', accession, '/', accession, '.project_seurat.rds')

#read in counts
counts = ReadMtx(mtx = paste0('scxa_input/', species, '/', accession ,'/matrix.mtx.gz'), 
        features= paste0('scxa_input/', species, '/', accession , '/genes.tsv.gz'), #genes should be ENS ids without .\d
        cells=   paste0('scxa_input/', species, '/', accession ,'/barcodes.tsv.gz'),  feature.column = 1)

#read in metadata for cells
metadata = read.csv(paste0('scxa_input/',  species, '/', accession,'/', accession, '.cell_metadata.tsv'), 
                        sep = '\t', row.names = 1) #row.names should be cellIDs

SeuratObject = CreateSeuratObject(counts = counts, meta.data = metadata)

# rename colum names. If "inferred_cell_type_._ontology_labels_ontology" is 
# not present in metadata, try C$inferred_cell_type_._authors_labels_ontology 
# we want to extract the UBERON ids from the links --> basename
first_try <- try(SeuratObject$cellType
                 <- as.factor(basename(as.character(SeuratObject$authors_cell_type_._ontology_labels_ontology))))
if("try-error" %in% class(first_try)){
    second_try <- try(SeuratObject$cellType
                     <- as.factor(basename(as.character(SeuratObject$authors_cell_type_ontology))))
   if("try-error" %in% class(second_try)){
         stop(paste(accession, 'does not seem to have ontology cell type labels, but this is a requirement to create references'))
   }
 }              
#add new columns
SeuratObject$tissue = as.factor(basename(as.character(SeuratObject$organism_part_ontology)))
SeuratObject$cellID = colnames(SeuratObject)

# Remove cells without CL and UBERON IDs
SeuratObject <- SeuratObject[, grepl("CL", SeuratObject$cellType) & grepl("UBERON", SeuratObject$tissue)]

# change 'exitatory Neuron' id to 'Neuron' as its obsolete
seurat$cellType = mapvalues(seurat$cellType, 'CL_0008030', 'CL_0000540')

#save SeuratObject
saveRDS(SeuratObject, SeuratObject_filename)
