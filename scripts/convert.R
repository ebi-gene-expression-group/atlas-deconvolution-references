## Format conversion script for anndata -> Seurat object conversions using sceasy.
##
## @zgr2788


suppressMessages(library(devtools))
devtools::install_github('cellgeni/sceasy')
suppressMessages(library(sceasy))
suppressMessages(library(Seurat))
suppressMessages(library(plyr))
suppressMessages(library(stringr))


args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]

# filename = paste0('scxa_input/',sample,'/', sample, '.project.h5ad')
# #convert h5ad to seuratBbject
# sceasy::convertFormat(filename, from = "anndata", to = "seurat", 
#         outFile = paste0(sub(".h5ad", "", filename), "_seurat.rds"))

#filename of seuratObject
filename_O = paste0('scxa_input/',sample,'/', sample,  '.project_seurat.rds')
#read in counts
#seuratObject = readRDS(filename_O)
counts = ReadMtx(mtx = paste0('scxa_input/',sample,'/matrix.mtx.gz'), 
        features= paste0('scxa_input/',sample, '/genes.tsv.gz'), 
        cells=   paste0('scxa_input/',sample,'/barcodes.tsv.gz'),  feature.column = 1)

#counts = GetAssayData(object = seuratObject, slot = "counts")
#read in metadata for cells
metadata = read.csv(paste0('scxa_input/',sample,'/', sample, '.cell_metadata.tsv'), 
                        sep = '\t', row.names = 1) #row.names should be cellIDs

#filter only cells that have metadata
#counts = counts[,colnames(counts) %in% rownames(metadata)]

#create SeuratObject
C = CreateSeuratObject(counts = counts, meta.data = metadata)

# C = AddMetaData(object = seuratObject, metadata = metadata)
print(C)
#rename colum names. If "inferred_cell_type_._ontology_labels_ontology" is 
# not present in metadata, try C$inferred_cell_type_._authors_labels_ontology 
#we want to extract the UBERON ids from the links --> basename
t <- try(C$cellType 
        <- as.factor(basename(as.character(C$authors_cell_type_._ontology_labels_ontology))))
if("try-error" %in% class(t)){
    C$cellType <- as.factor(basename(as.character(C$authors_cell_type_ontology)))
}
C$tissue = as.factor(basename(as.character(C$organism_part_ontology)))
C$cellID = colnames(C)

# Remove cells without CL and UBERON IDs
C <- C[, grepl("CL", C$cellType) & grepl("UBERON", C$tissue)]


#save SeuratObject
saveRDS(C, filename_O)
