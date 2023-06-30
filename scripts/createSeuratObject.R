#!/usr/bin/env Rscript
# create seuratObjects from mtx files for non ANND scxa experiments

suppressMessages(library(Seurat))
suppressMessages(library(plyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 6) {
   stop("Not correct number of arguments. Please supply 6 arguments")
}
accession = args[1]
species = args[2]
mtx_path = args[3]
cols_path = args[4]
rows_path = args[5]
sdrf_file = args[6]

get_cell_metadata <- function(sdrf_file) {
  # Read the condensed SDRF file
  data <- read.csv(sdrf_file, sep = "\t", header = FALSE)
  
  # Pivot the dataframe
  pivoted_df <- data.frame(pivot_wider(data, names_from = c(V4, V5), values_from = V7, id_cols = V3))
  
  # Set the row names
  row.names(pivoted_df) <- pivoted_df$V3
  
  # Remove the unnecessary column
  pivoted_df <- pivoted_df[, -1]
  
  # Function to check if a row contains URLs matching specific patterns
  contains_matching_urls <- function(row) {
    # Specify the regular expression pattern
    pattern <- "(http://purl.obolibrary.org/obo/CL_\\d+)|(http://purl.obolibrary.org/obo/UBERON_\\d+)"
    
    # Check if any element in the row matches the pattern
    any(grepl(pattern, row))
  }
  
  # Get the row indices that contain matching URLs
  matching_rows <- which(apply(pivoted_df, 1, contains_matching_urls))
  
  # Return the cell metadata dataframe
  return(pivoted_df[matching_rows, ])
}

# get cell metadata file from condensed sdrf file
cell_metadata <- get_cell_metadata(sdrf_file)

SeuratObject_filename = paste0('scxa_input/', species, '/', accession,'/', accession, '.project_seurat.rds')
#read in counts
counts = ReadMtx(mtx = mtx_path, 
                 features = rows_path, 
                 cells = cols_path)

#filter only cells that have metadata
counts = counts[,colnames(counts) %in% rownames(cell_metadata)]

#create SeuratObject
SeuratObject = CreateSeuratObject(counts = counts, meta.data = cell_metadata)
#rename colum names. If "inferred_cell_type_._ontology_labels_ontology"
# is not present in metadata, try SeuratObject$inferred_cell_type_._authors_labels_ontology 
#we want to extract the UBERON ids from the links --> basename
first_try <- try(SeuratObject$cellType
                 <- as.factor(basename(as.character(SeuratObject$factor_inferred.cell.type...ontology.labels))))
if("try-error" %in% class(first_try)){
    second_try <- try(SeuratObject$cellType
                     <- as.factor(basename(as.character(SeuratObject$characteristic_cell.type))))
   if("try-error" %in% class(second_try)){
         stop(paste(accession, 'does not seem to have ontology cell type labels, but this is a requirement to create references'))
   }
 }
SeuratObject$tissue = as.factor(basename(as.character(SeuratObject$characteristic_organism.part)))
SeuratObject$cellID = colnames(SeuratObject)

# Remove cells without CL and UBERON IDs
SeuratObject <- SeuratObject[, grepl("CL", SeuratObject$cellType) & grepl("UBERON", SeuratObject$tissue) ]

# change 'exitatory Neuron' id to 'Neuron' as its obsolete
SeuratObject$cellType = mapvalues(SeuratObject$cellType, 'CL_0008030', 'CL_0000540')

#save SeuratObject
saveRDS(SeuratObject, SeuratObject_filename)
