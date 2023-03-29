#!/usr/bin/env Rscript
## script to reduce the number of cell type labes (if there are to many)
## as this improves deconvolution analysis

# Load required packages
suppressMessages(library(devtools))
suppressMessages(devtools::install_github("YY-SONG0718/scOntoMatch"))
suppressMessages(library(scOntoMatch))
suppressMessages(library(ontologyIndex))
suppressMessages(library(Seurat))
suppressMessages(library(plyr))

# Define a function to reduce the number of cell types in a Seurat object
# this function uses the function scOntoMatch::getOntoMinimal() reduce the 
# granularity of the cell type labels if too many cell types are anotated for
# an organism part

reduce_number_of_cell_types <- function(seurat, cell_ontology, additional_ontologies) {
  onto_id_list <- unique(seurat$cellType)
  max_combine <- 2 # defines how many cell types can be mereged together
  i <- 0
  # we want to run this until we get to the desired number of cell types or
  # max_combine gets to high it would mapp cell types to too general terms
  while (length(unique(seurat$cellType)) > 8 && max_combine < 4) {
    for (ont_id in additional_ontologies) {
      i <- i + 1 
      if (i > length(additional_ontologies)) {
        # increase max_combine if we tried mapping to all the ontologies
        # and see if we find new mappings
        max_combine <- max_combine + 1
        i <- 0
        message(paste0("Maximum combine set to ", max_combine))
      }
      onto_id_list <- unique(seurat$cellType)
      # append onto_id_list with one ont_id from additional_ontologies list
      # we do this step by step to avoid dropping number of cell types to quick
      onto_id_list <- c(onto_id_list, ont_id)
      # mapp onto_ids
      mapping <- getOntoMinimal(cell_ontology, onts = onto_id_list)
      
      value_counts <- table(unlist(mapping))
      # remove CL terms that have more than max_combine or less than 2 CL ids 
      # mapping to them
      names_to_remove <- names(value_counts)[value_counts < 2 | value_counts > max_combine]
      mapping <- Filter(function(x) !x %in% names_to_remove, mapping)
      # if we find at least one CL id we can map to, decrease max_combine to 2 again
      # this is done to try more specific mappings first
      if (length(mapping) > 1) {
        i <- 0
        max_combine <- 2
        message("Maximum combine set back to 2")
      }
      old_names <- names(mapping)
      new_names <- unlist(mapping, use.names = FALSE)
      # Map old cell type names to new names
      seurat$cellType <- mapvalues(seurat$cellType, from = old_names, to = new_names)
      seurat$cellType[is.na(seurat$cellType)] <- 0
      # break loop if get desired number of cell tyes
      if (length(unique(seurat$cellType)) < 9) {
        break
      }
    }
  }
  # get cell type names from CL ids
  cell_type_names <- getOntologyName(unique(seurat$cellType), ont = cell_ontology)
  ids <- names(cell_type_names)
  cell_type_list <- unname(cell_type_names)
  # Map cell type IDs to cell type names
  seurat$cell_type_names <- mapvalues(seurat$cellType, from = ids , to = cell_type_names)
  # return SeuratObject with reduced cell type mappingså
  return(seurat)
}

# Load cell ontology index
obo_file <- "atlas-deconvolution-references/files/cl-basic.obo"
propagate_relationships <- c('is_a')
# create Ontology index
ont <- ontologyIndex::get_OBO(obo_file, propagate_relationships = propagate_relationships)

# get a list of CL ontologies that we can map CL ids predent in SeuratObject to 
additional_ontologies = names(ont$ancestors)
additional_ontologies <- unlist(ont$children, use.names = FALSE)

#sort ontology ids that highest ids are in front as they tend to be more specific
additional_ontologies <- additional_ontologies[order(-as.numeric(sub("CL:", "", additional_ontologies)))]
#remove 'precursor cell', 'progenitor', 'stuff accumulating cell', 'nucleate cell', 'single nucleate cell'
# id to avoid this term for immune cell
ontologies_to_remove = c("CL:0011115", "CL:0011026", "CL:0000325", "CL:0002242", "CL:0000226")
additional_ontologies = additional_ontologies[ !(additional_ontologies %in% ontologies_to_remove)]

args = commandArgs(trailingOnly=TRUE)
filename = args[1]

seurat = readRDS(filename)
seurat$cellType = as.character(seurat$cellType)
# make sure CL ids match CL ids in cl-basic.obo 
seurat$cellType = sub('_', ':', seurat$cellType)

# get cell type labels from list of CL ids
cell_type_names <- getOntologyName(unique(seurat$cellType), ont = ont)
ids <- names(cell_type_names)
cell_type_list <- unname(cell_type_names)
  
# store the old cell type labels for UMAP later
seurat$old_cell_type_names <- mapvalues(seurat$cellType, from = ids , to = cell_type_names)
  
# run function to reduce granularity of cell type labels
seurat = reduce_number_of_cell_types(seurat, ont, additional_ontologies)

#save curated seuratObject
saveRDS(seurat, sub("_seurat.rds", "_seurat_curated.rds", filename))
