#!/usr/bin/env Rscript
### script to create a summary file of all the organism part references we created

suppressMessages(library(devtools))
if (!require("scOntoMatch")) devtools::install_github("YY-SONG0718/scOntoMatch")
suppressMessages(library(scOntoMatch))
suppressMessages(library(ontologyIndex))
suppressMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)

# read in list of UMAP.pngs we created
sample_list <- args[1:(length(args)-1)]
outfile <- args[length(args)]

if (!file.exists(outfile)) {
write.table(data.frame(
  Accession = character(),
  UBERON_ID = character(),
  organism_part = character(),
  date = character()
  ), 
  outfile, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# read in UBERON ontologies to translate UBERON ids to organism part labels
obo_file = 'deconvolution-reference-library/files/basic.obo' #download here: http://purl.obolibrary.org/obo/uberon/basic.obo
propagate_relationships = c('is_a', 'part_of', 'relationship: part_of', 'synonym')
# create UBRERON ontology object
ont <- ontologyIndex::get_OBO(obo_file,
                              propagate_relationships = propagate_relationships, 
                              extract_tags = 'everything')

# append the summary file for all the references
for (filename in sample_list){

    base = basename(filename)
    uberon_id = str_extract(base, "UBERON_\\d{7}")
    organism_part = unlist(getOntologyName(ont = ont, sub('_', ':', uberon_id)))
    accession = str_extract(base, "E-[A-Z]+-\\d+")
    
    # Append a line with the accession and UBERON ID to the TSV file
    new_line <- data.frame(Accession = accession, 
                           UBERON_ID = uberon_id, 
                           organism_part = organism_part, 
                           date = Sys.Date())
    write.table(new_line, outfile, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)

}
