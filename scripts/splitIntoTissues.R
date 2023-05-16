#!/usr/bin/env Rscript

# split Seurat Object into seperate organism parts and downsample cell types
# to max 300 cells

library(Seurat)
library(SeuratObject)
library(R.utils)
library(argparse)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("input", help = "Input Seurat object file")
parser$add_argument("output", help = "Output file name")
parser$add_argument("tissue", help = "Tissue name")
args <- parser$parse_args()

# Load Seurat object
if (!file.exists(args$input)) {
  stop("Input file not found.")
}
sc <- readRDS(args$input)

sc$cellType <- as.factor(gsub(":", "_", sc$cellType))
sc$tissue <- as.factor(gsub(":", "_", sc$tissue))

# Subset Seurat object by tissue to get onne seurat object per tissue
if (sum(grepl(args$tissue, sc$tissue)) == 0) {
  stop("Tissue not found in Seurat object.")
}
seurat_tissue <- sc[, sc$tissue == args$tissue ]

# Remove cells without CL and UBERON IDs
seurat_tissue <- seurat_tissue[, grepl("CL", seurat_tissue$cellType) 
                            & grepl("UBERON", seurat_tissue$tissue) ]

# Downsampling to get max 300 cells per celltype
# this is done to reduce the size of the dataset
cell_types <- unique(seurat_tissue$cellType)
seurat_downsampled = list()
 for (ct in cell_types) {
   cells <- colnames(subset(seurat_tissue, cellType == ct))
   if (length(cells) > 300) {
     cells <- sample(cells, 300)
     print(cells)
   }
    append the downsamples cell types
   seurat_downsampled[ct] <- seurat_tissue[,cells]
 }
 seurat_downsampled = merge(
   seurat_downsampled[[1]],
   y = seurat_downsampled[2:length(cell_types)]
 )

# Save output
saveRDS(seurat_downsampled, args$output)
