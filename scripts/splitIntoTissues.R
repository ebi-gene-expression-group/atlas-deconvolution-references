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

# Save output
saveRDS(seurat_tissue, args$output)
