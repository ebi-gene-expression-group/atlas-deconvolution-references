#!/usr/bin/env Rscript

## Reference matrix type 1 (rows = genenames, cols = celltypes) generation script
##
## @zgr2788

suppressMessages(library(Seurat))
suppressMessages(library(sparseMatrixStats))

#Read data
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
   stop("Not correct number of arguments. Please supply 2 arguments")
}

filename <- args[1]
transform <- args[2]
C_0 <- readRDS(filename)
#filename <- sub("Input/Cell_splits", "Input/References", filename)

#normalize

#Write sc-counts matrix as separate object to save memory later on
C_counts <- C_0@assays$RNA@counts
C_metadata <- C_0@meta.data

#Initialize references
C_1 <- data.frame(row.names = rownames(C_0@assays$RNA@counts))

#Transform C matrix
switch(transform,
  "none" = {
    message('Warning: No transformation applied')
  },
  "log" = {
    C_0@assays$RNA@counts <- log1p(C_0@assays$RNA@counts)
  },
  "sqrt" = {
    C_0@assays$RNA@counts <- sqrt(C_0@assays$RNA@counts)
  },
  "vst" = {  # do not use
    message('Warning: varianceStabilizingTransformation requires high memory and is not recommended.')
    C_0@assays$RNA@counts <- round(C_0@assays$RNA@counts)
    C_0@assays$RNA@counts <- DESeq2::varianceStabilizingTransformation(C_0@assays$RNA@counts)
  }
)


#Split cells by type
cellSplits <- SplitObject(C_0, split.by = "cellType")

#Group all celltype rowsums
for (i in 1:length(unique(C_0@meta.data$cellType)))
{
    C_1[names(cellSplits[i])] <- rowMeans(cellSplits[[i]]@assays$RNA@counts)
}

#Debug
#write.csv(C_1, "C_1.csv")

#Write to RDS
saveRDS(C_1, file = sub("_seurat_curated", "_C1", filename ))
saveRDS(C_counts, file = sub("_seurat_curated", "_C0", filename ))
saveRDS(C_metadata, file = sub("_seurat_curated", "_phenData", filename))

