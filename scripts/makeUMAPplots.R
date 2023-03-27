#!/usr/bin/env Rscript
library(Seurat)
library(cowplot)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
filename <- args[1]
output_file_name <- args[2]

# Read in the Seurat object and store it 
seurat_obj <- readRDS(filename)
# run basic steps in order to gerneat UMAP plot from SeuratObject
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
Idents(seurat_obj) = seurat_obj$cell_type_names
#plot old cell type labels before running reduceCellTypes.R
p1 = DimPlot(seurat_obj, 
            reduction = "umap", 
            label = TRUE, 
            group.by = 'old_cell_type_names') + 
            NoLegend() + 
            ggtitle('orig. cell type labels')
# and plot cell type labels after running reduceCellTypes.R
p2 = DimPlot(seurat_obj, 
            reduction = "umap", 
            label = TRUE, 
            group.by = 'cell_type_names' ) + 
            NoLegend() +
            ggtitle('reduced cell type labels')
p = p1 + p2
ggsave(output_file_name, p, width = 24, height = 12, units = "cm")


