## Matrix scaling script for C matrices
##
## @zgr2788

suppressMessages(library(Matrix))
suppressMessages(library(future))


#Read data
args <- commandArgs(trailingOnly = TRUE)
filename_C <- args[1]
phenData_C <- args[2]
method <- args[3]
plan('multisession', workers = args[4]) #Paralellism


C <- readRDS(filename_C)
C_data <- readRDS(phenData_C)
filename_C<- sub("sc_References/Raw", "sc_References/Normalized_tables/", filename_C)


message(paste0("Scaling C with method: ", method))

#Preprocess matrices to avoid errors

#Only get rows where there is at least one cell with a count
C <- C[rowSums(C) != 0,]

#Same for cols
C <- C[,colSums(C) != 0]

#Only keep rows with different counts after transformations
C <- C[!apply(C, 1, function(x) var(x) == 0),]


#Switch modes
switch(method,


  "none" = {
    message('Warning: C not scaled')
  },


  "column" = {
    C <- apply(C, 2, function(x) x/sum(x))
  },


  "row" = {
    C <- t(apply(C, 1, function(x) x/sum(x)))
  },


  "mean" = {
    C <- apply(C, 2, function(x) x - mean(x))
  },


  "column_z-score" = {
    C <- scale(C, center = TRUE, scale = TRUE)
  },


  "global_z-score" = {
    C <- (as.matrix(C) - mean(as.matrix(C))) / sd(as.matrix(C))
  },


  "column_min-max" = {
    C <- apply(C, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  },


  "global_min-max" = {
    C <- (C - min(C))/(max(C) - min(C))
  },


  "LogNormalize" = {
    suppressMessages(library(Seurat))

    C <- as.matrix(expm1(Seurat::LogNormalize(C, verbose = FALSE)))
  },


  "QN" = {
    suppressMessages(library(preprocessCore))

    C_rownames <- rownames(C); C_colnames <- colnames(C)
    C <- preprocessCore::normalize.quantiles(as.matrix(C))
    rownames(C) <- C_rownames; colnames(C) <- C_colnames
  },


  "TMM" = { #Spikes? Data?
    suppressMessages(library(edgeR))

    cellType <- as.character(C_data$cellType[rownames(C_data) %in% colnames(C)])
    C <- edgeR::DGEList(counts = C, group = cellType)
    C <- edgeR::calcNormFactors(C, method = "TMM")
    C <- edgeR::cpm(C)
  },


 #sc-RNA specific
  "SCTransform" = {
    suppressMessages(library(sctransform))

    C <- sctransform::vst(C, return_corrected_umi=TRUE, verbosity = FALSE)$umi_corrected
  },


  "scran" =  {
    suppressMessages(library(scran))
    suppressMessages(library(SingleCellExperiment))


    C <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(C)))
    C <- scran::computeSumFactors(C, clusters = NULL)
    C <- scater::logNormCounts(C, exprs_values = 'counts', transform = "none")
    C <- C@assays@data$normcounts
  },


  "scater" = {
    suppressMessages(library(scater))

    sizeFactors <- scater::librarySizeFactors(C)
    C <- scater::normalizeCounts(C, sizeFactors, transform = "none")
  },


  "Linnorm" = {
    suppressMessages(library(Linnorm))

    C <- expm1(Linnorm::Linnorm(as.matrix(C)))
  }


)

#Only get rows where there is at least one cell with a count
C <- C[rowSums(C) != 0,]

#Same for cols
C <- C[,colSums(C) != 0]

#Only keep rows with different counts after transformations
C <- C[!apply(C, 1, function(x) var(x) == 0),]


C <- as(C, "matrix")


#Save transformed tables
saveRDS(C, sub(".rds", "_scaled.rds",filename_C))
