library(Seurat)
library(anndata)
library(Matrix)
library(stringr)
library(ggplot2)


data <- read.table("GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv", header=TRUE, row.names= 1, sep="\t")
class(data)

matrix_data <- as.matrix(data)

sparse_data <- as(matrix_data, "sparseMatrix")

seurat_object <- CreateSeuratObject(counts = sparse_data)
pbmc=seurat_object
table(pbmc$orig.ident)



#pbmc[["Sample"]][pbmc[["Sample"]] == "X105A"] <- "X105"
#table(pbmc$Sample)


saveRDS(pbmc, "GSE131928.rds")
