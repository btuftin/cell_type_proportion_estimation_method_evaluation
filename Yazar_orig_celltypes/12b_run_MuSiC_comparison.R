#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::available()
#BiocManager::install(c("TOAST"), update = FALSE)
#devtools::install_github('xuranw/MuSiC')
#install.packages("remotes")
#remotes::install_github("renozao/xbioc")
library(tidyverse)
library(data.table)
library(MuSiC)
library(xbioc)
library(Biobase)
library(SingleCellExperiment)

log_this <- function(x){
  cat(paste0(x,"\n"))
}

setwd( "/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/")

reference_file_name <- "refC.rds"
bulk_data_file_name <- "pseudoC_samples.rds"

reference_name <- "refC"
bulk_data_name <- "pseudoC"

#gene_subset <- "min"

reference_name <- commandArgs(trailingOnly = TRUE)[1]
reference_file_name <- paste0(reference_name, ".rds")
print(paste0(reference_file_name, " is the reference file"))
bulk_data_name <- commandArgs(trailingOnly = TRUE)[2]
print(paste0(bulk_data_file_name, " is the bulk data file"))
bulk_data_file_name <- paste0(bulk_data_name, "_samples.rds")
print(paste0(bulk_data_file_name, " is the bulk data file"))

#gene_subset <- commandArgs(trailingOnly = TRUE)[3]


bulkRNA <- readRDS(file.path("datasets_for_comparison", bulk_data_file_name))

bulkRNA.genes <- bulkRNA$RNA@Dimnames[[1]]

scRNA <- readRDS(file.path("datasets_for_comparison", reference_file_name))

# Getting the counts, this frame also has row and column names with
# genes and bar codes
scRNA.data <- scRNA@assays$RNA@data

# Getting cell types and deleting original object, as it is very large
# and I'm having trouble with RStudio requests with 40Gb memory allocation
scRNA.celltypes <- scRNA$cell_type


donor_id <- scRNA@meta.data$donor_id

rm(scRNA)
gc()


# Seurat
scRNA.genes <- rownames(scRNA.data)
scRNA.genes <- unique(scRNA.genes)

bulk_matrix <- bulkRNA$RNA

# filtering out null genes
bulk_rowsums <- rowSums(as.matrix(bulk_matrix))
not_null <- bulk_rowsums[bulk_rowsums > 0]

ref_rowsums <- rowSums(as.matrix(scRNA.data))
ref_rowsums <- ref_rowsums[ref_rowsums > 0]

gc()

common_ENSGID <- intersect(names(ref_rowsums),names(not_null))

log_this(paste(length(common_ENSGID), "genes in common"))

# 20k or so genes in common makes MuSiC run very, very slow
# The speed issues were mainly due to an error in the code
# when converting the Seurat data to a SingleCellExperiment,
# so dropping this section.
# testing two different subsets
# if (sum(gene_subset %in% c("min", "max"))) {
# 	common_ENSGID <- read_csv(paste0("results/selected_genes_for_music_", gene_subset, ".csv"))
# 	common_ENSGID <- common_ENSGID %>% pull(genes)
# } else {
# 	gene_subset <- "all"
# }
# 
# log_this(paste("using ", length(common_ENSGID), "genes from gene subset file"))

log_this("Filtering scRNA")
#seurat
scRNA.data <- scRNA.data[which(rownames(scRNA.data) %in% common_ENSGID),]
log_this("scRNA basic filter complete")
log_this(paste(nrow(scRNA.data), "rows remaining"))

scRNA.genes <- rownames(scRNA.data)
scRNA.genes <- unique(scRNA.genes)

bulkRNA <- bulkRNA$RNA[rownames(bulkRNA$RNA) %in% common_ENSGID,]
bulkRNA.genes <- rownames(bulkRNA)


# Creating bulkRNA data frame
pheno.data <- colnames(bulkRNA)
pheno.data.df <- data.frame(sampleID=pheno.data)
rownames(pheno.data.df) <- colnames(bulkRNA)
pd <- new('AnnotatedDataFrame', data = pheno.data.df)

# Creating Bulk data ExpressionSet
bulk.est <- ExpressionSet(assayData = as.matrix(bulkRNA), phenoData=pd)

log_this("Creating scRNA singleCellExperiment")
pheno.scRNA.df <- data.frame(sampleID=donor_id, 
                             cellType=scRNA.celltypes)
rownames(pheno.scRNA.df) <- colnames(scRNA.data)

gene.scRNA.df <- data.frame(gene.name = rownames(scRNA.data))
rownames(gene.scRNA.df) <- rownames(scRNA.data)

#pd <- new('AnnotatedDataFrame', data = pheno.scRNA.df)
#scRNA.est <- ExpressionSet(assayData = as.matrix(scRNA.data), phenoData=pd)

scRNA.sce <- SingleCellExperiment(list(counts=scRNA.data), colData = pheno.scRNA.df, 
                                  rowData = gene.scRNA.df)

bulk.mtx <- exprs(bulk.est)
# 
log_this("Running music_prop")
### Cell types proportion estimates:
Est.prop <- music_prop(bulk.mtx = bulk.mtx, sc.sce = scRNA.sce, clusters='cellType',
                       markers = rownames(bulk.mtx),
                      samples = 'sampleID', 
                      select.ct = levels(scRNA.celltypes))

log_this("saving results")
save(Est.prop, file = file.path("results", paste0("MuSiC_results_", reference_name, "_", 
bulk_data_name,  "_pseudobulk_samples.RData")))
cell_types <- colnames(Est.prop$Est.prop.weighted)
est_prop_weighted <- as_tibble(t(Est.prop$Est.prop.weighted))
est_prop_weighted <- cbind(cell_types, est_prop_weighted)
write_tsv(est_prop_weighted, file = file.path("results", paste0("MuSiC_results_", reference_name, "_", 
  bulk_data_name,  "_pseudobulk_samples_genes.tsv")))

#summary(Est.prop$Est.prop.weighted)

#load(file = file.path(current_wd, "JHS_RNAseq_MuSiC.RData"))

#log_this("Creating boxplot")
#png(file.path(current_wd, "prop_box_JHS_RNAseq_seurat_reference.png"), width = 1280, height = 256)
#boxplot(Est.prop$Est.prop.weighted, col = rainbow(ncol(Est.prop$Est.prop.weighted)),ylab="Proportion", main = "cell-type proportion from MuSiC")
#dev.off()
