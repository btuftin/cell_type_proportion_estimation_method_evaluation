library(Biobase)
library(BisqueRNA)
library(data.table)
library(tidyverse)
library(Seurat)

setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/")

log_this <- function(x){
  cat(paste0(x,"\n"))
}

reference_file_name <- "refC.rds"
bulk_data_file_name <- "pseudoC_samples.rds"

reference_name <- commandArgs(trailingOnly = TRUE)[1]
reference_file_name <- paste0(reference_name, ".rds")

bulk_data_name <- commandArgs(trailingOnly = TRUE)[2]
bulk_data_file_name <- paste0(bulk_data_name, "_samples.rds")


bulkRNA <- readRDS(file.path("datasets_for_comparison", bulk_data_file_name))

# This is 17GB all on its own, so size your virtual machine accordingly
# Running with a substantial fraction of the genes available requires
# tons more memory anyway
scRNA <- readRDS(file.path("datasets_for_comparison", reference_file_name))

# We're extracting the data from the Seurat object later, and need
# to save the donor_ids first
scRNA.individual_ids <- scRNA$donor_id

scRNA.celltypes <- scRNA$cell_type


# Getting the counts, this frame also has row and column names with
# genes and bar codes
scRNA.data <- scRNA@assays$RNA@data

gc(verbose=FALSE)
rm(scRNA)  

bulk_matrix <- bulkRNA$RNA

# filtering out null genes
bulk_rowsums <- rowSums(as.matrix(bulk_matrix))
not_null <- bulk_rowsums[bulk_rowsums > 0]

ref_rowsums <- rowSums(as.matrix(scRNA.data))
ref_rowsums <- ref_rowsums[ref_rowsums > 0]

gc()

common_ENSGID <- intersect(names(ref_rowsums),names(not_null))

log_this(paste(length(common_ENSGID), "genes in common"))

log_this("Filtering scRNA")
#seurat
scRNA.data <- scRNA.data[which(rownames(scRNA.data) %in% common_ENSGID),]
log_this("scRNA basic filter complete")
log_this(paste(nrow(scRNA.data), "rows remaining"))

scRNA.genes <- rownames(scRNA.data)
scRNA.genes <- unique(scRNA.genes)

bulkRNA <- bulkRNA$RNA[rownames(bulkRNA$RNA) %in% common_ENSGID,]
bulkRNA.genes <- rownames(bulkRNA)

  

  log_this("scRNA basic filter complete")
  log_this(paste(nrow(scRNA.data), "rows remaining"))

  bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulkRNA))
  
  scRNA.sample_ids <- colnames(scRNA.data)
  
  sc.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors = F,
                         row.names = scRNA.sample_ids,
                         SubjectName = scRNA.individual_ids,
                         cellType=scRNA.celltypes)
  sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"),
                        row.names=c("SubjectName","cellType"))
  sc.pdata <- new("AnnotatedDataFrame",data=sc.pheno, varMetadata=sc.meta)
  
  sc.eset<-Biobase::ExpressionSet(assayData = as.matrix(scRNA.data),phenoData=sc.pdata)
  
  log_this("analyzing")
  res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers= NULL,
                                                use.overlap=FALSE)  
  
  log_this("saving results")
  #  save(Est.prop, file = file.path(current_wd, paste0(test_prefix, "JHS_RNAseq_MuSiC_new_reference.RData")))
  write.table(res$bulk.props,
              file = paste0("results/Bisque_results_", reference_name, "_", 
                            bulk_data_name, "_pseudobulk_samples.tsv"),
              quote = F,sep="\t",row.names = T,col.names = T)
  
