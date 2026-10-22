setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/AIDA/")

reference_file_name <- "refD.rds"
bulk_data_file_name <- "pseudoD_samples.rds"

reference_name <- commandArgs(trailingOnly = TRUE)[1]
reference_file_name <- paste0(reference_name, ".rds")
print(paste0(reference_file_name, " is the reference file"))
bulk_data_name <- commandArgs(trailingOnly = TRUE)[2]
print(paste0(bulk_data_file_name, " is the bulk data file"))
bulk_data_file_name <- paste0(bulk_data_name, "_samples.rds")
print(paste0(bulk_data_file_name, " is the bulk data file"))

# If loading any of the packages fails, source 00_required_packages.R
source("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/ReCIDE/00_required_packages.R")

#library(DWLS)
source("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/ReCIDE/DWLS/Deconvolution_functions.R")
library(data.table)
library(tidyverse)
library(Seurat)
library(fastSave)
library(COSG)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggbiplot)
library(pbmcapply)
library(mclust)
library(gmodels)

source('/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/ReCIDE/ReCIDE/ReCIDE_all_function.R')

started <- now()
since_start <- function() {
  print(now() - started)
} 

SC_ref <- readRDS(file.path("datasets_for_comparison", reference_file_name))


##Obtaining reference cell type labels and subject labels
##Please note that the cell type labels should not contain any symbols other than "_".

celltype = SC_ref@meta.data[['author_cell_type']] ##In the demo data, it is major.celltype.
subject = SC_ref@meta.data[['donor_id']] ##In the demo data, it is subject.

# load bulk data
bulk_object <- readRDS(file.path("datasets_for_comparison", bulk_data_file_name))

bulk_matrix <- bulk_object$RNA

# filtering out null genes
bulk_rowsums <- rowSums(as.matrix(bulk_matrix))
not_null <- bulk_rowsums[bulk_rowsums > 0]

ref_rowsums <- rowSums(SC_ref@assays$RNA@layers$counts)
ref_rowsums <- ref_rowsums[ref_rowsums > 0]

common_genes <- intersect(names(ref_rowsums),names(not_null))


# 30916 genes in common
# AIDA reference doesn't have the counts slot
# testing if this breaks function, in which case
# I'll have to look at how to add it

# It does break the function, which expects an Assay object,
# not an Assay5 object. Creating a new object
#SC_ref@assays$RNA@counts <- SC_ref@assays$RNA@data

options(Seurat.object.assay.version = "v3")

counts <- SC_ref@assays$RNA@layers$counts
dimnames(counts) <- c(dimnames(SC_ref@assays$RNA@features)[1],
                      dimnames(SC_ref@assays$RNA@cells)[1])

new_SC <- CreateSeuratObject(counts = counts,
                             meta.data = SC_ref@meta.data)

SC_ref <- new_SC

rm(new_SC, counts)
SC_ref <- subset(SC_ref, features = common_genes)


# Filtering bulk data to common genes

bulk_matrix <- bulk_matrix[rownames(bulk_matrix) %in% common_genes,]

common_genes <- rownames(bulk_matrix)

# ReCIDE does an as.numeric on EXP_df[,i], which gives an error for tibbles
# but works for data.frames.
EXP_df <- as.data.frame(bulk_matrix)

# It's not described in the documentation, but EXP_df must have the gene names
# as row names to match with the gene names in the reference
row.names(EXP_df) <- common_genes
since_start()

ref_MGM.list = MGM_build(SC_ref = SC_ref,
                         EXP_df = EXP_df,
                         label_celltype = celltype,
                         label_subject = subject,
                         n_cores = 24)
since_start()

saveRDS(ref_MGM.list, file = paste0("results/temp_", bulk_data_name, "_",
                                    reference_name, "_ref_MGM.list.RDS"))
#ref_MGM.list <- readRDS("ref_MGM.list.RDS")

ReCIDE_results = ReCIDE_deconvolution(Sig_list = ref_MGM.list,
                                      EXP_df = EXP_df,
                                      Method = 'DWLS',
                                      n_cores = 24,
                                      n_celltype = 0)
since_start()
print(paste0(length(common_genes), " genes in common between ", reference_name,
             " and ", bulk_data_name))

final_results <- ReCIDE_results[["results_final_df"]]
cell_types <- data.frame(cell_types = row.names(final_results))
final_results <- cbind(cell_types, final_results)

write.csv(final_results, paste0("results/ReCIDE_result_", reference_name, "_", 
                                bulk_data_name, "_samples.csv"))
