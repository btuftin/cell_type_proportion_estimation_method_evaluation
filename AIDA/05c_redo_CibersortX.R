#sbatch --time=4:00:00 --mem=128GB --wrap="module load r/4.3.1 ; Rscript 05c_redo_CibersortX.R"

# In the update to AIDA I lost the renaming of the columns in the reference
# which caused a fatal error, and of course would have given nonsensical results

library(data.table)
library(tidyverse)
library(Matrix);
library(preprocessCore)
library(MASS)
library(magrittr)
library(Biobase)
library(Seurat)
library(SingleCellExperiment)

setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/AIDA/")


#for(ref in rev(c("refA", "refB", "refC", "refD", "refE"))) {
for(ref in c("refA")) {
    gc()
  sc_Ref = readRDS(paste0("datasets_for_comparison/", ref, ".rds"))

    sc_df <- as.matrix(sc_Ref@assays$RNA@layers$counts)
    gene_names <- dimnames(sc_Ref@assays$RNA@features)[[1]]
    colnames(sc_df) <- sc_Ref@meta.data$author_cell_type 
    sc_df <- as_tibble(sc_df) %>% mutate(GeneSymbol = gene_names, .before = 1)
    # add names after conversion, so tibble doesn't rename duplicate column names
    colnames(sc_df) <- c("GeneSymbol", as.character(sc_Ref@meta.data$author_cell_type))
    write_tsv(sc_df, paste0("datasets_for_comparison/", ref, "_for_cibersort.txt"))
    rm(sc_df, sc_Ref)
#  }
}
