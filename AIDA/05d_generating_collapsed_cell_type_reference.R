#sbatch --time=4:00:00 --mem=128GB --wrap="module load r/4.3.1 ; Rscript 05d_generating_collapsed_cell_type_reference.R"

library(data.table)
library(tidyverse)
library(Seurat)

setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/AIDA/")

#for testing ref_object <- readRDS("datasets_for_comparison/refE.rds")

#ref_object <- readRDS("datasets_for_comparison/refC.rds")
ref_object <- readRDS("datasets_for_comparison/refB.rds")

collapse_types <- function(cell_types) {
  df <- data.frame(ct = cell_types)
  df <- df %>% mutate(ct = case_when(
    ct %in% c("CD14p_Monocyte", "CD16p_Monocyte", "Monocyte")
    ~ "Monocyte",
    ct %in% c("cDC1", "cDC2", "cDC", "pDC", "DC") 
    ~ "DC",
    ct %in% c("naive_B", "IGHMhi_memory_B", "IGHMlo_memory_B",
              "atypical_B", "B")
    ~ "B",
    ct %in% c("dnT", "gdT", "CD4p_T_naive", "CD8p_T_naive", "CD8p_T_GZMKp", 
              "CD8p_T_GZMBp", "CD4p_T_cm", "CD4p_T_em", "CD4p_T", "CD8p_T",
              "MAIT", "Treg", "CD4p_T_cyt", "T")
    ~ "T",
    ct %in% c("CD16p_NK", "CD56p_NK", "NK")
    ~ "NK",
    TRUE ~ ct
  ))
  return(df$ct)
}

ref_object$author_cell_type <- collapse_types(ref_object$author_cell_type)
unique(ref_object$author_cell_type)

saveRDS(ref_object, "datasets_for_comparison/Bfewer_cell_types.rds")

sc_df <- as.matrix(ref_object@assays$RNA@layers$counts)
gene_names <- dimnames(ref_object@assays$RNA@features)[[1]]
colnames(sc_df) <- ref_object@meta.data$author_cell_type 
sc_df <- as_tibble(sc_df) %>% mutate(GeneSymbol = gene_names, .before = 1)
# add names after conversion, so tibble doesn't rename duplicate column names
colnames(sc_df) <- c("GeneSymbol", as.character(ref_object@meta.data$author_cell_type))
write_tsv(sc_df, paste0("datasets_for_comparison/Bfewer_cell_types_for_cibersort.txt"))
