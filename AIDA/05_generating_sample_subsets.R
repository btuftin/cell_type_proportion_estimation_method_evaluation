#sbatch --time=4:00:00 --mem=256GB --wrap="module load r/4.3.1 ; Rscript 05_generating_sample_subsets.R"

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

started <- now()
since_start <- function() {
  print(now() - started)
} 

# Matrix data is saved as float, even though it's a count
# double checking that we're not running into issues
# AIDA_matrix <- fread("datasets_for_comparison/matrix.mtx.gz", skip = 3,
#                      nrows = 1000000)

# Ran into issues with fread not able to read more than 2 billion rows
# Remember to check values in subsets instead

# Split original dataset in two, with 200 samples in ref set and
# remaining in test set for pseudobulk generation.

# These should both be small enough to load without hitting limits on
# int

AIDA_matrix <- Read10X("datasets_for_comparison/reference_all")
AIDA_metadata <- read_csv("datasets_for_comparison/reference_all/metadata.csv")

since_start()

sc_Ref <- CreateSeuratObject(counts = AIDA_matrix, meta.data = AIDA_metadata,
                             names.delim = "-")

rm(AIDA_metadata, AIDA_matrix)

rownames(sc_Ref@assays$RNA@layers$counts) <- dimnames(sc_Ref@assays$RNA@features)[[1]]
colnames(sc_Ref@assays$RNA@layers$counts) <- str_replace(dimnames(sc_Ref@assays$RNA@cells)[[1]],
                                                         "-", "_")

# Barcodes are not unique without suffix
length(unique(colnames(sc_Ref@assays$RNA@layers$counts)))

#Cell type names need to contain only letters and _, or ReCIDE will choke
#Some currently have +

sc_Ref$author_cell_type <- str_replace_all(sc_Ref$author_cell_type, "\\+", "p")

unique(sc_Ref$author_cell_type)

# Subset samples
# set seed to enable exact replication of random sampling!
set.seed(7702)
#sample_subjects <- 1:201 # test line
sample_subjects <- unique(sc_Ref@meta.data[["donor_id"]])

print(paste0("Checking number of subjects: ", length(sample_subjects), " sample IDS"))

refA <- sample(sample_subjects, 200)
refB <- sample(refA, 100)
refC <- sample(refB, 50)
refD <- sample(refC, 25)
refE <- sample(refD, 10)

ref <- "refE"



for(ref in c("refA", "refB", "refC", "refD", "refE")) {
#for(ref in c("refD", "refE")) {
    
  sc_Ref = subset(sc_Ref,donor_id %in% get(ref))
  
  rownames(sc_Ref@assays$RNA@layers$counts) <- dimnames(sc_Ref@assays$RNA@features)[[1]]
  colnames(sc_Ref@assays$RNA@layers$counts) <- str_replace(dimnames(sc_Ref@assays$RNA@cells)[[1]],
                                                           "-", "_")
  
  saveRDS(sc_Ref, file = paste0(ref, ".rds"))
#  if(ref %in% c("refD", "refE")) {
    # save cibersort matrix
    ref_rowsums <- rowSums(sc_Ref@assays$RNA@layers$counts)
    ref_rowsums <- ref_rowsums[ref_rowsums > 0]

    sc_Ref <- subset(sc_Ref, features = names(ref_rowsums))
    
    sc_df <- as.matrix(sc_Ref@assays$RNA@layers$counts)
    gene_names <- dimnames(sc_Ref@assays$RNA@features)[[1]]
    sc_df <- as_tibble(sc_df) %>% mutate(GeneSymbol = gene_names, .before = 1)
    colnames(sc_df) <- c("GeneSymbol", as.character(sc_Ref@meta.data$author_cell_type))
    write_tsv(sc_df, paste0("datasets_for_comparison/", ref, "_for_cibersort.txt"))
#  }
}

# Moved pseudobulk generation to 05b because it had an error and
# I wanted cleaner debugging. And also it would have run out of time

