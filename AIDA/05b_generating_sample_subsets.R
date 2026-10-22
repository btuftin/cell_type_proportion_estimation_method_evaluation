#sbatch --time=4:00:00 --mem=256GB --wrap="module load r/4.3.1 ; Rscript 05b_generating_sample_subsets.R"

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


AIDA_matrix <- Read10X("datasets_for_comparison/test_set_all/")
AIDA_metadata <- read_csv("datasets_for_comparison/test_set_all/metadata.csv")

since_start()

sc_Ref <- CreateSeuratObject(counts = AIDA_matrix, meta.data = AIDA_metadata)

rm(AIDA_metadata, AIDA_matrix)

rownames(sc_Ref@assays$RNA@layers$counts) <- dimnames(sc_Ref@assays$RNA@features)[[1]]
colnames(sc_Ref@assays$RNA@layers$counts) <- str_replace(dimnames(sc_Ref@assays$RNA@cells)[[1]],
                                                         "-", "_")

# Barcodes are not unique without suffix
length(unique(colnames(sc_Ref@assays$RNA@layers$counts)))

#Cell type names need to contain only letters and _, or ReCIDE will choke
#Some currently have +. This doesn't matter for the pseudobulk,
#but we're also saving the cell type truth and it's convenient if the names
#here match those in the reference matrix

sc_Ref$author_cell_type <- str_replace_all(sc_Ref$author_cell_type, "\\+", "p")

unique(sc_Ref$author_cell_type)

# create pseudo_bulk datasets from subset of remaining subjects
set.seed(7702)
pseudo_samples <- unique(sc_Ref@meta.data[["donor_id"]])

# Double check that these are distinct from reference samples

pseudoAdim_samples <- sample(pseudo_samples, 300)
pseudoB_samples <- sample(pseudoAdim_samples, 200)
pseudoC_samples <- sample(pseudoB_samples, 100)
pseudoD_samples <- sample(pseudoC_samples, 50)
pseudoE_samples <- sample(pseudoD_samples, 10)

# We end up doing this twice for the biggest sample set,
# but that's still some saving in time
sc_Ref = subset(sc_Ref,donor_id %in% pseudoAdim_samples)
#sc_Ref@assays$RNA@counts <- sc_Ref@assays$RNA@layers

for(ref in c("pseudoAdim_samples", "pseudoB_samples", "pseudoC_samples", 
             "pseudoD_samples", "pseudoE_samples")) {
#for(ref in c("pseudoD_samples", "pseudoE_samples")) {
  sc_Ref = subset(sc_Ref,donor_id %in% get(ref))
  pseudo_bulk <- AggregateExpression(sc_Ref, group.by = "donor_id")
  cell_type_truth <- tibble(donor = sc_Ref$donor_id, cell_type = sc_Ref$author_cell_type)
  saveRDS(pseudo_bulk, file = paste0("datasets_for_comparison/", ref, ".rds"))
  gene_names <- unlist(pseudo_bulk$RNA@Dimnames[1])
  pb_df <- as_tibble(pseudo_bulk$RNA) %>% mutate(GeneSymbol = gene_names, .before = 1)
  write_tsv(pb_df, paste0("datasets_for_comparison/", ref, "_for_cibersortX.txt"))
  
  write_csv(cell_type_truth, paste0("datasets_for_comparison/cell_type_truth_", ref, ".csv"))
}
