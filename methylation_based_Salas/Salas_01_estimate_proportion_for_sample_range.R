# Script comes from library creator Lucas A. Salas in email to 
# Laura "license to use software from [link to software]

## Step 0: added by me to install the required sample file

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("FlowSorted.Blood.EPIC")

#installation required additional memory, but completed with 20GB
#install.packages("/proj/raffieldlab/software/FlowSorted.BloodExtended.EPIC_1.1.2.tar.gz")

n_samples_start <- commandArgs(TRUE)[1]
n_samples_end <- commandArgs(TRUE)[2]

print(paste0("Analyzing samples ", n_samples_start, " to ",
             n_samples_end, ". Numbers refer to row in spreadsheet.\n"))

timestamp()

## Step 1: Extract the mix samples
library(FlowSorted.Blood.EPIC)
library(ExperimentHub)
#FlowSorted.Blood.EPIC <- libraryDataGet('FlowSorted.Blood.EPIC')

## Step 1 has a cache directory prompt if directory doesn't exist yet.
## Don't run whole file unless directory exists as every line from the
## previous code line will be taken ass input to the previous prompt

## Step 2 separate the reference from the testing dataset if you want to run
## examples for estimations for this function example

# RGsetTargets <- FlowSorted.Blood.EPIC[,
#                                       FlowSorted.Blood.EPIC$CellType == "MIX"]
# sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,
#                                    seq_len(dim(RGsetTargets)[2]), sep = "_")
# RGsetTargets

# Adding loading and transforming of our data

targets <- read.metharray.sheet("/proj/raffieldlab/projects/JHS/raw_methylation/Hazel_processed_methylation/IDAT_files/")

# Removed because of a very weird bug that only occured for ranges 
# covering some part of 400-500
# Just put the correct end range in and run again

# print(paste0("nrow(targets) result is : ", nrow(targets)))
 if(n_samples_end > nrow(targets)) {
   print("Range end was higher than row number. Adjust and rerun.")
   #   n_samples_end <- nrow(targets)
}
# print(paste0("After checking for OOR. Analyzing samples ", n_samples_start, " to ",
#              n_samples_end, ".\n"))

targets <- targets[n_samples_start:n_samples_end,]

print(paste0(nrow(targets), " sample subset from spreadsheet.\n"))

timestamp()
RGsetTargets <- read.metharray.exp(targets = targets,
                                   verbose = TRUE)
timestamp()
## Step 3: use your favorite package for deconvolution.
## Deconvolute the target data set 450K or EPIC blood DNA methylation.
## We recommend ONLY the IDOL method, the automatic method can lead to severe
## biases.

## We recommend using Noob processMethod = "preprocessNoob" in minfi for the
## target and reference datasets.
## Cell types included are "Bas", "Bmem", "Bnv", "CD4mem", "CD4nv",
## "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", and "Treg"
## Use estimateCellCounts2 from FlowSorted.Blood.EPIC.
## Do not run with limited RAM the normalization step requires a big amount
## of memory resources. Use the parameters as specified below for
## reproducibility.
#
prop_ext <- estimateCellCounts2(RGsetTargets,
                                compositeCellType = "BloodExtended",
                                processMethod = "preprocessNoob",
                                probeSelect = "IDOL",
                                cellTypes = c(
                                  "Bas", "Bmem", "Bnv",
                                  "CD4mem", "CD4nv",
                                  "CD8mem", "CD8nv", "Eos",
                                "Mono", "Neu", "NK", "Treg"),
                                CustomCpGs =if(RGsetTargets@annotation[1]=="IlluminaHumanMethylationEPIC"){
                                  IDOLOptimizedCpGsBloodExtended}else{IDOLOptimizedCpGsBloodExtended450k})
timestamp()
save(prop_ext, file = paste0("/proj/raffieldlab/users/btuftin/Various/methylation_cell_abundance_correlation/results/",
                             "prop_ext_object_samples_",n_samples_start, "-", 
                             n_samples_end, ".RData"))
timestamp()
