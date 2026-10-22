library(data.table)
#install.packages("snowfall")
library(snowfall)
#install.packages("NMF")
library(NMF)
library(gplots)
#install.packages("scran")
#BiocManager::install("scran")
library(scran)
library(BiocParallel)
library(Matrix)
#library(devtools)
# IF install_github fails due to API rate limit run the following commands
#usethis::create_github_token()

# Paste token here

#usethis::edit_r_environ()
#

library(BayesPrism)

setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/")

log_this <- function(x){
  cat(paste0(x,"\n"))
}
reference_name <- "refE"
bulk_data_name <- "pseudoE"

reference_name <- commandArgs(trailingOnly = TRUE)[1]
reference_file_name <- paste0(reference_name, "_for_cibersortX.txt")

bulk_data_name <- commandArgs(trailingOnly = TRUE)[2]
bulk_data_file_name <- paste0(bulk_data_name, "_samples_for_cibersortX.txt")

cores <- 1
cores <- as.numeric(commandArgs(trailingOnly = TRUE)[3])

myPrism <- readRDS(file.path("datasets_for_comparison", 
                             paste0(reference_name, "_", bulk_data_name, "_prism.rds")))

bp.res <- run.prism(prism = myPrism, n.cores=cores)

bp.res

theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")

theta_df <- as.data.frame(theta)

theta_df$sample_id <- rownames(theta_df)

write.csv(theta, 
          file.path("results",
          paste0("BayesPrism_", reference_name, "_", bulk_data_name, "_results.csv")))
