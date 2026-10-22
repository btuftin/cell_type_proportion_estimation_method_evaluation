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

setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/AIDA")

log_this <- function(x){
  cat(paste0(x,"\n"))
}

# refC pseudoB failed with 40GB on the plotting

reference_name <- "refC"
reference_name <- "fewer_cell_types"
bulk_data_name <- "pseudoC"

manual <- TRUE

# skip this line, and the assorted if statements for plots
# when testing the script manually
manual <- FALSE

reference_name <- commandArgs(trailingOnly = TRUE)[1]
reference_file_name <- paste0(reference_name, "_for_cibersort.txt")

bulk_data_name <- commandArgs(trailingOnly = TRUE)[2]
bulk_data_file_name <- paste0(bulk_data_name, "_samples_for_cibersortX.txt")


bk.dat <- fread(file.path("datasets_for_comparison", bulk_data_file_name))
gc()
sample_ids <- colnames(bk.dat)[2:length(colnames(bk.dat))]
genes <- bk.dat$GeneSymbol

bk.dat <- bk.dat[,GeneSymbol := NULL]

bk.dat <- t(as.matrix(bk.dat))
dimnames(bk.dat)[2] <- list(genes)

sc.dat <- fread(file.path("datasets_for_comparison", reference_file_name))
gc()
cell.type.labels <- colnames(sc.dat)[2:length(colnames(sc.dat))]
cell.state.labels <- cell.type.labels
genes <- sc.dat$GeneSymbol

sc.dat <- sc.dat[,GeneSymbol := NULL]

sc.dat <- t(as.matrix(sc.dat))
dimnames(sc.dat)[2] <- list(genes)

gc()

if(manual){
# Examine correlation between cell types

plot.cor.phi (input=sc.dat,
              input.labels=cell.type.labels,
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=0.5, cexCol=0.5,
)
gc() 

sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE, #return the data used for plotting.
  pdf.prefix="refA.gbm.sc.stat" #specify pdf.prefix if need to output to pdf
)
# Mt Ribosomal protein is the only one standing out for refE
# Running refB and sending it to pdf, I also see Ribosomal.
# The plot in RStudio may be hiding data with the legend.


gc()

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE,
  pdf.prefix="gbm.bk.stat" #specify pdf.prefix if need to output to pdf
)
# Similar
}
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  #gene.group=c("Mrp","chrX", "chrY"),
                                  exp.cells=5
)
rm(sc.dat)

#For refB
# Rb      Mrp other_Rb     chrM   MALAT1     chrX     chrY 
# 88       78        5       13        1      846       30 
# A total of  1055  genes from Rb Mrp other_Rb chrM MALAT1 chrX chrY  have been excluded 
# A total of  4022  gene expressed in fewer than  5  cells have been excluded

if (manual) {
plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)
}
# This filtering is suggested to speed up computation
# it seems fine for refE-bulkE, but should be revisited depending
# on QC for larger bulks and for real data
sc.dat.filtered.pc <- select.gene.type (sc.dat.filtered,
                                        gene.type = "protein_coding")
rm(sc.dat.filtered)
# number of genes retained in each category: 
#   
#   protein_coding 
# 14707 

gc()

# Object used by BayesPrism deconvolution algorithm

myPrism <- new.prism(
  reference=sc.dat.filtered.pc,
  mixture=bk.dat,
  input.type="count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

saveRDS(myPrism, file.path("datasets_for_comparison", 
                           paste0(reference_name, "_", bulk_data_name, "_prism.rds")))
