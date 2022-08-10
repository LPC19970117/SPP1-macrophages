if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()
# If your bioconductor version is previous to 3.9, see the section bellow

## Required
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "pheatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC")

options(stringsAsFactors = FALSE)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(SCENIC)

cell.info <- read.table("../data/Cell.Info.txt", sep = "\t", header = T, row.names = 1)
cell.type <- read.table("../data/CellType.Info.txt", sep = "\t", header = T)
cell.info$CellType <- mapvalues(
  x = cell.info$ClusterID,
  from = cell.type$Cluster,
  to = cell.type$Cell.Type
)
cell.info$CellState <- paste(cell.info$Tissue, cell.info$ClusterID, sep = "_")
head(cell.info)

header <- readLines("../data/Figure2-batch-removed.header.txt")
dge <- fread("../data/Figure2-batch-removed.txt.gz", sep = "\t", header = F, data.table = F)
rownames(dge) <- dge$V1
dge <- dge[, -1]
colnames(dge) <- header
dim(dge)
dge <- as.matrix(dge)

scenicOptions <- initializeScenic(
  org="mgi", # 物种名 or hgnc, or dmel
  dbDir="../cisTarget_databases", # RcisTarget databases location
  dbs="mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather", # file name of motif database
  datasetTitle="SCENIC on Mouse Cell Atlas", # choose a name for your analysis
  nCores=1
)
genesKept <- geneFiltering(
  exprMat=dge, 
  scenicOptions=scenicOptions,
  minCountsPerGene = 1,
  minSamples = 20
)