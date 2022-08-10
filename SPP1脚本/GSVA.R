library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)
rm(list=ls())
load(file = 'all_macro.RData')
#s <- load("sce.rds")
#DefaultAssay(scRNA) <- "RNA"
#scRNA <- NormalizeData(scRNA)
genesets <- msigdbr(species = "Homo sapiens", category = "H") #C2:KEGG，
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

#按分组取均值
Idents(sce) <- "subcelltype" 
expr <- AverageExpression(sce, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)
head(expr)

# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="gsva") 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)
pdf(paste0("./","hallmarker_GSVA.pdf"),width = 10,height = 9)
pheatmap::pheatmap(gsva.res, show_colnames = T, 
                   scale = "row",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

