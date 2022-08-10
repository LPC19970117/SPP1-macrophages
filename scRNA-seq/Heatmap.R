
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
## 2.重新读取数据#进行Mix_CRC数据巨噬细胞亚组差异分析
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = 'Mix_CRC_macro.RData')
View(sce@meta.data)
View(sce@meta.data)

Idents(sce) <- "subcelltype2"


####此处可计算细胞类型的marker用于展示，这里选择用于鉴定细胞类型的marker
#sce$subcelltype2 <- factor(x=sce$subcelltype2,
 #                          levels = c("Macro_C1QC","Macro_FCN1","Macro_SPP1","Macro_MKI67"))

heatmap_gene <- c("SEPP1", "C1QB","C1QA", "C1QC","SLC40A1", "RNASE1","HLA-DQA1", "APOE","FOLR2","LGMN",
                  "S100A8", "FCN1","S100A12", "S100A9","VCAN", "THBS1","EREG","APOBEC3A", "TIMP1","CD300E",
                  "CXCL5", "SPP1","INHBA", "CSTB","CXCL8", "MARCO","CXCL3","SDC2", "C15orf48","FBP1",
                  "MKI67", "STMN1","HIST1H4C", "H2AFZ","TUBA1B", "TUBB","PCLAF","HMGN2","TYMS", "CKS1B")
heatmap_AveE <- AverageExpression(sce, assays = "RNA", features = heatmap_gene,verbose = TRUE) %>% .$RNA

gene_num <- c(10,10,10,10)
gaps_row <- cumsum(gene_num)

cluster_num <- c(1,1,1,1)
gaps_col <- cumsum(cluster_num)
sce= SetIdent(sce,value="subcelltype2") 
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
#annotation_row <- data.frame(row.names = rownames(heatmap_AveE),
#                             `subcelltype2` = rep(factor(subcelltype2,levels = subcelltype2),gene_num))
#annotation_col <- data.frame(row.names = colnames(heatmap_AveE),
#                             `subcelltype2` = rep(factor(subcelltype2,levels = subcelltype2),cluster_num))
#annotation_colors = list(`subcelltype2` = subcelltype2_colors)
#names(annotation_colors$`subcelltype2`) = subcelltype2
pdf(paste0("./","heatmap20220508.pdf"),width = 3,height = 6)
 pheatmap(heatmap_AveE,cluster_cols = F,cluster_rows = F,show_colnames=T,show_rownames=T,
                border=F,#border_color = "white",
                color = c(colorRampPalette(colors = c("#2166ac","#f7fbff"))(length(bk)/2),
                          colorRampPalette(colors = c("#f7fbff","#b2182b"))(length(bk)/2)),
                breaks=bk,scale="row",legend_breaks=seq(-2,2,2),
                gaps_row = gaps_row,gaps_col = gaps_col,
                #annotation_row = annotation_row,annotation_col = annotation_col,
                #annotation_colors = annotation_colors,
                annotation_names_row = F,annotation_names_col = T)
dev.off()

