
## 2.重新读取数据#进行Mix_CRC数据巨噬细胞亚组差异分析
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = 'CRC_sce_after_harmony_2.RData')
View(sce@meta.data)
#确定用于细胞分群的PC,放置数据的位置"test1"
dim.use <- 1:50
sam.name <- "Mix_CRC数据巨噬细胞亚组差异分析"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

#### 2.1 计算marker基因 默认2000个高变gene
sce= SetIdent(sce,value="subcelltype") 
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","myloid_subcelltype_marker_genes",max(dim.use),"PC.txt"),sep="\t",quote = F)


#各亚组 全部基因的差异基因C
sce= SetIdent(sce,value="subcelltype") 
all.markers  <- FindAllMarkers(sce, only.pos = TRUE,
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(all.markers ,file=paste0("./",sam.name,"/","NCvsCRC",max(dim.use),"PC.txt"),sep="\t",quote = F)


#全部基因的差异NCvsCRC
markers <- FindMarkers(sce,  ident.1="NC", ident.2="CRC",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","NCvsCRC",max(dim.use),"PC.txt"),sep="\t",quote = F)


#全部基因的差异CRCvsLM
markers <- FindMarkers(sce,  ident.1="CRC", ident.2="LM",
                       assay = 'RNA',slot = 'counts',
                       logfc.threshold =0,min.pct = 0 )
write.table(markers,file=paste0("./",sam.name,"/","CRCvsLM",max(dim.use),"PC.txt"),sep="\t",quote = F)

