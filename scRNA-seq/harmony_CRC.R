##==harmony整合多样本==##
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(harmony)
View(sce@meta.data)
sam.name <- "test3"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
#数据输入
#scRNAlist <- readRDS("scRNAlist.rds")
load(file = 'Gao.RData')
Gao_sce<-sce
load(file = 'GSE132465.RData')
GSE132465_sce<-sce
load(file = 'GSE144735.RData')
GSE144735_sce<-sce
load(file = 'GSE146771.RData')
GSE146771_sce<-sce
load(file = 'GSE178318.RData')
GSE178318<-sce
load(file = 'GSE178341.RData')
GSE178341<-sce
rm(sce)
##PCA降维
sce_harmony <- merge(Gao_sce, y=c(GSE132465_sce, GSE144735_sce, GSE146771_sce, GSE178318, GSE178341))
rm(Gao_sce,GSE132465_sce, GSE144735_sce, GSE146771_sce, GSE178318, GSE178341)
View(sce_harmony@meta.data)
sce<-sce_harmony
rm(sce_harmony)
#4.存储数据
save(sce,file=paste0("./",sam.name,"/","Mix_all_myeloid_before_harmony.RData"))
#4.提取CRC亚群
sce_Macrophage<-subset(sce,tissue %in% c("CRC"))
sce<-sce_Macrophage
dim(sce)
rm(sce_Macrophage)
#4.存储数据
save(sce,file=paste0("./",sam.name,"/","Mix_CRC_myeloid_before_harmony.RData"))
#标准化和归一化，PCA
sce <- NormalizeData(sce) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
#sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000) 
#sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
##Harmony整合
system.time({sce <- RunHarmony(sce, group.by.vars = "sample")})
#   用户   系统   流逝 #9万细胞需运行约7100s
# 34.308  0.024 34.324
#降维聚类
ElbowPlot(sce, reduction = "pca",ndims = 50)
#接着就是常规聚类降维，都是基于Harmony的Embeddings矩阵
#1 同2
dim.use<- 1:30
sce <- sce %>% 
  RunUMAP(reduction = "harmony", dims = dim.use) %>% 
  FindNeighbors(reduction = "harmony", dims = dim.use) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
#################################################################

#2
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:20)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.5)
##作图umap
#group_by_cluster #默认为seurat_cluster
pdf(paste0("./",sam.name,"/seurat_cluster",max(dim.use),".pdf"),width = 10,height = 9)
DimPlot(sce, reduction = "umap", label=T) 
dev.off()
####计算myloid.marker基因 ####费时，选做
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","myloid_marker_genes_uamp_",max(dim.use),".txt"),sep="\t",quote = F)
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","CRC_subcelltype2_myloid_marker_genes_uamp_",max(dim.use),".txt"),sep="\t",quote = F)
#再次重命名髓系亚群
cluster5celltype <- c("0"="Macro-C1QC",
                      "1"="Macro-FCN1",
                      "2"="Macro-FCN1",
                      "3"="Macro-SPP1",
                      "4"="cDC2-CD1C",
                      "5"="T-cell doublets",
                      "6"="Macro-SPP1",
                      "7"="Macro-MKI67",
                      "8"="Macro-SPP1",
                      "9"="Macro-CXCL10",
                      "10"="Macro-C1QC",
                      "11"="pDC-LILRA4",
                      "12"="cDC3-LAMP3",
                      "13"="Macro-FCN1",
                      "14"="Macro-ATP",
                      "15"="B-cell doublets",
                      "16"="Macro-FCN1",
                      "17"="cDC1-CLEC9A",
                      "18"="Epithelial-cell doublets"
)
sce[['subcelltype2']] = unname(cluster5celltype[sce@meta.data$seurat_clusters])
View(sce@meta.data) 

pdf(paste0("./",sam.name,"/myeloid_subcelltype2_umap_",max(dim.use),".pdf"),width = 10,height = 9)
DimPlot(sce, reduction = "umap",group.by = 'subcelltype2',label = T)
dev.off()
#group_by_sample
pdf(paste0("./",sam.name,"/","sample",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce , reduction = "umap", group.by='sample') 
dev.off()
pdf(paste0("./",sam.name,"/","tissue",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce , reduction = "umap", group.by='tissue') 
dev.off()
pdf(paste0("./",sam.name,"/","subcelltype",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce , reduction = "umap", group.by='subcelltype') 
dev.off()
pdf(paste0("./",sam.name,"/","sourse",max(dim.use),".pdf"),width = 10,height = 9)
DimPlot(sce , reduction = "umap", group.by='sourse') 
dev.off()
pdf(paste0("./",sam.name,"/","MMR",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce , reduction = "umap", group.by='MMR') 
dev.off()
pdf(paste0("./",sam.name,"/","MSI",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce , reduction = "umap", group.by='MSI') 
dev.off()
pdf(paste0("./",sam.name,"/","sex",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce , reduction = "umap", group.by='sex') 
dev.off()
pdf(paste0("./",sam.name,"/","site",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce , reduction = "umap", group.by='site') 
dev.off()

sce= SetIdent(sce,value="seurat_clusters") 
pdf(paste0("./",sam.name,"/","sourse_split",max(dim.use),".pdf"),width = 25,height = 9)
DimPlot(sce , reduction = "umap", split.by='sourse') 
dev.off()
#单基因umap可视化
pdf(paste0("./",sam.name,"/myloid.dim40.umap",max(dim.use),"PC.pdf"),width = 15,height = 14)
FeaturePlot(sce, features = c("LYZ","LILRA4","CLEC9A","CD1C","LAMP3","CD14","FCGR3A","SPP1","APOE","C1QC","NLRP3","CD163","S100A9","FCN1","MKI67","MARCO"),cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,reduction = "umap")
dev.off()
## 气泡图
features= c("LILRA4","CLEC9A","CLEC10A","CD1C","LAMP3","SPP1","C1QC","FCN1","MKI67",
            "ISG15","CD3D","MZB1","CD68","MRC1","MARCO","CD14","FCGR3A",
            "LYZ","APOE","NLRP3","S100A9","FCGR3B","TPSAB1","EPCAM","CXCL10")
pdf(paste0("./",sam.name,"/ALL_myeloid_DotPlot_re1.5",max(dim.use),".pdf"),width = 8,height = 7)
DotPlot(sce, features = unique(features)) + RotatedAxis()
dev.off()
#聚类树
sce= SetIdent(sce,value="seurat_clusters")               
sce <- BuildClusterTree(
  sce,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/all_myeloid_ClusterTree_",max(dim.use),".pdf"),width = 5,height = 4)
PlotClusterTree(sce)
dev.off()
#聚类树
sce= SetIdent(sce,value="subcelltype2")               
sce <- BuildClusterTree(
  sce,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/all_myeloid_subcelltype2_ClusterTree_",max(dim.use),".pdf"),width = 5,height = 4)
PlotClusterTree(sce)
dev.off()
## 保存harmony后的数据
save(sce,file=paste0("./",sam.name,"/","sce_after_harmony.RData"))

#saveRDS(sce,'sce_harmony.rds')
#sce <- readRDS("sce_harmony.rds")
#细胞构成表
table(sce@meta.data$subcelltype2,sce@meta.data$tissue)
table(sce@meta.data$subcelltype2,sce@meta.data$subcelltype)
table(sce@meta.data$subcelltype2,sce@meta.data$site)
table(sce@meta.data$BRAF)

#4.提取指定单细胞亚群
sce_Macrophage<-subset(sce,tissue %in% c("CRC"))
sce<-sce_Macrophage
dim(sce)
rm(sce_Macrophage)

#4.存储数据myloid
save(sce,file=paste0("./",sam.name,"/","CRC_sce_after_harmony.RData"))





