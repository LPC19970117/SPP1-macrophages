library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
rm(list = ls())
# 读取多个csv
data_LN <- read.csv("GSE164522_CRLM_LN_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_LN <- CreateSeuratObject(
  data_LN,
  project = "GSE164522", 
  min.cells = 10,
  min.features = 200)
save(sce_LN,file=paste0("./","sce_LN.RData"))
rm(data_LN)

data_MN <- read.csv("GSE164522_CRLM_MN_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_MN <- CreateSeuratObject(
  data_MN,
  project = "GSE164522", 
  min.cells = 10,
  min.features = 200)
save(sce_MN,file=paste0("./","sce_MN.RData"))
rm(data_MN)

data_MT <- read.csv("GSE164522_CRLM_MT_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_MT <- CreateSeuratObject(
  data_MT,
  project = "GSE164522", 
  min.cells = 10,
  min.features = 200)
save(sce_MT,file=paste0("./","sce_MT.RData"))
rm(data_MT)

data_PBMC <- read.csv("GSE164522_CRLM_PBMC_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_PBMC <- CreateSeuratObject(
  data_PBMC,
  project = "GSE164522", 
  min.cells = 10,
  min.features = 200)
save(sce_PBMC,file=paste0("./","sce_PBMC.RData"))
rm(data_PBMC)

data_PN<- read.csv("GSE164522_CRLM_PN_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_PN <- CreateSeuratObject(
  data_PN,
  project = "GSE164522", 
  min.cells = 10,
  min.features = 200)
save(sce_PN,file=paste0("./","sce_PN.RData"))
rm(data_PN)

data_PT <- read.csv("GSE164522_CRLM_PT_expression.csv.gz", header = T,stringsAsFactors = F,row.names= 1)
sce_PT <- CreateSeuratObject(
  data_PT,
  project = "GSE164522", 
  min.cells = 10,
  min.features = 200)
save(sce_PT,file=paste0("./","sce_PT.RData"))
rm(data_PT)

# 合并seurat除LN外的5种组织
sce<- merge(sce_MN, y=c(sce_MT, sce_PBMC, sce_PN, sce_PT))
save(sce,file=paste0("./","sce0_of_5_tissue.RData"))
rm(sce_MN,sce_MT, sce_PBMC, sce_PN, sce_PT)
dim(sce)

#这里文件夹的名字可以修改，但最好只用英文字母
sam.name <- "test1"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
### QC
#############################################################################################################
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
#画图的代码需要选中从pdf一直到dev.off()的所有代码，一起运行
pdf(paste0("./",sam.name,"/QC-VlnPlot.pdf"),width = 8,height = 4.5)
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
dev.off()

plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pdf(paste0("./",sam.name,"/QC-FeatureScatter.pdf"),width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)
#### 5. 筛选细胞 ####
cat("Before filter :",nrow(sce@meta.data),"cells\n")
sce <- subset(sce, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
cat("After filter :",nrow(sce@meta.data),"cells\n")

# ## Normalization
# ############################################################################################################
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)  # Checkpoint
# ## Feature selection
# ############################################################################################################
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sce), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = paste0(sam.name,"/Norm-feature_variable_plot.pdf"),width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)
#均一化（需要一点时间）
#这里做均一化的目的是用于消除特定变量对于整体数据的影响，如批次效应，线粒体数量的差异等，可以同时给出多个因素
#主要变更vars.to.regress即可，具体因素需要在sce@meta.data中存在
#这是所有基因的归一化，内存不够就选2000高变基因的归一化
#all.genes <- rownames(sce) 
#sce <- ScaleData(sce, features = all.genes, vars.to.regress = c("percent.mt"))
#这是2000高变基因的归一化
sce <- ScaleData(sce,
                 vars.to.regress = c("percent.mt"))
#############################################################################################################
# ## Reduction
#############################################################################################################
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
DimPlot(sce, reduction = "pca")
dev.off()
pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 5,height = 4)
ElbowPlot(sce, reduction = "pca",ndims = 50)
dev.off()


### Cluster
dim.use <- 1:50
# ############################################################################################################
sce <- RunUMAP(sce, dims = dim.use) # umap tsne
sce <- FindNeighbors(sce, dims = 1:50)  # louvain cluster, graph based
sce <- FindClusters(sce, resolution = 0.8)

pdf(paste0("./",sam.name,"/CellCluster-DimPlot_umap_",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce, reduction = "umap",group.by = 'seurat_clusters',label = T)
dev.off()

pdf(paste0("./",sam.name,"/umap.",max(dim.use),"PC.pdf"),width = 15,height = 14)
FeaturePlot(sce, features = c("CD3D","KLRF1","MS4A1","MZB1","LYZ","TPSAB1","EPCAM","DCN","TNFRSF17","COL1A1","COL1A2","PLVAP","LILRA4","FCGR3B"),
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()
pdf(paste0("./",sam.name,"/umap2.",max(dim.use),"PC.pdf"),width = 15,height = 14)
FeaturePlot(sce, features = c("CD3D", "KLRF1", "MS4A1",'MZB1', "LYZ","TPSAB1", "FCGR3B", "CD14", "CD1C"),
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()
#GSE164522的大类分型
pdf(paste0("./",sam.name,"/GSE164522的大类分型.",max(dim.use),"PC.pdf"),width = 15,height = 14)
FeaturePlot(sce, features = c("CD3D", "CD3E", "NKG7",'CD79A', "CD19","CST3", "CD68", "LYZ", "CAP3"),
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()

#SPP1
pdf(paste0("./",sam.name,"/SPP1_o.",max(dim.use),"PC.pdf"),width = 4,height = 4)
FeaturePlot(sce, features = c("SPP1"),
            cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,
            reduction = "umap")
dev.off()

## 气泡图
features= c("CD3D","KLRF1","MS4A1","MZB1","LYZ","TPSAB1","EPCAM","DCN","FCGR3B",
            "LILRA4","CD68","COL1A2","FCGR3B","PLVAP","CD14")

pdf(paste0("./",sam.name,"/DotPlot_umap.",max(dim.use),"PC.pdf"),width = 10,height = 10)
DotPlot(sce, features = unique(features)) + RotatedAxis()
dev.off()
#细胞类群相似树
sce <- BuildClusterTree(
  sce,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/CellCluster-ClusterTree_",max(dim.use),"PC.pdf"),width = 10,height = 8)
PlotClusterTree(sce)
dev.off()
#以下重命名能在meta.data中显示
#方法二：使用unname函数配合向量：
cluster2celltype <- c( "0"="T",
                       "1"="T", 
                       "2"="T", 
                       "3"= "NK", 
                       "4"= "T", 
                       "5"= "Plasma",
                       "6"= "B", 
                       "7"= "T", 
                       "8"= "NK",
                       "9"= "Myeloid",
                       "10"="Myeloid",
                       "11"="NK", 
                       "12"="T", 
                       "13"= "B", 
                       "14"= "T", 
                       "15"= "T",
                       "16"= "Plasma", 
                       "17"= "T", 
                       "18"= "T",
                       "19"= "Myeloid",
                       "20"="T",
                       "21"="T", 
                       "22"="Plasma", 
                       "23"= "B",
                       "24"= "T", 
                       "25"= "T",
                       "26"= "T", 
                       "27"= "T", 
                       "28"= "Epithelial",
                       "29"= "NK",
                       "30"="Myeloid",
                       "31"="B", 
                       "32"="Myeloid", 
                       "33"= "NK", 
                       "34"= "Plasma", 
                       "35"= "Mast",
                       "36"= "Myeloid", 
                       "37"= "T", 
                       "38"= "T",
                       "39"= "Myeloid",
                       "40"="T",
                       "41"="Stromal", 
                       "42"="T", 
                       "43"= "Myeloid",
                       "44"= "Stromal", 
                       "45"= "B",
                       "46"= "T", 
                       "47"= "Stromal"
)
sce[['celltype']] = unname(cluster2celltype[sce@meta.data$seurat_clusters])
View(sce@meta.data) 

pdf(paste0("./",sam.name,"/celltype-DimPlot_umap_",max(dim.use),".pdf"),width = 5,height = 4)
DimPlot(sce, reduction = "umap",group.by = 'celltype',label = T)
dev.off()
## 气泡图
sce= SetIdent(sce,value="celltype") 
pdf(paste0("./",sam.name,"/celltype_DotPlot_umap.",max(dim.use),"PC.pdf"),width = 10,height = 10)
DotPlot(sce, features = unique(features)) + RotatedAxis()
dev.off()
#细胞类群相似树
sce= SetIdent(sce,value="celltype") 
sce <- BuildClusterTree(
  sce,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/celltype-ClusterTree_",max(dim.use),"PC.pdf"),width = 5,height = 4)
PlotClusterTree(sce)
dev.off()
# ############################################################################################################
# ## DE analysis
# ############################################################################################################

sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#sce.markers %>%
#  group_by(cluster) %>%
# slice_max(n = 2, order_by = avg_log2FC)
write.table(sce.markers,file=paste0("./",sam.name,"/","_marker_genes_umap_",max(dim.use),"PC.txt"),sep="\t",quote = F)
#4.存储数据大类sce1
save(sce,file=paste0("./",sam.name,"/","sce1.RData"))
#然后将sample和tissue加入到seurat对象的metadata里
#cell.name<-row.names(sce@meta.data)
#cell.name<-strsplit(cell.name,"_")
#sample<-c()
#tissue<-c()
#for (i in 1:length(cell.name)){sample<-c(sample,cell.name[[i]][3]);tissue<-c(tissue,cell.name[[i]][2])}
#sce@meta.data[["sample"]]<-sample
#sce@meta.data[["tissue"]]<-tissue
#rm(cell.name,sample,tissue,i)
table(sce@meta.data$celltype)
View(sce@meta.data)

#输出meta.data数据，在excel中修改补充
write.table(sce@meta.data,file=paste0("./",sam.name,"/","sce@meta.data",max(dim.use),"PC.txt"),sep="\t",quote = F)
#补充meta.data数据
ldf=read.table("patient.txt",header = T)[,1]
sce@meta.data$patient <- ldf
ldf=read.table("tissue.txt",header = T)[,1]
sce@meta.data$tissue <- ldf
ldf=read.table("age.txt",header = T)[,1]
sce@meta.data$age <- ldf
ldf=read.table("sex.txt",header = T)[,1]
sce@meta.data$sex<- ldf
ldf=read.table("sample.txt",header = T)[,1]
sce@meta.data$sample <- ldf
ldf=read.table("TNM_T.txt",header = T)[,1]
sce@meta.data$TNM_T <- ldf
ldf=read.table("TNM_N.txt",header = T)[,1]
sce@meta.data$TNM_N <- ldf
ldf=read.table("TNM_M.txt",header = T)[,1]
sce@meta.data$TNM_M <- ldf
ldf=read.table("stage.txt",header = T)[,1]
sce@meta.data$stage <- ldf
ldf=read.table("lymph_node_metastasis.txt",header = T)[,1]
sce@meta.data$lymph_node_metastasis <- ldf

sce@meta.data$sourse<-"GSE164522"
rm(ldf)
#4.存储数据大类sce1
save(sce,file=paste0("./",sam.name,"/","sce1.RData"))


## 2.重新读取数据
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = 'sce1.RData')
View(sce@meta.data)
#确定用于细胞分群的PC,放置数据的位置"test2"
dim.use <- 1:50
sam.name <- "test2"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

#按照数据来源分组展示细胞异同-sample
#sce= SetIdent(sce,value="subname")
pdf(paste0("./",sam.name,"/umap_patient.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="patient",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_tissue.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="tissue",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_TNM_T.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="TNM_T",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_TNM_N.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="TNM_N",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_TNM_M.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="TNM_M",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/umap_sex.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="sex",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/lymph_node_metastasis.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="lymph_node_metastasis",reduction='umap')
dev.off()


write.csv(table(sce@meta.data$celltype,sce@meta.data$tissue),file="celltype_tissue.csv")
write.csv(table(sce@meta.data$sample,sce@meta.data$celltype),file="sample_celltype.csv")
table(sce@meta.data$sample,sce@meta.data$celltype)
table(sce@meta.data$tissue)
table(sce@meta.data$celltype)
#4.提取指定单细胞亚群
sce_Macrophage<-subset(sce,celltype %in% c("Myeloid"))
sce<-sce_Macrophage
dim(sce)
table(sce@meta.data$celltype)
rm(sce_Macrophage)

#4.存储数据myloid
save(sce,file=paste0("./",sam.name,"/","sce2.RData"))


## 2.重新读取数据
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = 'sce2.RData')
View(sce@meta.data)
#确定用于细胞分群的PC,放置数据的位置"test1"
dim.use <- 1:50
sam.name <- "test2"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
## 重新聚类分组
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
#all.genes <- rownames(sce) 
#sce <- ScaleData(sce, features = all.genes, vars.to.regress = c("percent.mt"))
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
#pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
#DimPlot(sce, reduction = "pca")
#dev.off()
#pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 5,height = 4)
ElbowPlot(sce, reduction = "pca",ndims = 50)
#dev.off()
dim.use<-1:50
sce <- FindNeighbors(sce, dims = 1:50)
sce <- FindClusters(sce, resolution = 0.8 )
sce <- RunUMAP(sce, dims = 1:50)
pdf(paste0("./",sam.name,"/myeloid_umap",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(sce, reduction = 'umap',label=T)
dev.off()
pdf(paste0("./",sam.name,"/myeloid_marker.umap",max(dim.use),"PC.pdf"),width = 15,height = 14)
FeaturePlot(sce, features = c("LILRA4","CLEC9A","CD1C","LAMP3","SPP1","C1QC","FCN1","MKI67","ISG15","CD3D","MZB1","CD68","MRC1","MARCO","CD14","FCGR3A"),cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,reduction = "umap")
dev.off()

pdf(paste0("./",sam.name,"/myeloid_1.umap",max(dim.use),"PC.pdf"),width = 15,height = 14)
FeaturePlot(sce, features = c("LYZ","LILRA4","CLEC9A","CD1C","LAMP3","CD14","FCGR3A","SPP1","APOE","C1QC","NLRP3","CD163","S100A9","KIT","MKI67","CSF1R"),cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,reduction = "umap")
dev.off()

features= c("LILRA4","CLEC9A","CLEC10A","CD1C","LAMP3","SPP1","C1QC","FCN1","MKI67",
            "ISG15","CD3D","MZB1","EPCAM","CD68","MRC1","MARCO","CD14","FCGR3A",
            "LYZ","APOE","NLRP3","S100A9","FCGR3B","TPSAB1","CCR7","CXCL10")
## 气泡图
pdf(paste0("./",sam.name,"/DotPlot_myeloid_marker.umap",max(dim.use),"PC.pdf"),width = 10,height = 10)
DotPlot(sce, features = unique(features)) + RotatedAxis()
dev.off()

table(sce@meta.data$seurat_clusters,sce@meta.data$tissue)


#按照数据来源分组展示细胞异同-sample
#sce= SetIdent(sce,value="subname")
pdf(paste0("./",sam.name,"/myeloid_umap_patient.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="patient",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/myeloid_umap_tissue.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="tissue",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/myeloid_umap_TNM_T.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="TNM_T",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/myeloid_umap_TNM_N.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="TNM_N",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/myeloid_umap_sex.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="sex",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/myeloid_umap_TNM_M.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="TNM_M",reduction='umap')
dev.off()
#### 2.1 计算marker基因
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","myloid_marker_genes",max(dim.use),"PC.txt"),sep="\t",quote = F)

#细胞类群相似树
sce= SetIdent(sce,value="seurat_clusters") 
sce <- BuildClusterTree(
  sce,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/myeloid_seurat_resolution0.8_clustersClusterTree_",max(dim.use),"PC.pdf"),width = 5,height = 4)
PlotClusterTree(sce)
dev.off()

#以下重命名髓系亚群
#方法二：使用unname函数配合向量：
cluster5celltype <-c("0"="Doublets_T",
                      "1"="cDC2_CD1C",
                      "2"="cDC2_CD1C",
                      "3"="Macro_SPP1",
                      "4"="Mono",
                      "5"="Macro_FCN1",
                      "6"="Macro_FCN1",
                      "7"="cDC2_CD1C",
                      "8"="cDC1_CLEC9A",
                      "9"="cDC2_CD1C",
                      "10"="Mono",
                      "11"="Macro_FCN1",
                      "12"="Macro_C1QC",
                      "13"="Macro_C1QC",
                      "14"="Macro_FCN1",
                      "15"="pDC_LILRA4",
                      "16"="Macro_FCN1",
                      "17"="Mono",
                      "18"="Macro_FCN1",
                      "19"="Macro_MKI67",
                      "20"="Doublets_T",
                      "21"="Kupffer_cell",
                      "22"="cDC1_CLEC9A",
                      "23"="Doublets_B",
                      "24"="Doublets_T",
                      "25"="Doublets_T",
                      "26"="cDC3_LAMP3"
)
sce[['subcelltype']] = unname(cluster5celltype[sce@meta.data$seurat_clusters])
View(sce@meta.data) 

pdf(paste0("./",sam.name,"/myeloid_subcelltype_umap_",max(dim.use),".pdf"),width = 8,height = 7)
DimPlot(sce, reduction = "umap",group.by = 'subcelltype',label = T)
dev.off()
#### 2.1 计算marker基因
sce= SetIdent(sce,value="subcelltype") 
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","myloid_subcelltype_marker_genes",max(dim.use),"PC.txt"),sep="\t",quote = F)

table(sce@meta.data$subcelltype,sce@meta.data$tissue)
table(sce@meta.data$subcelltype)
table(sce@meta.data$tissue)
#4.存储数据大类sce3
save(sce,file=paste0("./",sam.name,"/","sce3.RData"))
## 2.重新读取数据
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = 'sce2.RData')
View(sce@meta.data)
#确定用于细胞分群的PC,放置数据的位置"test1"
dim.use <- 1:20
sam.name <- "test1"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
## 重新聚类分组
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
#all.genes <- rownames(sce) 
#sce <- ScaleData(sce, features = all.genes, vars.to.regress = c("percent.mt"))
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
#pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
#DimPlot(sce, reduction = "pca")
#dev.off()
#pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 5,height = 4)
ElbowPlot(sce, reduction = "pca",ndims = 50)
#dev.off()
dim.use<-1:20
sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce, resolution = 1 )
sce <- RunUMAP(sce, dims = 1:20)
pdf(paste0("./",sam.name,"/myeloid_umap",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(sce, reduction = 'umap',label=T)
dev.off()
pdf(paste0("./",sam.name,"/myeloid_marker.umap",max(dim.use),"PC.pdf"),width = 15,height = 14)
FeaturePlot(sce, features = c("LILRA4","CLEC9A","CD1C","LAMP3","SPP1","C1QC","FCN1","MKI67","ISG15","CD3D","MZB1","CD68","MRC1","MARCO","CD14","FCGR3A"),cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,reduction = "umap")
dev.off()

pdf(paste0("./",sam.name,"/myeloid_1.umap",max(dim.use),"PC.pdf"),width = 15,height = 14)
FeaturePlot(sce, features = c("LYZ","LILRA4","CLEC9A","CD1C","LAMP3","CD14","FCGR3A","SPP1","APOE","C1QC","NLRP3","CD163","S100A9","KIT","MKI67","CSF1R"),cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25") ,reduction = "umap")
dev.off()

features= c("LILRA4","CLEC9A","CLEC10A","CD1C","LAMP3","SPP1","C1QC","FCN1","MKI67",
            "ISG15","CD3D","MZB1","CD68","MRC1","MARCO","CD14","FCGR3A",
            "LYZ","APOE","NLRP3","S100A9","FCGR3B","TPSAB1","CCR7","CXCL10")
## 气泡图
pdf(paste0("./",sam.name,"/DotPlot_myeloid_marker.umap",max(dim.use),"PC.pdf"),width = 10,height = 10)
DotPlot(sce, features = unique(features)) + RotatedAxis()
dev.off()

table(sce@meta.data$seurat_clusters,sce@meta.data$tissue)


#按照数据来源分组展示细胞异同-sample
#sce= SetIdent(sce,value="subname")
pdf(paste0("./",sam.name,"/myeloid_umap_patient.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="patient",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/myeloid_umap_tissue.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="tissue",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/myeloid_umap_TNM_T.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="TNM_T",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/myeloid_umap_TNM_N.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="TNM_N",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/myeloid_umap_sex.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="sex",reduction='umap')
dev.off()
pdf(paste0("./",sam.name,"/myeloid_umap_MMR.",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = sce, group.by="MMR",reduction='umap')
dev.off()
#### 2.1 计算marker基因
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","myloid_marker_genes",max(dim.use),"PC.txt"),sep="\t",quote = F)

#细胞类群相似树
sce= SetIdent(sce,value="seurat_clusters") 
sce <- BuildClusterTree(
  sce,
  dims = dim.use,
  reorder = F,
  reorder.numeric = F)
pdf(paste0("./",sam.name,"/myeloid_seurat_resolution0.5_clustersClusterTree_",max(dim.use),"PC.pdf"),width = 5,height = 4)
PlotClusterTree(sce)
dev.off()

#以下重命名髓系亚群
#方法二：使用unname函数配合向量：
cluster5celltype <- c("0"="T-cell doublets",
                      "1"="Macro-SPP1",
                      "2"="Macro-FCN1",
                      "3"="Macro-SPP1",
                      "4"="Mono",
                      "5"="Macro-C1QC",
                      "6"="Macro-C1QC",
                      "7"="Macro-C1QC",
                      "8"="Macro-C1QC",
                      "9"="cDC2-CD1C",
                      "10"="Macro-C1QC",
                      "11"="Macro-SPP1",
                      "12"="cDC2-CD1C",
                      "13"="Macro-FCN1",
                      "14"="Mono",
                      "15"="pDC-LILRA4",
                      "16"="Macro-FCN1",
                      "17"="Macro-MKI67",
                      "18"="cDC2-CD1C",
                      "19"="Macro-C1QC",
                      "20"="Mono",
                      "21"="cDC1-CLEC9A",
                      "22"="cDC3-LAMP3",
                      "23"="Mono",
                      "24"="Mono",
                      "25"="Epithelial-cell doublets",
                      "26"="B-cell doublets"
)
sce[['subcelltype']] = unname(cluster5celltype[sce@meta.data$seurat_clusters])
View(sce@meta.data) 

pdf(paste0("./",sam.name,"/myeloid_subcelltype_umap_",max(dim.use),".pdf"),width = 8,height = 7)
DimPlot(sce, reduction = "umap",group.by = 'subcelltype',label = T)
dev.off()
#### 2.1 计算marker基因
sce= SetIdent(sce,value="subcelltype") 
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","myloid_subcelltype_marker_genes",max(dim.use),"PC.txt"),sep="\t",quote = F)

table(sce@meta.data$subcelltype,sce@meta.data$tissue)
table(sce@meta.data$subcelltype)
table(sce@meta.data$tissue)
write.csv(table(sce@meta.data$subcelltype,sce@meta.data$tissue),file="subcelltype_tissue.csv")
write.csv(table(sce@meta.data$sample,sce@meta.data$subcelltype),file="sample_subcelltype.csv")
table(sce@meta.data$sample,sce@meta.data$subcelltype)
#4.存储数据大类sce3
save(sce,file=paste0("./",sam.name,"/","sce3.RData"))

#####提取CRC髓系细胞

#4.提取指定单细胞亚群
sce_Macrophage<-subset(sce,tissue %in% c("CRC"))
sce<-sce_Macrophage
dim(sce)
rm(sce_Macrophage)

#4.存储数据myloid
save(sce,file=paste0("./",sam.name,"/","sce4.RData"))

## 2.重新读取数据
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = 'sce3.RData')
View(sce@meta.data)
#确定用于细胞分群的PC,放置数据的位置"test1"
dim.use <- 1:50
sam.name <- "test3"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
#4.提取指定单细胞亚群
sce= SetIdent(sce,value="subcelltype") 
sce_Macrophage<-subset(sce,subcelltype %in% c("Macro_C1QC","Macro_FCN1",
                                              "Macro_MKI67","Macro_SPP1",
                                              "Mono","Kupffer_cell"))
sce<-sce_Macrophage
dim(sce)
rm(sce_Macrophage)

#4.存储数据myloid
save(sce,file=paste0("./",sam.name,"/","5tissue_mono_macro.RData"))
#### 2.1 计算marker基因
sce= SetIdent(sce,value="tissue") 
all.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/","myloid_subcelltype_marker_genes",max(dim.use),"PC.txt"),sep="\t",quote = F)
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




