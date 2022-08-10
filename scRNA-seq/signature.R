
#4.计算signature
library(tidyverse)
library(Matrix)
library(cowplot)
inflammatory_gene <- readxl::read_xlsx("inflammatory_gene.xlsx")
View(inflammatory_gene)
#转换成list
gene <- as.list(inflammatory_gene)
sce <- AddModuleScore(
  object = sce,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'Inflammatory_Score'
)
sce= SetIdent(sce,value="subcelltype")
colnames(sce@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"   
# [4] "percent.mt"      "RNA_snn_res.0.5" "seurat_clusters"
# [7] "cell_type"       "CD_Features1" 
#colnames(sce@meta.data)["26"] <- 'Inflammatory_Score' 
pdf(paste0("./","Inflammatory_Score",".pdf"),width = 5,height = 4)
VlnPlot(sce,features = 'Inflammatory_Score1')
dev.off()

library(tidyverse)
library(Matrix)
library(cowplot)
M1_score <- readxl::read_xlsx("M1_score.xlsx")
View(M1_score)
#转换成list
gene <- as.list(M1_score)
sce <- AddModuleScore(
  object = sce,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'M1_score'
)
sce= SetIdent(sce,value="subcelltype")
pdf(paste0("./","M1_score1",".pdf"),width = 5,height = 4)
VlnPlot(sce,features = 'M1_score1')
dev.off()

library(tidyverse)
library(Matrix)
library(cowplot)
M2_score <- readxl::read_xlsx("M2_score.xlsx")
View(M1_score)
#转换成list
gene <- as.list(M2_score)
sce <- AddModuleScore(
  object = sce,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'M2_score'
)
#sce= SetIdent(sce,value="subcelltype")

sce= SetIdent(sce,value="subcelltype")
pdf(paste0("./","M2_score1",".pdf"),width = 5,height = 4)
VlnPlot(sce,features = 'M2_score1')
dev.off()


anti_inflammatory_gene <- readxl::read_xlsx("anti-inflammatory_gene.xlsx")
#转换成list
gene <- as.list(anti_inflammatory_gene)
sce <- AddModuleScore(
  object = sce,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'anti_inflammatory_gene'
)
#sce= SetIdent(sce,value="subcelltype")
pdf(paste0("./","anti_inflammatory_gene",".pdf"),width = 5,height = 4)
VlnPlot(sce,features = 'anti_inflammatory_gene1')
dev.off()


library(tidyverse)
library(Matrix)
library(cowplot)
proliferation_score <- readxl::read_xlsx("proliferation_score.xlsx")
View(proliferation_score)
#转换成list
#转换成list
gene <- as.list(proliferation_score)
sce <- AddModuleScore(
  object = sce,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'proliferation_score'
)
#sce= SetIdent(sce,value="subcelltype")
pdf(paste0("./","proliferation_score",".pdf"),width = 5,height = 4)
VlnPlot(sce,features = 'proliferation_score1')
dev.off()

library(tidyverse)
library(Matrix)
library(cowplot)
Angiogenesis <- readxl::read_xlsx("Angiogenesis.xlsx")
View(Angiogenesis)
#转换成list
#转换成list
gene <- as.list(Angiogenesis)
sce <- AddModuleScore(
  object = sce,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'Angiogenesis'
)
#sce= SetIdent(sce,value="subcelltype")
pdf(paste0("./","Angiogenesis",".pdf"),width = 5,height = 4)
VlnPlot(sce,features = 'Angiogenesis1')
dev.off()

Phageocytosis <- readxl::read_xlsx("Phageocytosis.xlsx")
View(Phageocytosis)
#转换成list
#转换成list
gene <- as.list(Phageocytosis)
sce <- AddModuleScore(
  object = sce,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'Phageocytosis'
)
#sce= SetIdent(sce,value="subcelltype")
pdf(paste0("./","Phageocytosis",".pdf"),width = 5,height = 4)
VlnPlot(sce,features = 'Phageocytosis1')
dev.off()

Liver_metastasis_genes
Liver_metastasis_genes <- readxl::read_xlsx("Liver_metastasis_genes.xlsx")
View(Liver_metastasis_genes)
#转换成list
#转换成list
gene <- as.list(Liver_metastasis_genes)
sce <- AddModuleScore(
  object = sce,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'Liver_metastasis_genes'
)
#sce= SetIdent(sce,value="subcelltype")
pdf(paste0("./","Liver_metastasis_genes",".pdf"),width = 5,height = 4)
VlnPlot(sce,features = 'Liver_metastasis_genes1')
dev.off()

view(sce@meta.data)
