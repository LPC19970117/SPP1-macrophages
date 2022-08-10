
# monocle2
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
exp.matrix<-as(as.matrix(cg_sce@assays$RNA@data), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann)<-rownames(exp.matrix)
exp_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-cg_sce@meta.data
rownames(sample_ann)<-colnames(exp.matrix)
exp_pd<-new("AnnotatedDataFrame", data =sample_ann)

#生成monocle对象
exp.monocle<-newCellDataSet(exp.matrix,phenoData =exp_pd,featureData =exp_fd,expressionFamily=negbinomial.size())
head(pData(exp.monocle))
head(fData(exp.monocle))

#计算sizefactor
exp.monocle <- estimateSizeFactors(exp.monocle)
exp.monocle <- estimateDispersions(exp.monocle, cores=8)


#根据seurat cluster计算差异表达基因并挑选用于构建拟时序轨迹的基因
table(cg_sce$subcelltype)
if(F){
  ex1 <- c(cc.genes$s.genes, cc.genes$g2m.genes)
  ex2 <- grep('^MT-', rownames(cg_sce), value = T)
  ex3 <- grep('^RP[LS]', rownames(cg_sce), value = T)
  ex <- unique(c(ex1, ex2, ex3))
  #确定需要加入的基因
  #include <- c('CD3D','CD4','CD8A')
}
#order.genes <-  VariableFeatures()
#order.genes <- order.genes[!order.genes %in% ex]
#exp.monocle <- setOrderingFilter(exp.monocle, order.genes)
#plot_ordering_genes(exp.monocle)

if(F){
  disp_table <- dispersionTable(exp.monocle)
  order.genes <- subset(disp_table, mean_expression >= 0.01 & 
                          dispersion_empirical >= 1 * dispersion_fit) %>% 
    pull(gene_id) %>% as.character()
  #order.genes <- order.genes[!order.genes %in% ex]
  exp.monocle <- setOrderingFilter(exp.monocle, order.genes)
  plot_ordering_genes(exp.monocle)
}
#DDRTree方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
save(exp.monocle, file = "exp.monocle.Rdata")
#画图
plot_cell_trajectory(exp.monocle,color_by = "subcelltype")
ggsave("celltype.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(exp.monocle,color_by = "tissue")
ggsave("tissue.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(exp.monocle,color_by = "State")
ggsave("State.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(exp.monocle,color_by = "Pseudotime")
ggsave("Pseudotime.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(exp.monocle,color_by = "subcelltype")+facet_wrap(~subcelltype,nrow=2)
ggsave("subcelltypeb.pdf",device = "pdf",width = 21,height = 9,units = c("cm"))

s.genes <- c("FCN1","C1QC","SPP1","MKI67","CD14","MACRO")
plot_genes_violin(exp.monocle[s.genes,], grouping = "State", color_by = "State")
plot5 <-plot_genes_in_pseudotime(exp.monocle[s.genes,], color_by = "State")
pdf("plot5_2 ---.pdf",width = 10,height = 5)
print(plot5 )
