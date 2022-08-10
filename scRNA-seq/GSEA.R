library(clusterProfiler) 
library(GSVA)
library(pheatmap)
library(gplots)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
expr <- read.csv("easy_input_expr.csv", check.names = F, row.names = 1)
expr[1:3,1:3]
dim(expr)
expr<-as.matrix(expr)
subt <- read.table("easy_input_subtype.txt", sep = "\t", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
head(subt)
table(subt$TCGA_Subtype)
#########如果是直接来自单细胞的数据
library(Seurat)
load(file = 'sce4.RData')
write.table(cg_sce@meta.data,file=paste0("./","sce@meta.data",".txt"),sep="\t",quote = F)
## 重新聚类分组
cg_sce <- NormalizeData(cg_sce, normalization.method = "LogNormalize", scale.factor = 10000) 
#sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
all.genes <- rownames(cg_sce) 
#sce <- ScaleData(cg_sce, features = all.genes, vars.to.regress = c("percent.mt"))
expr <-as.matrix(cg_sce@assays$RNA@data)
cg_sce
expr<-as.matrix(expr)#加快分析速度
subt <- read.table("easy_input_subtype.txt", sep = "\t", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
head(subt)
table(subt$TCGA_Subtype)

### 自定义函数显示进度 ###
display.progress = function ( index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
} 

n.sub <- length(table(subt$TCGA_Subtype)) # 亚型个数
n.sub.label <- unique(subt$TCGA_Subtype) # 亚型名称

# expr <- log2(expr + 1) 
treat_list <- ctrl_list <- degs.list <- list() # 初始化列表
for (i in 1:n.sub) {
  cat(paste0(n.sub.label[i], " vs. Others starts!\n"))
  treat_list[[i]] <- rownames(subt)[which(subt$TCGA_Subtype == n.sub.label[i])] # 选取某亚型为treat组
  ctrl_list[[i]] <- rownames(subt)[-which(subt$TCGA_Subtype == n.sub.label[i])] # 选取剩余亚型为control组
  
  meanA <- meanB <- p <- fc <- lgfc <- c() #初始化向量
  for (k in 1:nrow(expr)) {
    display.progress(index = k,totalN = nrow(expr)) #显示进程
    a <- as.numeric(expr[k,treat_list[[i]]])
    b <- as.numeric(expr[k,ctrl_list[[i]]])
    p <- c(p,t.test(a,b,na.rm=T)$p.value) # 检验表达差异
    meanA <- c(meanA,mean(a)) # treat组基因k表达均值
    meanB <- c(meanB,mean(b)) # control组基因k表达均值
    fc <- c(fc,mean(a)/mean(b)) # 计算foldchange
    lgfc <- c(lgfc,log2(mean(a)/mean(b))) # 计算log2foldchange
  }
  fdr <- p.adjust(p,method = "fdr") # 校正p值
  
  # 生成差异表达结果，其中log2FoldChange, pvalue, padj模仿DESeq2结果格式
  # 由于差异表达分析不是本众筹的目的，所以这里采用简单的两样本t检验寻找显著差异表达基因。差异表达分析可根据实际数据情况换用limma（例文）, DESeq, DESeq2, edgeR等方法。
  degs <- data.frame(mean_treat=meanA,
                     mean_ctrl=meanB,
                     FoldChange=fc,
                     log2FoldChange=lgfc,
                     pvalue=p,
                     padj=fdr,
                     row.names = rownames(expr),
                     stringsAsFactors = F)
  
  write.table(degs,paste0(n.sub.label[[i]],"_degs.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
  
  degs.list[[n.sub.label[i]]] <- as.data.frame(na.omit(degs))
  
  cat("\n")
}


subtype_specific_gsea <- function(msigdb=NULL,n.top=10,mode=c("up","down"),degs.list=NULL,subtype.label=NULL,nPerm.gsea=1000,minGSSize.gsea=10,maxGSSize.gsea=500,pvalueCutoff.gsea=1){
  
  MSigDB <- read.gmt(msigdb)
  GSEA.list <- top.gs <- list() #初始化结果列表
  
  if(!is.element(mode, c("up", "dn"))) { stop("mode must be up or dn!\n") }
  
  for (i in 1:n.sub) {
    degs <- degs.list[[n.sub.label[i]]]
    ### 去除无穷值
    degs <- degs[!is.infinite(degs$log2FoldChange),]
    geneList <- degs$log2FoldChange; names(geneList) <- rownames(degs)
    geneList <- sort(geneList,decreasing = T) # ranked gene set
    
    # 由于GSEA不可重复，所以保存GSEA对象入列表，方便下次调用
    cat(paste0("GSEA for ",subtype.label[i]," starts!\n"))
    GSEA.list[[subtype.label[i]]] <- GSEA(geneList = geneList,
                                          TERM2GENE=MSigDB,
                                          nPerm = nPerm.gsea,
                                          minGSSize = minGSSize.gsea,
                                          maxGSSize = maxGSSize.gsea,
                                          seed = T,
                                          verbose = F,
                                          pvalueCutoff = pvalueCutoff.gsea) # 输出全部的GESA结果
    
    GSEA.dat <- as.data.frame(GSEA.list[[subtype.label[i]]])
    
    if(mode == "up") {
      GSEA.dat <- GSEA.dat[order(GSEA.dat$NES,decreasing = T),] # 根据NES降序排列，也就是找特异性上调通路
    } else {
      GSEA.dat <- GSEA.dat[order(GSEA.dat$NES,decreasing = F),] # 根据NES升序排列，也就是找特异性下调通路
    }
    
    # 输出每一次GSEA结果
    write.table(GSEA.dat,paste0(subtype.label[[i]],"_degs_",mode,"_gsea.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
    
    # 亚型特异性top基因集保存入列表
    top.gs[[subtype.label[i]]] <- rownames(GSEA.dat)[1:n.top] 
  }
  
  # 构建GSVA分析需要的gene list
  gs <- list()
  for (i in as.character(unlist(top.gs))) {
    gs[[i]] <- MSigDB[which(MSigDB$term %in% i),"gene"]
  }
  
  return(list(mode=mode,top.gs=top.gs,gs=gs))
}

msigdfFile = "c5.all.v7.5.1.symbols.gmt"
n.top = 5
mode = "up" #"up"和"dn"二选一
gs.up <- subtype_specific_gsea(msigdb = msigdfFile,
                               n.top = n.top,
                               degs.list = degs.list,
                               subtype.label = n.sub.label,
                               mode = mode)

# 计算GSVA得分
gsva_gs.up <- gsva(as.matrix(expr), gs.up$gs, method="gsva") 
dim(gsva_gs.up)

# 这里是30条通路，说明top通路无重叠

# 每个亚型取均值（也可以换用其他统计量，比如中位数等等）
gsva_gs.up_mean <- data.frame(row.names = rownames(gsva_gs.up)) 
for (i in n.sub.label) {
  gsva_gs.up_mean <- cbind.data.frame(gsva_gs.up_mean,
                                      data.frame(rowMeans(gsva_gs.up[,rownames(subt)[which(subt$TCGA_Subtype == i)]])))
}
colnames(gsva_gs.up_mean) <- n.sub.label
#自定义分组的颜色
jco <- c("#F2CCCC","#E6D8CF","#D5E3F0","#FDE7DA","#E2D6EC", "#CCEFDB","#FFDAB9","#6495ED")

annRows <- data.frame(subtype=rep(n.sub.label,each=n.top), names = unlist(gs.up$top.gs), stringsAsFactors = F)
annRows <- annRows[!duplicated(annRows$names),]; rownames(annRows) <- annRows$names # 倘若出现一条通路在>=2个亚型中上调，去掉重复值，这种情况在亚型较多的时候会发生

#示例数据是3个分组，有更多组就继续往后添加
annColors <- list(subtype=c("Macro-FCN1"=jco[1],
                            "Macro-SPP1"=jco[2],
                            "Mono"=jco[3],
                            "Macro-C1QC"=jco[4],
                            "Macro-MKI67"=jco[5],
                            "Kupffer cells"=jco[6]#,
                            #"Macro-MKI67"=jco[7],
                            #"pDC-LILRA4"=jco[8]
                            ))

filename <- paste0("subtype_specific_top_",mode,"_gsea.pdf")
pheatmap(gsva_gs.up_mean[rownames(annRows),],
         cellwidth = 10, cellheight = 10,
         #color = bluered(64), #自定义颜色
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA, #如果想要边框，就去掉这行
         annotation_row = annRows[,"subtype",drop = F],
         annotation_colors = annColors,
         filename = filename)

sessionInfo()

