#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma", version = "3.8")


library(limma)
library(estimate)
setwd("C:\\Users\\lexb4\\Desktop\\tcgaEstimate\\05.estimate")           #???ù???Ŀ¼
inputFile="normalize.txt"                                                  #?????ļ?????

#??ȡ?ļ?
rt=read.table(inputFile,sep="\t",header=T,check.names=F)

#????һ??????ռ?˶??У?ȡ??ֵ
#rt=as.matrix(rt)
#rownames(rt)=rt[,1]
#exp=rt[,2:ncol(rt)]
#dimnames=list(rownames(exp),colnames(exp))
#data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#data=avereps(data)
#group=sapply(strsplit(colnames(data),"\\-"),"[",4)
#group=sapply(strsplit(group,""),"[",1)
#group=gsub("2","1",group)
#ata=data[,group==0]
#out=data[rowMeans(data)>0,]
#out=rbind(ID=colnames(out),out)
#???????????ľ????ļ?
#write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

#????estimate??
filterCommonGenes(input.f="normalize.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct", 
              platform="illumina")

#????ÿ????Ʒ?Ĵ???
scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="scores.txt",sep="\t",quote=F,col.names=F)


######??????ѧ??: http://study.163.com/u/biowolf
######??????ѧ??: https://shop119322454.taobao.com
######??????ѧ??: http://www.biowolf.cn/
######???????䣺2740881706@qq.com
######????΢??: seqBio
######QQȺ:  259208034
