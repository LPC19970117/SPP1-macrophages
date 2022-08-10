#终端部分
# 创建名为cellphonedb的虚拟环境
conda create -n cellphonedb python=3.7 
# 激活虚拟环境
conda activate cellphonedb 
# 在虚拟环境中下载软件
pip install cellphonedb
# 假如网络不行，就加上  -i http://mirrors.aliyun.com/pypi/simple/   
# 很简单就安装成功， 试试看运行它 获取帮助信息
cellphonedb --help

#R部分由于细胞类型已经注释好，接下来准备cellphonedb的文件：
#表达谱文件cellphonedb_count.txt和细胞分组注释文件cellphonedb_meta.txt。
write.table(as.matrix(sce@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(sce@meta.data), sce@meta.data[,'cell_type', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA

write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

#终端
cellphonedb method statistical_analysis  cellphonedb_meta.txt  cellphonedb_count.txt      --counts-data=gene_name


#如果我们count的基因是基因名格式，需要添加参数--counts-data=gene_name，如果行名为ensemble名称的话，可以不添加这个参数，使用默认值即可。

#下面没卵用
conda activate cellphonedb
# 必须要保证当前路径下面有前面的步骤输出的out文件夹哦 
cellphonedb plot dot_plot 
cellphonedb plot heatmap_plot cellphonedb_meta.txt 


# 做完后为了跟前面的区分，我把 out文件夹，修改名字为 out-by-symbol 文件夹啦 

#R 热图 可视化
library(tidyverse)
library(RColorBrewer)
library(scales)

pvalues=read.table("./out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
pvalues=pvalues[,12:dim(pvalues)[2]] #此时不关注前11列
statdf=as.data.frame(colSums(pvalues < 0.05)) #统计在某一种细胞pair的情况之下，显著的受配体pair的数目；阈值可以自己选
colnames(statdf)=c("number")

#排在前面的分子定义为indexa；排在后面的分子定义为indexb
statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
#设置合适的细胞类型的顺序
rankname=sort(unique(statdf$indexa)) 
#转成因子类型，画图时，图形将按照预先设置的顺序排列
statdf$indexa=factor(statdf$indexa,levels = rankname)
statdf$indexb=factor(statdf$indexb,levels = rankname)

statdf%>%ggplot(aes(x=indexa,y=indexb,fill=number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,65))+
  scale_x_discrete("cluster 1 produces molecule 1")+
  scale_y_discrete("cluster 2 produces molecule 2")+
  theme_minimal()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45),
    panel.grid = element_blank()
  )
ggsave(filename = "interaction.num.1.pdf",device = "pdf",width = 12,height = 10,units = c("cm"))

#对称热图
library(tidyverse)
library(RColorBrewer)
library(scales)

pvalues=read.table("./out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
colnames(statdf)=c("number")

statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
statdf$total_number=0

for (i in 1:dim(statdf)[1]) {
  tmp_indexb=statdf[i,"indexb"]
  tmp_indexa=statdf[i,"indexa"]
  if (tmp_indexa == tmp_indexb) {
    statdf[i,"total_number"] = statdf[i,"number"]
  } else {
    statdf[i,"total_number"] = statdf[statdf$indexb==tmp_indexb & statdf$indexa==tmp_indexa,"number"]+
      statdf[statdf$indexa==tmp_indexb & statdf$indexb==tmp_indexa,"number"]
  }
}

rankname=sort(unique(statdf$indexa)) 
statdf$indexa=factor(statdf$indexa,levels = rankname)
statdf$indexb=factor(statdf$indexb,levels = rankname)

statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,120))+#limits=c(-5,100))能显著改变颜色
  scale_x_discrete("cluster 1")+
  scale_y_discrete("cluster 2")+
  theme_minimal()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45),
    panel.grid = element_blank()
  )
ggsave(filename = "interaction.num.2——.pdf",device = "pdf",width = 12,height = 10,units = c("cm"))

#气泡图
source("CCC.bubble.R")
CCC(
  pfile="./out/pvalues.txt",
  mfile="./out/means.txt",
  #neg_log10_th= -log10(0.05),
  #means_exp_log2_th=1,
  #neg_log10_th2=3,
  #means_exp_log2_th2=c(-4,6),
  #notused.cell=c("Bcell","Gcell"),
  #used.cell=c("Mcell"),
  #cell.pair=c("Mcell.Scell","Mcell.NKcell","Mcell.Tcell","Scell.Mcell","NKcell.Mcell","Tcell.Mcell"),#这里是自定义的顺序，若是可选细胞对的子集，则只展示子集，若有交集则只展示交集；空值情况下，会根据可选细胞对自动排序
  #gene.pair=c("MIF_TNFRSF14","FN1_aVb1 complex","EGFR_MIF")#作用同上
)
ggsave(filename = "interaction.detail.1.pdf",device = "pdf",width =20,height = 12,units = "cm")
#再细化
CCC(
  pfile="./out/pvalues.txt",
  mfile="./out/means.txt",
  cell.pair=c("Mcell.Scell","Mcell.NKcell","Mcell.Tcell","Scell.Mcell","NKcell.Mcell","Tcell.Mcell"),#这里是自定义的顺序，若是可选细胞对的子集，则只展示子集，若有交集则只展示交集；空值情况下，会根据可选细胞对自动排序
  gene.pair=c("MIF_TNFRSF14","FN1_aVb1 complex","EGFR_MIF")#作用同上
)
ggsave(filename = "interaction.detail.2.pdf",device = "pdf",width =16,height = 10,units = "cm")
