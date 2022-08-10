###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######??????ѧ??: http://www.biowolf.cn/
######QQ??2749657388
######????QȺ??219795300
######΢??: 18520221056

####install.packages("sm")
####install.packages("vioplot")
library(vioplot)                                                    #???ð?
setwd("C:\\Users\\lexb4\\Desktop\\TCGAimmune\\11.vioplot")          #???ù???Ŀ¼
normal=110
tumor=60                                                           #??????Ʒ??Ŀ

rt=read.table("CIBERSORT.filter.txt",sep="\t",header=T,row.names=1,check.names=F)   #??ȡ?????ļ?

pdf("vioplot.pdf",height=8,width=15)              #????ͼƬ???ļ?????
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

#??ÿ??????ϸ??ѭ????????vioplot????????��ɫ??ʾ???????ú?ɫ??ʾ
for(i in 1:ncol(rt)){
  normalData=rt[1:normal,i]
  tumorData=rt[(normal+1):(normal+tumor),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'RoyalBlue')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'Firebrick3')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.001,paste0("p<0.001"),paste0("p=",p)),cex = 0.8)
  text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
}
dev.off()

###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######??????ѧ??: http://www.biowolf.cn/
######QQ??2749657388
######????QȺ??219795300
######΢??: 18520221056