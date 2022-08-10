###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######QQ：2749657388
######交流Q群：219795300
######微信: 18520221056

#install.packages('e1071')
#install.packages('parallel')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore", version = "3.8")


setwd("C:\\Users\\lexb4\\Desktop\\GEOimmune\\06.CIBERSORT")
source("GEOimmune.CIBERSORT.R")
results=CIBERSORT("ref.txt", "normalize.txt", perm=100, QN=TRUE)

###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######QQ：2749657388
######交流Q群：219795300
######微信: 18520221056