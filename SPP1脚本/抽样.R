# 每个细胞亚群抽50 
sce= SetIdent(sce,value="subcelltype")
allCells=names(Idents(sce))
allType = levels(Idents(sce))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce)== x ]
  cg=sample(cgCells,300)
  cg
}))
cg_sce = sce[, allCells %in% choose_Cells]
cg_sce
as.data.frame(table(Idents(cg_sce)))