library(systemPipeR);library(ShortRead);library(data.table)
mouse_edgeDF <- read.table("./mouseDESeq2comp.xls")
# plot DEG Counts
DEG_list <- filterDEGs(degDF = mouse_edgeDF, filter = c(Fold = 2, FDR = 20))
# venn diagram upregulated
vennsetup <- overLapper(DEG_list$Up[c(13,4,3,8)], type = "vennsets")
vennPlot(vennsetup, mymain = "", mysub = "", colmode = 1)
# venn diagram downregulated
vennsetup <- overLapper(DEG_list$Down[c(13,4,3,8)], type = "vennsets")
vennPlot(vennsetup, mymain = "", mysub = "", colmode = 1)

