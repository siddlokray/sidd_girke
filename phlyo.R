##
### INSTALL/LOAD PACKAGES
##
install.packages('poolr')
library(systemPipeR);library(ShortRead);library(data.table);library(poolr);library(ape);library(dendextend);library(gplots)
homedir <- "/rhome/slokr001/shared/slokray"
setwd(homedir)

##
### SET ARGS
##
parampath <- "/rhome/slokr001/shared/slokray/param/hisat2PE.param"
targets <- read.delim("/rhome/slokr001/shared/slokray/SraRunInfo - pairedTargets.txt", comment.char ="#")
args <- systemArgs(sysma=parampath, mytargets="/rhome/slokr001/shared/slokray/SraRunInfo - pairedTargets.txt")

##
### PVALUE
##

##
### PHYLO FOR FPKM
##
factor <- factor(targetsin(args)$Factor)
rpkmDFeBygpath <- "./mousefpkmDFeByg.xls"
rpkmDFeByg <- read.table(rpkmDFeBygpath, check.names=FALSE)
tableEdgeR <- read.table('./mouseedgeRcomp.xls')
pvalColumns <- c("CR.M-CR.F_PValue", "CR.M-SnellDwarfMice.M_PValue", "CR.M-MethionineRestriction.M_PValue", "CR.F-SnellDwarfMice.M_PValue", 
                 "CR.F-MethionineRestriction.M_PValue", "SnellDwarfMice.M-MethionineRestriction.M_PValue", "CR.M-Rapamycin.M_PValue", 
                 "CR.M-GHRKO.M_PValue", "Rapamycin.M.GHRKO.M_PValue")
tablePval <- tableEdgeR[, seq(4, ncol(tableEdgeR), by = 5)]
rpkmDFeByg <- rpkmDFeByg[complete.cases(tableEdgeR), ]
tablePval <- na.omit(tablePval)
for (i in 1:nrow(tablePval)) {
  rowPval <- tablePval[i, ]
  rowPval <- c(as.numeric(rowPval))
  tablePval$fisherPval[i] <- as.numeric(min(rowPval))
}
filteredCounts <- rpkmDFeByg[tablePval$fisherPval < 0.05,]
d <- cor(filteredCounts, method="spearman")
hc <- hclust(as.dist(1-d))
dend <- as.dendrogram(hc)
rank <- rank_branches(dend)
# branches 12
colorCode <- c(Acarbose.M='cadetblue4', Acarbose.F='cadetblue1', alphaestradiol.M='deepskyblue4', alphaestradiol.F='deepskyblue1', 
                  Control.M='gray30', Control.F='gray70', CR.M='steelblue4', CR.F='steelblue1', GHRKO.M='darkgreen', 
                  GHRKOControl.M='olivedrab4', MethionineRestriction.M='purple4', MRControl.M='orchid4', Protandim.M='blue4', 
                  Protandim.F='lightcyan1', Rapamycin.M='cyan4', Rapamycin.F='lightcyan2', SnellDwarfMice.M='brown4', 
                  SnellDwarfMiceControl.M='coral4')
tipColor <- colorCode[gsub('.{2}$', '', labels(dend))]
labels_colors(dend) <- tipColor
plot(dend, horiz=TRUE)

##
### PHYLO FOR RLOG COUNTS
##
tableEdgeR <- read.table('./mouseedgeRcomp.xls')
rlogCounts <- read.table('./rlog.xls')
pvalColumns <- c("CR.M-CR.F_PValue", "CR.M-SnellDwarfMice.M_PValue", "CR.M-MethionineRestriction.M_PValue", "CR.F-SnellDwarfMice.M_PValue", 
                 "CR.F-MethionineRestriction.M_PValue", "SnellDwarfMice.M-MethionineRestriction.M_PValue", "CR.M-Rapamycin.M_PValue", 
                 "CR.M-GHRKO.M_PValue", "Rapamycin.M.GHRKO.M_PValue")
tablePval <- tableEdgeR[,pvalColumns]
rlogCounts <- rlogCounts[complete.cases(tableEdgeR), ]
tablePval <- na.omit(tablePval)
for (i in 1:nrow(tablePval)) {
  rowPval <- tablePval[i, ]
  rowPval <- c(as.numeric(rowPval))
  tablePval$fisherPval[i] <- as.numeric(min(rowPval))
}
filteredCounts <- rlogCounts[tablePval$fisherPval < 0.05,]
d <- cor(filteredCounts, method="spearman")
hc <- hclust(as.dist(1-d))
dend <- as.dendrogram(hc)
rank <- rank_branches(dend)
# branches 13
colorCode <- c(Acarbose.M='cadetblue4', Acarbose.F='cadetblue1', alphaestradiol.M='deepskyblue4', alphaestradiol.F='deepskyblue1', 
               Control.M='gray30', Control.F='gray70', CR.M='steelblue4', CR.F='steelblue1', GHRKO.M='darkgreen', 
               GHRKOControl.M='olivedrab4', MethionineRestriction.M='purple4', MRControl.M='orchid4', Protandim.M='blue4', 
               Protandim.F='lightcyan1', Rapamycin.M='cyan4', Rapamycin.F='lightcyan2', SnellDwarfMice.M='brown4', 
               SnellDwarfMiceControl.M='coral4')
tipColor <- colorCode[gsub('.{2}$', '', labels(dend))]
labels_colors(dend) <- tipColor
plot(dend, horiz=TRUE)

##
### PHYLO FOR SCALED COUNTS
##
tableEdgeR <- read.table('./mouseedgeRcomp.xls')
scaledCounts <- read.table('./scale.xls')
pvalColumns <- c("CR.M-CR.F_PValue", "CR.M-SnellDwarfMice.M_PValue", "CR.M-MethionineRestriction.M_PValue", "CR.F-SnellDwarfMice.M_PValue", 
                 "CR.F-MethionineRestriction.M_PValue", "SnellDwarfMice.M-MethionineRestriction.M_PValue", "CR.M-Rapamycin.M_PValue", 
                 "CR.M-GHRKO.M_PValue", "Rapamycin.M.GHRKO.M_PValue")
tablePval <- tableEdgeR[,pvalColumns]
scaledCounts <- scaledCounts[complete.cases(tableEdgeR), ]
tablePval <- na.omit(tablePval)
for (i in 1:nrow(tablePval)) {
  rowPval <- tablePval[i, ]
  rowPval <- c(as.numeric(rowPval))
  tablePval$fisherPval[i] <- as.numeric(min(rowPval))
}
filteredCounts <- scaledCounts[tablePval$fisherPval < 0.05,]
d <- cor(filteredCounts, method="spearman")
hc <- hclust(as.dist(1-d))
dend <- as.dendrogram(hc)
rank <- rank_branches(dend)
# branches 14
colorCode <- c(Acarbose.M='cadetblue4', Acarbose.F='cadetblue1', alphaestradiol.M='deepskyblue4', alphaestradiol.F='deepskyblue1', 
               Control.M='gray30', Control.F='gray70', CR.M='steelblue4', CR.F='steelblue1', GHRKO.M='darkgreen', 
               GHRKOControl.M='olivedrab4', MethionineRestriction.M='purple4', MRControl.M='orchid4', Protandim.M='blue4', 
               Protandim.F='lightcyan1', Rapamycin.M='cyan4', Rapamycin.F='lightcyan2', SnellDwarfMice.M='brown4', 
               SnellDwarfMiceControl.M='coral4')
tipColor <- colorCode[gsub('.{2}$', '', labels(dend))]
labels_colors(dend) <- tipColor
plot(dend, horiz=TRUE)
