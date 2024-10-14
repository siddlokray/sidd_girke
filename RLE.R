library(systemPipeR);library(ShortRead);library(data.table)
library(edgeR)
homedir <- "/rhome/slokr001/shared/slokray"
setwd(homedir)

parampath <- "/rhome/slokr001/shared/slokray/param/hisat2PE.param" ##system.file("extdata", "hisathuman.param", package="systemPipeR")
read.delim(parampath, comment.char = "#")
targets <- read.delim("/rhome/slokr001/shared/slokray/sidd_girke/SraRunInfo - pairedTargets.txt", comment.char ="#")
args <- systemArgs(sysma=parampath, mytargets="/rhome/slokr001/shared/slokray/sidd_girke/SraRunInfo - pairedTargets.txt")
args
modules(args)
cores(args)
names(args)
sysargs(args)[1]
outpaths(args)[1]
normalizePath(results(args))


rawCounts = read.table('./mousecountDFeByg.xls')

## dge <- DGEList(counts=rawCounts)
## dge = calcNormFactors(dge, method='RLE')
## rleCounts = cpm(dge, normalized.lib.sizes=TRUE)

mouse_cmp <- readComp(args, format="matrix", delim="-")
mouse_edgeDF <- run_edgeR(countDF=rawCounts, 
                          targets=targetsin(args), 
                          cmp=mouse_cmp[[1]], 
                          independent=FALSE, 
                          mdsplot="")
write.table(mouse_edgeDF, "./mouseedgeRcounts.csv", 
            quote=FALSE, sep="\t", col.names = NA)


library(biomaRt)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
deseq2_data = fread('mousecount_edgeR_UDPATED.xls')
rownames = deseq2_data[,'V1']
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'ensembl_gene_id',
                   values = rownames,
                   mart = ensembl)
gene_info <- as.data.table(gene_info)

for (i in deseq2_data$V1) {
  index1 = which(deseq2_data$V1 == i)
  index2 = which(gene_info$ensembl_gene_id == i)
  if (length(index2) != 0) {
    deseq2_data$V1[index1] = gene_info$external_gene_name[index2]
  }
}

fwrite(deseq2_data, 'mouseedgeRcounts_UPDATED.csv')



library(GOstats); library(GO.db)
library(DESeq2); library(clusterProfiler); library(org.Mm.eg.db); library(GO.db); library(pathview); library(readr);library(fgsea)
library(tidyverse)
deseq <- read_csv('mouseedgeRcounts_UPDATED.csv', show_col_types = FALSE)
deseq <- deseq[complete.cases(deseq), ]
pvalColumns <- c("Control.F.12-Acarbose.F.12_PValue", "Control.M.12-Acarbose.M.12_PValue", "Control.F.12-CR.F.12_PValue", "Control.M.12-CR.M.12_PValue", 
                  "GHRKOControl.M-GHRKO.M_PValue", "Control.F.12-Rapamycin.F.12_PValue",  "Control.M.12-Rapamycin.M.12_PValue", 
                  "SnellDwarfMiceControl.M-SnellDwarfMice.M_PValue", "Control.F.6-Acarbose.F.6_PValue", "Control.M.6-Acarbose.M.6_PValue", 
                  "Control.F.6-CR.F.6_PValue", "Control.M.6-CR.M.6_PValue", "Control.F.6-Rapamycin.F.6_PValue",  "Control.M.6-Rapamycin.M.6_PValue", 
                  "Control.F.6-Protandim.F.6_PValue", "Control.F.6-Protandim.M.6_PValue", "Control.F.6-alphaestradiol.F.6_PValue", 
                  "Control.M.6-alphaestradiol.M.6_PValue", "MRControl.M-MethionineRestriction.M_PValue")
statColumns <- c("Control.F.12-Acarbose.F.12_logFC", "Control.M.12-Acarbose.M.12_logFC", "Control.F.12-CR.F.12_logFC", "Control.M.12-CR.M.12_logFC", 
                 "GHRKOControl.M-GHRKO.M_logFC", "Control.F.12-Rapamycin.F.12_logFC",  "Control.M.12-Rapamycin.M.12_logFC", 
                 "SnellDwarfMiceControl.M-SnellDwarfMice.M_logFC", "Control.F.6-Acarbose.F.6_logFC", "Control.M.6-Acarbose.M.6_logFC", 
                 "Control.F.6-CR.F.6_logFC", "Control.M.6-CR.M.6_logFC", "Control.F.6-Rapamycin.F.6_logFC",  "Control.M.6-Rapamycin.M.6_logFC", 
                 "Control.F.6-Protandim.F.6_logFC", "Control.F.6-Protandim.M.6_logFC", "Control.F.6-alphaestradiol.F.6_logFC", 
                 "Control.M.6-alphaestradiol.M.6_logFC", "MRControl.M-MethionineRestriction.M_logFC")
tablePval <- deseq[,pvalColumns]
for (i in 1:nrow(tablePval)) {
  rowPval <- tablePval[i, ]
  rowPval <- c(as.numeric(rowPval))
  tablePval$minPval[i] <- min(p.adjust(rowPval, method = "BH"))
}
deseq <- deseq[tablePval$minPval < 0.05, ]

load_keggList <- function(org="mmu") {
  suppressMessages(suppressWarnings(library(KEGG.db))) 
  kegg_gene_list <- as.list(KEGGPATHID2EXTID) # All organisms in kegg
  kegg_gene_list <- kegg_gene_list[grepl(org, names(kegg_gene_list))] # Only human
  kegg_name_list <- unlist(as.list(KEGGPATHID2NAME)) # All organisms in kegg
  kegg_name_list <- kegg_name_list[gsub(paste0("^", org), "", names(kegg_gene_list))]
  names(kegg_gene_list) <- paste0(names(kegg_gene_list), " (", names(kegg_name_list), ") - ", kegg_name_list)
  return(kegg_gene_list)
}
keggdb <- load_keggList(org="mmu")


load_reacList <- function(org="R-MMU") {
  library(reactome.db)
  reac_gene_list <- as.list(reactomePATHID2EXTID) # All organisms in reactome
  reac_gene_list <- reac_gene_list[grepl(org, names(reac_gene_list))] # Only human
  reac_name_list <- unlist(as.list(reactomePATHID2NAME)) # All organisms in reactome
  reac_name_list <- reac_name_list[names(reac_gene_list)]
  names(reac_gene_list) <- paste0(names(reac_gene_list), " (", names(reac_name_list), ") - ", gsub("^.*: ", "", reac_name_list))
  return(reac_gene_list)
}
reacdb <- load_reacList(org="R-MMU")

load_goList <- function(ont="BP", org="mmu") {
  go_gene_mapping <- bitr(keys(org.Mm.eg.db, keytype = "ENTREZID"),
                          fromType = "ENTREZID", 
                          toType = "GOALL", 
                          OrgDb = org.Mm.eg.db)
  go_gene_mapping <- go_gene_mapping[go_gene_mapping$ONTOLOGYALL == ont, ]
  go_gene_list <- split(go_gene_mapping$ENTREZID, go_gene_mapping$GOALL)
  go_names <- Term(GO.db::GOTERM[names(go_gene_list)])
  names(go_gene_list) <- paste0(names(go_gene_list), " - ", go_names)
  
  return(go_gene_list)
}

godb <- load_goList(ont="BP", org="mmu")

gene_ids <- bitr(deseq$V1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
entrez_gene_list <- gene_ids$ENTREZID

result_table <- merge(deseq, gene_ids, by.x = "V1", by.y = "SYMBOL", all.x = TRUE)
result_table$V1 <- ifelse(is.na(result_table$ENTREZID), result_table$V1, result_table$ENTREZID)
result_table <- result_table[complete.cases(result_table), ]

out_table <- list()

for (i in statColumns) {
  stats <- setNames(result_table[[i]], result_table$ENTREZID)
  fgseaResKegg <- fgsea(pathways=keggdb, stats=stats, minSize=15, maxSize=500)
  print(fgseaResKegg)
  out_table[[i]] <- fgseaResKegg
}

combined_results <- do.call(rbind, lapply(names(out_table), function(name) {
  res <- out_table[[name]]
  res$Comparison <- name
  return(as.data.frame(res))
}))

combined_results <- combined_results %>%
  mutate(across(where(is.list), ~map_chr(., toString)))

combined_results <- combined_results %>%
  unnest(cols = everything())

write.csv(combined_results, file="fgsea_kegg_edgeR.csv", row.names=FALSE)
