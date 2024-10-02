library(data.table);library(biomaRt)
homedir <- "/rhome/slokr001/shared/slokray"
setwd(homedir)

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

deseq2_data = fread('mouseedgeRcomp.xls')
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

fwrite(deseq2_data, 'mouseedgeRcomp_updated.xls')

write_csv(deseq2_data, 'mouseedgeRcomp_updated.csv')

deseq2_data = fread('mouseDESeq2comp.xls')
deseq2_data_updated = fread('mouseDESeq2comp_updated.xls')

##### FEA

library(GOstats); library(GO.db)
library(DESeq2); library(clusterProfiler); library(org.Mm.eg.db); library(GO.db); library(pathview); library(readr);library(fgsea)
library(tidyverse)



deseq <- read_csv('mouseedgeRcomp_updated.csv', show_col_types = FALSE)
edgeR <- read_csv('mouseedgeRcomp.xls', show_col_types = FALSE)
# deseq <- read_tsv('mouseDESeq2comp.tsv', show_col_types = FALSE)
deseq <- deseq[complete.cases(deseq), ]

pvalColumns <- c("Control.F-Acarbose.F_PValue", "Control.M-Acarbose.M_PValue", "Control.F-CR.F_PValue", 
                 "Control.M-CR.M_PValue", "GHRKOControl.M-GHRKO.M_PValue", "Control.F-Rapamycin.F_PValue", 
                 "Control.M-Rapamycin.M_PValue", "SnellDwarfMiceControl.M-SnellDwarfMice.M_PValue", 
                 "Control.F-Protandim.F_PValue", "Control.M-Protandim.M_PValue", "Control.F-alphaestradiol.F_PValue", 
                 "Control.M-alphaestradiol.M_PValue", "MRControl.M-MethionineRestriction.M_PValue")
statColumns <- c("Control.F-Acarbose.F_logCPM", "Control.M-Acarbose.M_logCPM", "Control.F-CR.F_logCPM", 
                 "Control.M-CR.M_logCPM", "GHRKOControl.M-GHRKO.M_logCPM", "Control.F-Rapamycin.F_logCPM", 
                 "Control.M-Rapamycin.M_logCPM", "SnellDwarfMiceControl.M-SnellDwarfMice.M_logCPM", 
                 "Control.F-Protandim.F_logCPM", "Control.M-Protandim.M_logCPM", "Control.F-alphaestradiol.F_logCPM", 
                 "Control.M-alphaestradiol.M_logCPM", "MRControl.M-MethionineRestriction.M_logCPM")
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
  fgseaResKegg <- fgsea(pathways=godb, stats=stats, minSize=15, maxSize=500)
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

write.csv(combined_results, file="fgsea_go_results_edgeR_BH.csv", row.names=FALSE)











combined_results <- data.frame()

for (i in statColumns) {
  stats <- setNames(result_table[[i]], result_table$ENTREZID)
  stats <- stats[order(stats, decreasing = TRUE)]
  
  go_enrich <- enrichGO(gene = names(stats),
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENTREZID",
                        ont = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
  
  print(go_enrich)
  go_enrich_df <- as.data.frame(go_enrich)
  go_enrich_df$Comparison <- i
  
  combined_results <- rbind(combined_results, go_enrich_df)
}

write.csv(combined_results, file="fgsea_go_results.csv", row.names=FALSE)
