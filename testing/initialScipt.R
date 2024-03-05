library(systemPipeRdata)
genWorkenvir(workflow = "rnaseq")
setwd("rnaseq")
# load packages
library(systemPipeR)
library(DT)
# sets path
targetspath <- "~/shared/slokray/sidd_girke/testing/RNA_Seq_targetsPE.txt - targetsPE.tsv"
# loads tsv file
targets <- read.delim(targetspath, comment.char = "#")
# creates interactive datatable
DT::datatable(targets, options = list(scrollX = TRUE, autoWidth = TRUE))
# define SYSargsList
library(systemPipeR)
sal <- SPRproject()
# build and run workflow
sal <- importWF(sal, file_path = "~/shared/slokray/sidd_girke/testing/systemPipeRNAseq.Rmd")
sal
sal <- runWF(sal)
sal <- renderReport(sal)
rmarkdown::render("systemPipeRNAseq.Rmd", clean = TRUE, output_format = "BiocStyle::html_document")
# required packages
appendStep(sal) <- LineWise(code = {
  library(systemPipeR)
}, step_name = "loadSPR")
# read preprocessing
appendStep(sal) <- SYSargsList(step_name = "preprocessing_", targets = "targetsPE.txt",
                               dir = TRUE, wf_file = "preprocessReads/preprocessReads-pe.cwl",
                               input_file = "preprocessReads/preprocessReads-pe.yml", dir_path = "param/cwl",
                               inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_",
                                             SampleName = "_SampleName_"), dependency = c("loadSPR"))
yamlinput(sal, "preprocessing")$Fct
# fastq report
appendStep(sal) <- LineWise(code = {
  fastq <- getColumn(sal, step = "preprocessing_", "targetsWF",
                     column = 1)
  fqlist <- seeFastq(fastq = fastq, batchsize = 10000, klength = 8)
  pdf("./results/fastqReport.pdf", height = 18, width = 4 *
        length(fqlist))
  seeFastqPlot(fqlist)
  dev.off()
}, step_name = "fastq_report_", dependency = "preprocessing_")
# read mapping with hisat
appendStep(sal) <- SYSargsList(step_name = "hisat2index", dir = FALSE,
                               targets = NULL, wf_file = "hisat2/hisat2-index.cwl", input_file = "hisat2/hisat2-index.yml",
                               dir_path = "param/cwl", dependency = "loadSPR")
appendStep(sal) <- SYSargsList(step_name = "hisat2_mapping_",
                               dir = TRUE, targets = "preprocessing", wf_file = "workflow-hisat2/workflow_hisat2-pe.cwl",
                               input_file = "workflow-hisat2/workflow_hisat2-pe.yml", dir_path = "param/cwl",
                               inputvars = c(preprocessReads_1 = "_FASTQ_PATH1_", preprocessReads_2 = "_FASTQ_PATH2_",
                                             SampleName = "_SampleName_"), rm_targets_col = c("FileName1",
                                                                                              "FileName2"), dependency = c("preprocessing_", "hisat2index"))
cmdlist(sal, step = "hisat2_mapping_", targets = 1)
# read and alignment stats
appendStep(sal) <- LineWise(code = {
  fqpaths <- getColumn(sal, step = "preprocessing_", "targetsWF",
                       column = "FileName1")
  bampaths <- getColumn(sal, step = "hisat2_mapping_", "outfiles",
                        column = "samtools_sort_bam")
  read_statsDF <- alignStats(args = bampaths, fqpaths = fqpaths,
                             pairEnd = TRUE)
  write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE,
              quote = FALSE, sep = "\t")
}, step_name = "align_stats_", dependency = "hisat2_mapping_")
read.table("results/alignStats.xls", header = TRUE)[1:4, ]
#creating symbolic links for viewing
appendStep(sal) <- LineWise(code = {
  bampaths <- getColumn(sal, step = "hisat2_mapping_", "outfiles",
                        column = "samtools_sort_bam")
  bampaths <- setNames(normalizePath(bampaths), names(bampaths))
  symLink2bam(sysargs = bampaths, htmldir = c("~/.html/", "somedir/"),
              urlbase = "https://cluster.hpcc.ucr.edu/~slokr001/",
              urlfile = "./results/IGVurl.txt")
}, step_name = "bam_urls", dependency = "hisat2_mapping_", run_step = "optional")
#create database for gene annotation
appendStep(sal) <- LineWise(code = {
  library(GenomicFeatures)
  txdb <- suppressWarnings(makeTxDbFromGFF(file = "data/tair10.gff",
                                           format = "gff", dataSource = "TAIR", organism = "Arabidopsis thaliana"))
  saveDb(txdb, file = "./data/tair10.sqlite")
}, step_name = "create_db_", dependency = "hisat2_mapping_")
#read counting
appendStep(sal) <- LineWise(code = {
  library(GenomicFeatures)
  library(BiocParallel)
  txdb <- loadDb("./data/tair10.sqlite")
  outpaths <- getColumn(sal, step = "hisat2_mapping_", "outfiles",
                        column = "samtools_sort_bam")
  eByg <- exonsBy(txdb, by = c("gene"))
  bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())
  multicoreParam <- MulticoreParam(workers = 4)
  register(multicoreParam)
  registered()
  counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg,
                                                           x, mode = "Union", ignore.strand = TRUE, inter.feature = FALSE,
                                                           singleEnd = FALSE, BPPARAM = multicoreParam))
  countDFeByg <- sapply(seq(along = counteByg), function(x) assays(counteByg[[x]])$counts)
  rownames(countDFeByg) <- names(rowRanges(counteByg[[1]]))
  colnames(countDFeByg) <- names(bfl)
  rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts = x,
                                                             ranges = eByg))
  write.table(countDFeByg, "results/countDFeByg.xls", col.names = NA,
              quote = FALSE, sep = "\t")
  write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names = NA,
              quote = FALSE, sep = "\t")
  #creating a SummarizedExperiment object
  colData <- data.frame(row.names = SampleName(sal, "hisat2_mapping_"),
                        condition = getColumn(sal, "hisat2_mapping_", position = "targetsWF",
                                              column = "Factor"))
  colData$condition <- factor(colData$condition)
  countDF_se <- SummarizedExperiment::SummarizedExperiment(assays = countDFeByg,
                                                           colData = colData)
  #add results as SummarizedExperiment to the workflow object
  SE(sal, "read_counting") <- countDF_se
}, step_name = "read_counting_", dependency = "create_db")
countDF <- read.delim("results/countDFeByg.xls", row.names = 1,
                      check.names = FALSE)[1:200, ]
DT::datatable(countDF, options = list(scrollX = TRUE, autoWidth = TRUE))
read.delim("results/rpkmDFeByg.xls", row.names = 1, check.names = FALSE)[1:4,1:4]
#samplewise correlation analysis
appendStep(sal) <- LineWise(code = {
  library(DESeq2, quietly = TRUE)
  library(ape, warn.conflicts = FALSE)
  ## Extracting SummarizedExperiment object
  se <- SE(sal, "read_counting")
  dds <- DESeqDataSet(se, design = ~condition)
  d <- cor(assay(rlog(dds)), method = "spearman")
  hc <- hclust(dist(1 - d))
  pdf("results/sample_tree.pdf")
  plot.phylo(as.phylo(hc), type = "p", edge.col = "darkblue", edge.width = 2, show.node.label = TRUE, no.margin = TRUE)
  dev.off()
}, step_name = "sample_tree_", dependency = "read_counting_")
#run edgeR
appendStep(sal) <- LineWise(code = {
  library(edgeR)
  countDF <- read.delim("results/countDFeByg.xls", row.names = 1,
                        check.names = FALSE)
  cmp <- readComp(stepsWF(sal)[["hisat2_mapping_"]], format = "matrix",
                  delim = "-")
  edgeDF <- run_edgeR(countDF = countDF, targets = targetsWF(sal)[["hisat2_mapping_"]],
                      cmp = cmp[[1]], independent = FALSE, mdsplot = "")
}, step_name = "run_edger_", dependency = "read_counting_")
appendStep(sal) <- LineWise(code = {
  library("biomaRt")
  m <- useMart("plants_mart", dataset = "athaliana_eg_gene",
               host = "https://plants.ensembl.org")
  desc <- getBM(attributes = c("tair_locus", "description"),
                mart = m)
  desc <- desc[!duplicated(desc[, 1]), ]
  descv <- as.character(desc[, 2])
  names(descv) <- as.character(desc[, 1])
  edgeDF <- data.frame(edgeDF, Desc = descv[rownames(edgeDF)],
                       check.names = FALSE)
  write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote = FALSE,
              sep = "\t", col.names = NA)
}, step_name = "custom_annot_", dependency = "run_edger_")
appendStep(sal) <- LineWise(code = {
  edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names = 1,
                       check.names = FALSE)
  pdf("results/DEGcounts.pdf")
  DEG_list <- filterDEGs(degDF = edgeDF, filter = c(Fold = 2,
                                                    FDR = 20))
  dev.off()
  write.table(DEG_list$Summary, "./results/DEGcounts.xls",
              quote = FALSE, sep = "\t", row.names = FALSE)
}, step_name = "filter_degs_", dependency = "custom_annot_")
appendStep(sal) <- LineWise(code = {
  vennsetup <- overLapper(DEG_list$Up[6:9], type = "vennsets")
  vennsetdown <- overLapper(DEG_list$Down[6:9], type = "vennsets")
  pdf("results/vennplot.pdf")
  vennPlot(list(vennsetup, vennsetdown), mymain = "", mysub = "",
           colmode = 2, ccol = c("blue", "red"))
  dev.off()
}, step_name = "venn_diagram_", dependency = "filter_degs_")
