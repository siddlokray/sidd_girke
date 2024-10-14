#######################################################################################################################################
#### RNA-Seq analysis #################################################################################################################
#######################################################################################################################################
library(systemPipeR);library(ShortRead);library(data.table)
homedir <- "/rhome/slokr001/shared/slokray"
setwd(homedir)
# fastq-dump --split-files --gzip --outdir /rhome/hmang007/bigdata/girke_lab_analysis/data SRR5927754


#### https://pubmed.ncbi.nlm.nih.gov/31353263/
# GSE131901 in GEO website

#######################
#### Load function ####
#######################
parallelFASTqDownloadR <- function(Files, Outpath = path, regdir = getwd(), partition="girkelab", WallTime=172800, Memory=80000, Ntask = 1, Ncpus = 1, DeleteRegistry = TRUE){ # ,batch,intel, girkelab
#### Load function ####
parallelFASTqDownload <- function(x, files = Files, outpath = Outpath){
  library('RenvModule'); module('load','slurm/21.08.5'); module('load','sratoolkit/3.0.0')#; module('load','openmpi/4.1.2_slurm-21.08.5')
  system(paste("fastq-dump --split-files --gzip --outdir", outpath, files[x]))
  } #### download files ####
#### Submit jobs from R to cluster
reg <- makeRegistry(file.dir=regdir, conf.file=".batchtools.conf.R")
Njobs <- 1:length(SRRnames) # Define number of jobs
ids <- batchMap(fun=parallelFASTqDownload, x=Njobs,
        more.args = list(files = SRRnames, outpath = Outpath))
FQdownload <- submitJobs(ids, reg=reg, resources=list(partition=partition, walltime=WallTime, ntasks=Ntask, ncpus=Ncpus, memory=Memory)) #; getStatus()
waitForJobs()
if(DeleteRegistry == TRUE){
clearRegistry(); removeRegistry(wait=0, reg=reg) # Clear registry in R session # Delete registry directory
} }
#############################
#### Download fastq files ###
#############################
dt <- fread(file.path(homedir, "data/SraRunInfo.csv"))
SRRnames <- dt$Run
SRRnames
regdir <- file.path(homedir,"myregdir2")
outpath <- file.path(homedir,"results")
setwd(file.path(homedir,"myregdir"))
library('RenvModule'); module('load','slurm/21.08.5'); module('load','openmpi/4.1.2_slurm-21.08.5'); module('load','sratoolkit/3.0.0'); library(batchtools)
parallelFASTqDownloadR(Files = SRRnames, Outpath = outpath, regdir = regdir, partition="girkelab", WallTime=6000, Memory=6000, Ntask = 1, Ncpus = 1, DeleteRegistry = FALSE) # ,batch,intel
getStatus()
getOption('timeout')
options(timeout=2500)

system(paste("fastq-dump --split-files --gzip --outdir", file.path(outpath, "SRR9119201")))

system("vdb-config --interactive")
##############################################################
#### Download GFF/GTF and genome (fasta) file of interest ####
##############################################################
URL <- "https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz"
download.file(URL, destfile = "/rhome/slokr001/shared/slokray/data/Mus_musculus.GRCm39.110.gtf.gz")

URLGEN <- "https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
download.file(URLGEN, destfile = "/rhome/slokr001/shared/slokray/data/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz", timeout = 300)
      
#######################################################
#### Set up targets.txt file and .param file ##########
#######################################################
parampath <- "/rhome/slokr001/shared/slokray/param/hisat2PE.param" ##system.file("extdata", "hisathuman.param", package="systemPipeR")
read.delim(parampath, comment.char = "#")

targets <- read.delim("/rhome/slokr001/shared/slokray/sidd_girke/SraRunInfo - pairedTargets-2.txt", comment.char ="#")
targets#[1:3,]

#### Read .param file and targets.txt file into args object ####
################################################################
args <- systemArgs(sysma=parampath, mytargets="/rhome/slokr001/shared/slokray/sidd_girke/SraRunInfo - pairedTargets-2.txt")
args
modules(args)
cores(args)
names(args)
sysargs(args)[1]
outpaths(args)[1]
normalizePath(results(args))

############################################
#### Generate a quality report #############
############################################
fqlist <- seeFastq(fastq=infile1(args), batchsize=10000, klength=8)
pdf("./fastqReport_human.pdf", height=18, width=4*length(fqlist))
seeFastqPlot(fqlist); dev.off()

##################################
#### 1 Index reference genome ####
##################################
moduleload(modules(args)) # Skip if module system is not available
system("hisat2-build ./data/Mus_musculus.GRCm39.dna.primary_assembly.fa ./data/hisat_musculus") # the .fasta file should contain the genome of interest to align reads to.

########################################################
#### Align all FASTq files with hisat2 on a cluster ####
########################################################
moduleload(modules(args)) # Skip if module system is not available
library(batchtools)
setwd(file.path(homedir,"gl_analysis"))
resources <- list(walltime = 3024, ntasks = 78, ncpus = 8, memory = 70000, partition="girkelab") # , ntasks = 10, ncpus = 5
reg <- clusterRun(args, FUN = runCommandline, more.args = list(args = args, make_bam = TRUE, dir = FALSE),
                  conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl",
                  Njobs = 1, runid = "01", resourceList = resources)
getStatus(reg = reg)
waitForJobs(reg = reg)

########################################
#### Generate an allignment summary ####
########################################
(read_statsDF <- alignStats(args=args))#[1:8,]
write.table(read_statsDF, "./humanalignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")

###########################################################################
#### Counting reads per Feature by Storing Annotations in TranscriptDb ####
###########################################################################
library(GenomicFeatures)
mouse_txdb <- makeTxDbFromGFF(file="/rhome/slokr001/shared/slokray/data/Mus_musculus.GRCm39.110.gtf",
                     format = "gtf",
                     dataSource= "ensembl",
                     organism ="Mus musculus")
saveDb(mouse_txdb, file="/rhome/slokr001/shared/slokray/data/mouseGTF.sqlite")
mouse_txdb <- loadDb("/rhome/slokr001/shared/slokray/data/mouseGTF.sqlite")
mouse_eByg <- exonsBy(mouse_txdb, by="gene")

##############################################
#### Read Counting with summarizeOverlaps ####
##############################################
homedir <- "/rhome/slokr001/shared/slokray"
setwd(homedir)
library("GenomicFeatures"); library(BiocParallel)
mouse_txdb <- loadDb("/rhome/slokr001/shared/slokray/data/mouseGTF.sqlite")
mouse_eByg <- exonsBy(mouse_txdb, by=c("gene"))
mouse_bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
mouse_multicoreParam <- MulticoreParam(workers=7); register(mouse_multicoreParam)
mouse_counteByg <- bplapply(mouse_bfl, function(x) summarizeOverlaps(mouse_eByg, x, mode="Union", ignore.strand=FALSE, inter.feature=TRUE, singleEnd=FALSE, fragments=TRUE)) # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'
saemouse_countDFeByg <- sapply(seq(along=mouse_counteByg), function(x) assays(mouse_counteByg[[x]])$counts)
rownames(mouse_countDFeByg) <- names(rowRanges(mouse_counteByg[[1]])); colnames(mouse_countDFeByg) <- names(mouse_bfl)
mouse_rpkmDFeByg <- apply(mouse_countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=mouse_eByg))
write.table(mouse_countDFeByg, "./mousecountDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(mouse_rpkmDFeByg, "./mouserpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")

################################
### DEG analysis with edgeR ####
################################
#countDF <- read.table("./humancountDFeByg.xls")
mouse_countDF <- read.table("./mousefpkmDFeByg.xls")
mouse_cmp <- readComp(args, format="matrix", delim="-")
mouse_edgeDF <- run_edgeR(countDF=mouse_countDF, targets=targetsin(args), cmp=mouse_cmp[[1]], independent=FALSE, mdsplot="")
write.table(mouse_edgeDF, "./mouseedgeRcomp.xls", quote=FALSE, sep="\t", col.names = NA)

################################
### DEG analysis with DESeq2 ###
################################
mouse_countDF <- read.table("mousecountDFeByg.xls")
mouse_cmp <- readComp(args, format="matrix", delim="-")
mouse_degseqDF <- run_DESeq2(countDF=mouse_countDF, targets=targetsin(args), cmp=mouse_cmp[[1]], independent=FALSE)
write.table(mouse_degseqDF, "./mouseDESeq2comp.xls", quote=FALSE, sep="\t", col.names = NA)
mouse_degseqDF <- read.table("./mouseDESeq2comp.xls")
human_DEG_list2 <- filterDEGs(degDF=mouse_degseqDF, filter=c(Fold=2, FDR=10))

#########################################
#### Correlation analysis of samples ####
#########################################
library(DESeq2, warn.conflicts=FALSE, quietly=TRUE); library(ape, warn.conflicts=FALSE)
countDFpath <- "./mousecountDFeByg.xls"
countDF <- as.matrix(read.table(countDFpath))
colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
dds_rlog <- rlog(dds)
d <- cor(assay(dds_rlog), method="spearman")
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)

pdf("./phylotree.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
dev.off()

#####################################################################
#### Correlation analysis with RPKM normalized expression values ####
#####################################################################
rpkmDFeBygpath <- "./mousefpkmDFeByg.xls"  # system.file("extdata", "rpkmDFeByg.xls", package="systemPipeR")
rpkmDFeByg <- read.table(rpkmDFeBygpath, check.names=FALSE)
rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
d <- cor(rpkmDFeByg, method="spearman")
hc <- hclust(as.dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)

pdf("./rpkmphylotree.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
dev.off()

# boxplots to compare mean distributions
mouse_countDF <- read.table("mousecountDFeByg.xls")
meansCounts <- colMeans(countDF)
rlog <- rlog(countDF)
meansLog <- colMeans(rlog)
meansScale <- apply(countDF, 2, scale)
ms <- as.matrix(meansScale)
meansScale <- colMeans(ms)
meansCountsSTD <- colSds(countDF)
meansLogSTD <- colSds(rlog)
meansScaleSTD <- colSds(ms)
par(mfrow = c(2, 1))
boxplot(meansCounts,meansCountsSTD, horizontal = TRUE, main = "raw counts", xlab = "Mean Count", col = 'skyblue', border = 'brown')
boxplot(meansCountsSTD, horizontal = TRUE, main = "raw counts", xlab = "Mean Count", col = 'skyblue', border = 'brown')
boxplot(meansLog, horizontal = TRUE, main = "rlog() applied", xlab = "Mean Count", col = 'orange', border = 'brown')
boxplot(meansLogSTD, horizontal = TRUE, main = "rlog() applied", xlab = "Mean Count", col = 'orange', border = 'brown')
boxplot(meansScale, horizontal = TRUE, main = "scale() applied", xlab = "Mean Count", col = 'pink', border = 'brown')
boxplot(meansScaleSTD, horizontal = TRUE, main = "scale() applied", xlab = "Mean Count", col = 'pink', border = 'brown')

# load edgeR results
mouse_edgeDF <- read.table("./mouseedgeRcomp.xls")
# plot DEG Counts
DEG_list <- filterDEGs(degDF = mouse_edgeDF, filter = c(Fold = 2, FDR = 20))
# venn diagram setup
vennsetup <- overLapper(DEG_list$Up[c(1,2,4,6,7)], type = "vennsets")
vennsetdown <- overLapper(DEG_list$Down[c(1,2,4,6,7)], type = "vennsets")
vennPlot(list(vennsetup, vennsetdown), mymain = "", mysub = "", colmode = 2, ccol = c("blue", "red"))
# heat map
library(pheatmap)
mat <- as.matrix(mouse_edgeDF)
mat <- mat[rowSums(is.na(mat)) != ncol(mat),]
mat <- mat[is.na(mat)] = 0
pheatmap(mat)

library(ape)
d <- cor(mat, method='spearman')
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type = "p", edge.col = "blue")






library(factoextra)
library(corrplot)
#load counts
counts <- read.table('mousefpkmDFeByg.xls')
counts <- data.frame(counts)
#normalize data
counts <- as.data.frame(sapply(counts, function(x) scale(x)))
#run pca analysis
counts_pca <- prcomp(rpkmDFeByg)
summary(counts_pca)
eigval <- get_eigenvalue(counts_pca)
#plot eigenvalues
fviz_eig(counts_pca, col.var="blue")
## ???
var <- get_pca_var(counts_pca)
corrplot(var$cos2, is.corr=FALSE)
