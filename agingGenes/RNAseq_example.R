#######################################################################################################################################
#### RNA-Seq analysis #################################################################################################################
#######################################################################################################################################
library(systemPipeR);library(ShortRead);library(data.table)
homedir <- "/rhome/slokr001/shared/slokray/sidd_girke/agingGenes/"
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
dt <- fread(file.path(homedir, "", "SraRunInfo.csv file"))
SRRnames <- dt$Run
regdir <- file.path(homedir,"myregdir")
outpath <- file.path(homedir,"data")
setwd(file.path(homedir,"gl_analysis"))
library('RenvModule'); module('load','slurm/21.08.5'); module('load','openmpi/4.1.2_slurm-21.08.5'); module('load','sratoolkit/3.0.0'); library(batchtools)
parallelFASTqDownloadR(Files = SRRnames, Outpath = outpath, regdir = regdir, partition="girkelab", WallTime=6000, Memory=6000, Ntask = 1, Ncpus = 1, DeleteRegistry = FALSE) # ,batch,intel
getStatus()
getOption('timeout')
options(timeout=2500)

##############################################################
#### Download GFF/GTF and genome (fasta) file of interest ####
##############################################################
URL <- "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
download.file(URL, destfile = "/rhome/slokr001/shared/slokray/sidd_girke/agingGenes/data/humanGTF.gtf.gz")

URLGEN <- "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
download.file(URLGEN, destfile = "/rhome/slokr001/shared/slokray/sidd_girke/agingGenes/data/humanGTFprim.gtf.gz", timeout = 120)

#######################################################
#### Set up targets.txt file and .param file ##########
#######################################################
parampath <- "/rhome/slokr001/shared/slokray/sidd_girke/agingGenes/gl_analysis/hisathuman.param" ##system.file("extdata", "hisathuman.param", package="systemPipeR")
read.delim(parampath, comment.char = "#")

targets <- read.delim("/rhome/hmang007/bigdata/girke_lab_analysis/gl_analysis/targethuman.txt", comment.char ="#")
targets#[1:3,]

#### Read .param file and targets.txt file into args object ####
################################################################
args <- systemArgs(sysma=parampath, mytargets="/rhome/hmang007/bigdata/girke_lab_analysis/gl_analysis/targethuman.txt")
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
system("hisat2-build ./data/humangenome.fasta ./data/humangenome.fasta") # the .fasta file should contain the genome of interest to align reads to.

########################################################
#### Align all FASTq files with hisat2 on a cluster ####
########################################################
moduleload(modules(args)) # Skip if module system is not available
library(batchtools)
setwd(file.path(homedir,"gl_analysis"))
resources <- list(walltime = 3024, ntasks = 1, ncpus = 1, memory = 6000, partition="girkelab") # , ntasks = 10, ncpus = 5
reg <- clusterRun(args, FUN = runCommandline, more.args = list(args = args, make_bam = TRUE, dir = FALSE),
                  conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl",
                  Njobs = 10, runid = "01", resourceList = resources)
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
human_txdb <- makeTxDbFromGFF(file="/rhome/hmang007/bigdata/girke_lab_analysis/data/humanGTF.gtf",
                     format = "gtf",
                     dataSource= "ensembl",
                     organism ="Homo sapiens")
saveDb(human_txdb, file="/rhome/hmang007/bigdata/girke_lab_analysis/data/humanGTF.sqlite")
human_txdb <- loadDb("/rhome/hmang007/bigdata/girke_lab_analysis/data/humanGTF.sqlite")
human_eByg <- exonsBy(human_txdb, by="gene")

##############################################
#### Read Counting with summarizeOverlaps ####
##############################################
library("GenomicFeatures"); library(BiocParallel)
human_txdb <- loadDb("/rhome/hmang007/bigdata/girke_lab_analysis/data/humanGTF.sqlite")
human_eByg <- exonsBy(human_txdb, by=c("gene"))
human_bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
human_multicoreParam <- MulticoreParam(workers=7); register(human_multicoreParam)
human_counteByg <- bplapply(human_bfl, function(x) summarizeOverlaps(human_eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=TRUE, singleEnd=TRUE)) # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'
human_countDFeByg <- sapply(seq(along=human_counteByg), function(x) assays(human_counteByg[[x]])$counts)
rownames(human_countDFeByg) <- names(rowRanges(human_counteByg[[1]])); colnames(human_countDFeByg) <- names(human_bfl)
human_rpkmDFeByg <- apply(human_countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=human_eByg))
write.table(human_countDFeByg, "./humancountDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(human_rpkmDFeByg, "./humanrpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")

################################
### DEG analysis with edgeR ####
################################
#countDF <- read.table("./humancountDFeByg.xls")
human_countDF <- read.table("./humanrpkmDFeByg.xls")
human_cmp <- readComp(args, format="matrix", delim="-")
human_edgeDF <- run_edgeR(countDF=human_countDF, targets=targetsin(args), cmp=human_cmp[[1]], independent=FALSE, mdsplot="")
write.table(human_edgeDF, "./humanedgeRcomp.xls", quote=FALSE, sep="\t", col.names = NA)

################################
### DEG analysis with DESeq2 ###
################################
human_countDF <- read.table("./humancountDFeByg.xls")
human_cmp <- readComp(args, format="matrix", delim="-")
human_degseqDF <- run_DESeq2(countDF=human_countDF, targets=targetsin(args), cmp=human_cmp[[1]], independent=FALSE)
write.table(human_degseqDF, "./humanDESeq2comp.xls", quote=FALSE, sep="\t", col.names = NA)
human_DEG_list2 <- filterDEGs(degDF=human_degseqDF, filter=c(Fold=2, FDR=10))










#########################################
#### Correlation analysis of samples ####
#########################################
library(DESeq2, warn.conflicts=FALSE, quietly=TRUE); library(ape, warn.conflicts=FALSE)
countDFpath <- "./results/3 countDFeByg.xls"    # system.file("extdata", "countDFeByg.xls", package="systemPipeR")
countDF <- as.matrix(read.table(countDFpath))
colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
d <- cor(assay(rlog(dds)), method="spearman")
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)

pdf("./results/5 phylotree.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
dev.off()

#####################################################################
#### Correlation analysis with RPKM normalized expression values ####
#####################################################################
rpkmDFeBygpath <- "./results/4 rpkmDFeByg.xls"  # system.file("extdata", "rpkmDFeByg.xls", package="systemPipeR")
rpkmDFeByg <- read.table(rpkmDFeBygpath, check.names=FALSE)
rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
d <- cor(rpkmDFeByg, method="spearman")
hc <- hclust(as.dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)

pdf("./results/6 rpkmphylotree.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
dev.off()
