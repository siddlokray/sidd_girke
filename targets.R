targets <- read.delim("/rhome/slokr001/shared/slokray/SraRunInfo - pairedTargets.txt", comment.char ="#")
sampleName <- targets$SampleName
sampleNameUnique <- unique(sampleName)
long <- c(sampleNameUnique, sampleName)
long <- make.unique(as.character(long), sep = ".")
sampleName <- long[19:96]
targets$SampleName <- sampleName
write.table(targets, '/rhome/slokr001/shared/slokray/SraRunInfo - pairedTargets.xls')


sampleName
