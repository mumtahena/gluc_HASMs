library(Rsubread)
library(limma)
library(edgeR)
library(tools)
options(digits=2)

referenceGenomeFastaFilePath = commandArgs()[7]
inFilePath1 = commandArgs()[8]
inFilePath2 = commandArgs()[9] # NULL for single-end analyses or when a BAM file has been specified
gtfFilePath = commandArgs()[10]
tempFilePrefix = commandArgs()[11]
outFpkmFilePath = commandArgs()[12]
outTpmFilePath = commandArgs()[13]
outCountsFilePath = commandArgs()[14]
outStatsFilePath = commandArgs()[15]

print(referenceGenomeFastaFilePath) 

memory = 4000
nthreads = 1

input_format = "gzFASTQ"
if (file_ext(inFilePath1) == "bam")
  input_format = "BAM"
if (file_ext(inFilePath1) %in% c("fastq", "fq"))
  input_format = "FASTQ"

outBamFilePath = paste(tempFilePrefix, "bam", sep=".")
print(outBamFilePath)
referenceGenomeIndexFilePrefix = paste(referenceGenomeFastaFilePath, "__reference_index", sep="")

if (!file.exists(paste(referenceGenomeIndexFilePrefix, ".reads", sep="")))
  buildindex(basename=referenceGenomeIndexFilePrefix, reference=referenceGenomeFastaFilePath, memory=memory)

if (inFilePath2 == "NULL")
  inFilePath2 = NULL

if (!file.exists(outBamFilePath))
  align(index=referenceGenomeIndexFilePrefix, readfile1=inFilePath1, readfile2=inFilePath2, type="rna", output_file=outBamFilePath, nthreads=nthreads, input_format=input_format, unique=TRUE, indels=5)
print(outBamFilePath)

#align(index=referenceGenomeIndexFilePrefix, readfile1=inFilePath1, readfile2=inFilePath2, output_file=outBamFilePath, nthreads=nthreads, input_format=input_format,  unique=TRUE, indels=5)
fCountsList = featureCounts(outBamFilePath, annot.inbuilt="hg38",annot.ext=NULL, isGTFAnnotationFile=FALSE, nthreads=nthreads, isPairedEnd=!is.null(inFilePath2), verbose=TRUE)
print(fCountsList)
dgeList = DGEList(counts=fCountsList$counts, genes=fCountsList$annotation)
fpkm = rpkm(dgeList, dgeList$genes$Length)
tpm = exp(log(fpkm) - log(sum(fpkm)) + log(1e6))

write.table(fCountsList$stat, outStatsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

featureCounts = cbind(fCountsList$annotation[,1], fCountsList$counts)
write.table(featureCounts, outCountsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
outFpkmLogFilePath=paste(outFpkmFilePath,"log",sep="_")
outTpmLogFilePath= paste(outTpmFilePath,"log", sep="_")
write.table(cbind(fCountsList$annotation[,1], fpkm), outFpkmFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(cbind(fCountsList$annotation[,1], log2(fpkm + 1)), outFpkmLogFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(cbind(fCountsList$annotation[,1], tpm), outTpmFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(cbind(fCountsList$annotation[,1], log2(tpm + 1)), outTpmLogFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
print(outFpkmLogFilePath)
#unlink(outBamFilePath)
#unlink(paste(outBamFilePath, ".indel", sep=""))
