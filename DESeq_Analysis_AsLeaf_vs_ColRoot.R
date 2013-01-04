setwd("C:/Dropbox/HMM_Duplication_Search/Data/RNA-Seq/AsLeaf_vs_ColRoot")
myPath = file.path("C:/Dropbox/HMM_Duplication_Search/Data/RNA-Seq/AsLeaf_vs_ColRoot/mRNA_Expression.tsv")
counts = read.table(myPath, header=TRUE, row.names=1)
metadata = data.frame(
	row.names = colnames(counts),
	condition = c("colroot","sueleaf"),
	libType = c("single-end", "single-end") )
singleSamples = metadata$libType == "single-end"
countTable = counts[,singleSamples]
condition = metadata$condition[singleSamples]
library("DESeq")
cds = newCountDataSet(countTable, condition)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
#head( counts( cds, normalized=TRUE) ) # Prints library-normalized read counts
cds = estimateDispersions(cds, method='blind', sharingMode="fit-only")
#plotDispEsts(cds) # Plot dispersion estimates

# Examine differential expression between A. thaliana (Col) root tissue and A. suecica leaf tissue
res = nbinomTest(cds, "colroot", "sueleaf")
#head(res)
#plotMA(res)
#hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="") # Prints histogram of p-values
#hist(res$padj, breaks=100, col="skyblue", border="slateblue", main="") # Prints histogram of adjusted p-values
#resSig2wk = res[res$padj < 0.1,] # With only one biological sample, all gene-by-gene comparisons will have padj=1. Must use pval instead
resSig = res[res$pval < 0.1,]
resSig = na.omit(resSig)
resSigUp = resSig[resSig$foldChange > 1,]
resSigDown = resSig[resSig$foldChange < 1,]
head(resSig[order(resSig$pval),]) # Most significantly differentially expressed genes
head(resSigUp[order(resSigUp$foldChange, -resSigUp$baseMean ), ] ) # Most strongly down-regulated significant genes
head(resSigDown[order(-resSigDown$foldChange, -resSigDown$baseMean ), ] ) # Most strongly up-regulated significant genes
write.table(resSig[order(resSig$pval), ], file="Col_to_Sue_SigDifExp.tsv", quote=FALSE, sep="\t", row.names=FALSE )
write.table(resSigUp[order(-resSigUp$foldChange, -resSigUp$baseMean ), ], file="Col_to_Sue_SigMostUpReg.tsv", quote=FALSE, sep="\t", row.names=FALSE )
write.table(resSigDown[order(resSigDown$foldChange, -resSigDown$baseMean ), ], file="Col_to_Sue_SigMostDownReg.tsv", quote=FALSE, sep="\t", row.names=FALSE )
