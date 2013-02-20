setwd("C:/Dropbox/SuecicaDupSearch/Data/RNA-Seq/AsLeaf_vs_ColLeaf")
myPath = file.path("C:/Dropbox/SuecicaDupSearch/Data/RNA-Seq/AsLeaf_vs_ColLeaf/mRNA_Expression_Technicals_Combined.tsv")
counts = read.table(myPath, header=TRUE, row.names=1)
metadata = data.frame(
	row.names = colnames(counts),
	condition = c("col_leaf_rep", "col_leaf_rep", "sue_leaf_rep", "sue_leaf_rep"),
	libType = c("single-end", "single-end", "single-end", "single-end") )
singleSamples = metadata$libType == "single-end"
countTable = counts[,singleSamples]
condition = metadata$condition[singleSamples]
library("DESeq")
cds = newCountDataSet(countTable, condition)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
#head( counts( cds, normalized=TRUE) ) # Prints library-normalized read counts
cds = estimateDispersions(cds, method='blind', sharingMode="fit-only")
plotDispEsts(cds) # Plot dispersion estimates

# Examine differential expression between A. thaliana (Col) and A. suecica leaf tissue for each species' two biological replicates
res = nbinomTest(cds, "col_leaf_rep", "sue_leaf_rep")
#head(res)
#plotMA(res)
#hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="") # Prints histogram of p-values
#hist(res$padj, breaks=100, col="skyblue", border="slateblue", main="") # Prints histogram of adjusted p-values
resSig = res[res$padj < 0.1,]
resSig = na.omit(resSig)
resSigUp = resSig[resSig$foldChange > 1,]
resSigDown = resSig[resSig$foldChange < 1,]
head(resSig[order(resSig$pval),]) # Most significantly differentially expressed genes
head(resSigUp[order(-resSigUp$foldChange, -resSigUp$baseMean ), ] ) # Most strongly up-regulated significant genes
head(resSigDown[order(resSigDown$foldChange, -resSigDown$baseMean ), ] ) # Most strongly down-regulated significant genes
write.table(resSig[order(resSig$pval), ], file="Col_to_Sue_SigDifExp.tsv", quote=FALSE, sep="\t", row.names=FALSE )
write.table(resSigUp[order(-resSigUp$foldChange, -resSigUp$baseMean ), ], file="Col_to_Sue_SigMostUpReg.tsv", quote=FALSE, sep="\t", row.names=FALSE )
write.table(resSigDown[order(resSigDown$foldChange, -resSigDown$baseMean ), ], file="Col_to_Sue_SigMostDownReg.tsv", quote=FALSE, sep="\t", row.names=FALSE )

# Examine differential expression between A. thaliana (Col) leaf tissue for the two biological replicates with no p-value cutoff applied
rm(resSig,resSigUp,resSigDown)
resSig = res[res$padj < 1,]
resSig = na.omit(resSig)
resSigUp = resSig[resSig$foldChange > 1,]
resSigDown = resSig[resSig$foldChange < 1,]
head(resSig[order(resSig$pval),]) # Most significantly differentially expressed genes
head(resSigUp[order(-resSigUp$foldChange, -resSigUp$baseMean ), ] ) # Most strongly up-regulated significant genes
head(resSigDown[order(resSigDown$foldChange, -resSigDown$baseMean ), ] ) # Most strongly down-regulated significant genes
write.table(resSig[order(resSig$pval), ], file="Col_to_Sue_SigDifExp_NOCUTOFF.tsv", quote=FALSE, sep="\t", row.names=FALSE )
write.table(resSigUp[order(-resSigUp$foldChange, -resSigUp$baseMean ), ], file="Col_to_Sue_SigMostUpReg_NOCUTOFF.tsv", quote=FALSE, sep="\t", row.names=FALSE )
write.table(resSigDown[order(resSigDown$foldChange, -resSigDown$baseMean ), ], file="Col_to_Sue_SigMostDownReg_NOCUTOFF.tsv", quote=FALSE, sep="\t", row.names=FALSE )
