library(edgeR)


Counts <- read.table("count.txt", header=TRUE, row.names=1)


meta <- read.table("meta.txt", header=TRUE, row.names=1)


dgList <- DGEList(counts=Counts)


dim(Counts)
head(Counts)


sampleType <- meta$sampletype
sampleType

countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]


dgList <- calcNormFactors(dgList, method="TMM")


plotMDS(dgList, labels=colnames(dgList$counts), col=as.numeric(as.factor(sampleType)))


sampleReplicate <- c("S1", "S1", "S2", "S2")
designMat <- model.matrix(~sampleReplicate + sampleType)


dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)


plotBCV(dgList)


fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef=3)  # 


edgeR_result <- topTags(lrt, n=15000)$table


significant_genes_edgeR <- edgeR_result[edgeR_result$FDR < 0.05, ]
significant_genes_edgeR



















Counts <- read.table("Book2.txt", header=TRUE, row.names=1)

dim(Counts)
head(Counts)
dgList <- DGEList(counts=Counts, genes=rownames(Counts))
dgList
dgList$samples
head(dgList$counts) #Many rows!
head(dgList$genes)
countsPerMillion <- cpm(dgList)
summary(countsPerMillion)

countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]
summary(cpm(dgList))

?calcNormFactors
dgList <- calcNormFactors(dgList, method="TMM")
plotMDS(dgList)

sampleType <- rep("N", ncol(dgList)) # N=normal; T=tumour
sampleType
sampleType[grep("SARS", colnames(dgList))] <- "S"
sampleType
sampleReplicate <- paste("S", rep(1:2, each=2), sep="")
sampleReplicate
designMat <- model.matrix(~sampleReplicate + sampleType)
designMat
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)
plotBCV(dgList)
fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef=3)
edgeR_result <- topTags(lrt,Inf)

edgeR_result
sig_infected <- edgeR_result$table[edgeR_result$table$FDR < 0.05 & abs(edgeR_result$table$logFC) > 0.58, ]
sig_infected
sigOE$gene
common_elements <- intersect(sig_infected$genes, sigOE$gene)
print(common_elements)
edgeR_result$table[edgeR_result$table$genes < 0.1 & abs(edgeR_result$table$logFC) > 0.58, ]

deGenes <- decideTests(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags=deGenes)
abline(h=c(-1, 1), col=2)
deGenes
