library(edgeR)

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
