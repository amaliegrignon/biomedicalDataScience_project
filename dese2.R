## Setup
### Bioconductor and CRAN libraries used
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
## Load in data
data <- read.table("count.txt", header=T, row.names=1) 

meta <- read.table("meta.txt", header=T, row.names=1)
### Check classes of the data we just brought in
class(meta)
class(data)

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

dds <- estimateSizeFactors(dds)                              
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
dds
### Plot PCA 
plotPCA(rld, intgroup="sampletype")


# Input is a matrix of log transformed values
rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = sampletype))


### Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames
### Plot heatmap
pheatmap(rld_cor)


dds <- DESeq(dds)
plotDispEsts(dds)



dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

## Run analysis
dds <- DESeq(dds)
#BiocManager::install("apeglm")
#install.packages('ashr')

#### Define contrasts, extract results table, and shrink the log2 fold changes
contrast_oe <- c("sampletype", "SARS", "MOCK")

res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)

res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken, type="ashr")

plotMA(res_tableOE_unshrunken, ylim=c(-2,2))
plotMA(res_tableOE, ylim=c(-2,2))

class(res_tableOE)
mcols(res_tableOE, use.names=T)
res_tableOE %>% data.frame() %>% View()




res_tableKD <- results(dds, contrast=contrast_kd, alpha = 0.05)

res_tableKD <- lfcShrink(dds, contrast=contrast_kd, res=res_tableKD, type="ashr")
summary(res_tableOE)
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sigOE <- res_tableOE_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sigOE


library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
res_tableOE_tb <- res_tableOE_tb %>% 
  mutate(threshold_OE = padj < 0.1 & log2FoldChange <= -0.58)

## Volcano plot
ggplot(res_tableOE_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  
## Create a column to indicate which genes to label
res_tableOE_tb <- res_tableOE_tb %>% arrange(padj) %>% mutate(genelabels = "")

res_tableOE_tb$genelabels[1:10] <- res_tableOE_tb$gene[1:10]

View(res_tableOE_tb)


ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_OE)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 



DEGreport::degPlot(dds = dds, res = res_tableOE_tb, n = 20, xs = "type", group = "condition") # dds object is output from DESeq2

DEGreport::degVolcano(
  data.frame(res[,c("log2FoldChange","padj")]), # table - 2 columns
  plot_text = data.frame(res[1:10,c("log2FoldChange","padj","id")])) # table to add names

# Available in the newer version for R 3.4
DEGreport::degPlotWide(dds = dds, genes = row.names(res)[1:5], group = "condition")


sessionInfo()
