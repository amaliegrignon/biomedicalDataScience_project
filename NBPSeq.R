library(NBPSeq)
library(edgeR)

counts_data <- read.table("Book2.txt", header=TRUE, row.names=1, sep="\t")


meta_data <- data.frame(
  sampletype = c("SARS", "MOCK", "SARS", "MOCK"),
  row.names = c("S_2_SARS_5dpi_S70003", "S_2_mock_5dpi_S70002", "S_3_SARS_5dpi_S69996", "S_3_mock_5dpi_S69997")
)
counts_data <- as.matrix(counts_data)


dge_list <- DGEList(counts = counts_data)
dge_list

cpm_data <- cpm(dge_list)


keep <- rowSums(cpm_data > 1) >= 2
filtered_counts_data <- counts_data[keep, ]




condition <- factor(meta_data$sampletype)

dge_list <- DGEList(counts = filtered_counts_data, group = condition)


dge_list <- calcNormFactors(dge_list, method = "TMM")
counts_data = as.matrix(counts_data)
#dge_list <- estimate.norm.factors(counts_data)
#norm_factors = estimate.norm.factors(counts_data);
norm_factors <- as.vector(dge_list$samples$norm.factors)
filtered_counts_data

is.matrix(filtered_counts_data)
set.seed(1)
nbp_test <- nbp.test(filtered_counts_data, grp.ids = condition,
                     grp1 = "SARS", grp2 = "MOCK", norm.factors = norm_factors,model.disp="NBP")

alpha = 0.05;
sig.res = nbp_test$q.values < alpha;
table(sig.res)


p_values <- nbp_test$p.values
adj_p_values <- nbp_test$q.values
lfc <- nbp_test$log.fc

results <- data.frame(
  gene_id = rownames(filtered_counts_data),
  adj_p_values = adj_p_values,
  lfcs = lfc
)
results <- na.omit(results)

significant_genes_NBP <- results[results$adj_p_values < 0.05, ]
print(significant_genes_NBP)

#sigOE$gene
#gene_ids_edgeR <- rownames(significant_genes_edgeR)
#significant_genes_NBP$gene_id

#common_elements <- intersect(gene_ids_edgeR, sigOE$gene)
#print(common_elements)
#common_elements_2 <- intersect(common_elements, significant_genes_NBP$gene_id)
#common_elements_2

