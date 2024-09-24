# Install and load VennDiagram package
library(VennDiagram)
library(corrplot)
library(ggplot2)



# Extract gene names (assuming row names or 'gene' column contain gene IDs)
genes1 <- rownames(brca_res_sig)
genes2 <- rownames(ret_res_sig)
genes3 <- rownames(kas_res_sig)
genes4 <- rownames(thp_res_sig)

genes11 <- brca_res_sig$symbol
genes22 <- ret_res_sig$symbol
genes33 <- kas_res_sig$symbol
genes44 <- thp_res_sig$symbol

# Find common DEGs across all four cell lines
common_genes <- Reduce(intersect, list(genes1, genes2, genes3, genes4))


common_genez <- Reduce(intersect, list(genes11, genes22, genes33, genes44))


# Create Venn Diagram for four sets
draw.quad.venn(
  area1 = length(genes1),
  area2 = length(genes2),
  area3 = length(genes3),
  area4 = length(genes4),
  n12 = length(intersect(genes1, genes2)),
  n13 = length(intersect(genes1, genes3)),
  n14 = length(intersect(genes1, genes4)),
  n23 = length(intersect(genes2, genes3)),
  n24 = length(intersect(genes2, genes4)),
  n34 = length(intersect(genes3, genes4)),
  n123 = length(Reduce(intersect, list(genes1, genes2, genes3))),
  n124 = length(Reduce(intersect, list(genes1, genes2, genes4))),
  n134 = length(Reduce(intersect, list(genes1, genes3, genes4))),
  n234 = length(Reduce(intersect, list(genes2, genes3, genes4))),
  n1234 = length(Reduce(intersect, list(genes1, genes2, genes3, genes4))),
  category = c("MCF-7", "Y79", "Kasumi-1", "THP-1"),
  fill = c("skyblue", "pink1", "mediumorchid", "orange"),
  lty = "blank",  # No lines around circles
  cex = 1.5,
  cat.cex = 1.2,
  cat.col = c("darkblue", "darkred", "darkorchid", "darkorange")
)


# correlation studies
# Subset each result to only keep the common genes
brca_res_common <- brca_res_sig[rownames(brca_res_sig) %in% common_genes, ]
brcdf <- as.data.frame(brca_res_common)

ret_res_common <- ret_res_sig[rownames(ret_res_sig) %in% common_genes, ]
rrcdf <- as.data.frame(brca_res_common)

kas_res_common <- kas_res_sig[rownames(kas_res_sig) %in% common_genes, ]
krcdf <- as.data.frame(brca_res_common)

thp_res_common <- thp_res_sig[rownames(thp_res_sig) %in% common_genes, ]
trcdf <- as.data.frame(brca_res_common)

# Assuming rownames contain the gene IDs, merge the data
merged_results <- merge(brca_res_common[, c("log2FoldChange", "padj")], 
                        ret_res_common[, c("log2FoldChange", "padj")], 
                        by = "row.names", suffixes = c("_MCF7", "_Y79"))

rownames(merged_results) <- merged_results$Row.names
merged_results$Row.names <- NULL

# Merge with the third cell line results
merged_results <- merge(merged_results, 
                        kas_res_common[, c("log2FoldChange", "padj")], 
                        by = "row.names", suffixes = c("_Kasumi1", "_THP1"))

rownames(merged_results) <- merged_results$Row.names
merged_results$Row.names <- NULL

# Finally, merge with the fourth cell line results
merged_results <- merge(merged_results, 
                        thp_res_common[, c("log2FoldChange", "padj")], 
                        by = "row.names", suffixes = c("_Kasumi1", "_THP1"))

rownames(merged_results) <- merged_results$Row.names
merged_results$Row.names <- NULL

# View the merged results
head(merged_results)

# Save to a CSV file
write.csv(merged_results, "common_genes_logFC_padj.csv", row.names = TRUE)


# Subset only the log2FoldChange columns
lfc_data <- merged_results[, grep("log2FoldChange", colnames(merged_results))]
colnames(lfc_data) <- c("MCF7", "Y79", "Kasumi1", "THP1")

# View the first few rows
head(lfc_data)

# Calculate Pearson correlation matrix
corr_matrix <- cor(lfc_data, method = "pearson")

# Calculate Spearman correlation matrix (optional)
spearman_corr_matrix <- cor(lfc_data, method = "spearman")


# Visualize the Pearson correlation matrix
corrplot(corr_matrix, method = "color", addCoef.col = "black", tl.col = "black", 
         title = "Pearson Correlation Matrix", mar = c(0,0,2,0))

# Optional: Visualize the Spearman correlation matrix
corrplot(spearman_corr_matrix, method = "color", addCoef.col = "black", tl.col = "black", 
         title = "Spearman Correlation Matrix", mar = c(0,0,2,0))


# add gene symbols
lfc_data$symbol = mapIds(org.Hs.eg.db, 
                           keys = row.names(lfc_data), 
                           keytype = "ENSEMBL", 
                           column = "SYMBOL")

# Example scatter plot between MCF7 and Y79
ggplot(lfc_data, aes(x = MCF7, y = Y79)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  ggtitle("MCF7 vs Y79 Log2 Fold Change Correlation") +
  xlab("MCF7 Log2FoldChange") +
  ylab("Y79 Log2FoldChange")


# Example scatter plot between Kasumi-1 and THP-1
ggplot(lfc_data, aes(x = Kasumi1, y = THP1)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  ggtitle("Kasumi-1 vs THP-1 Log2 Fold Change Correlation") +
  xlab("Kasumi-1 Log2FoldChange") +
  ylab("THP-1 Log2FoldChange")


# Example scatter plot between Kasumi-1 and MCF7
ggplot(lfc_data, aes(x = Kasumi1, y = MCF7)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  ggtitle("Kasumi-1 vs MCF7 Log2 Fold Change Correlation") +
  xlab("Kasumi-1 Log2FoldChange") +
  ylab("MCF7 Log2FoldChange")

# Example scatter plot between THP1 and MCF7
ggplot(lfc_data, aes(x = THP1, y = MCF7)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  ggtitle("THP1 vs MCF7 Log2 Fold Change Correlation") +
  xlab("THP1 Log2FoldChange") +
  ylab("MCF7 Log2FoldChange")
