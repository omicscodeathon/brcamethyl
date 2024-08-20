#!/usr/bin/env R

#loading packages required for the analysis
library(DESeq2)
library(ggpubr)
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
library(KEGG.db)
library(KEGGREST)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)

#setting working directory
setwd("/Users/jkalami/Desktop/ich/rnabrca/counts/")

countdata <- read.table("feature_Counts.txt", header = TRUE, skip = 1, row.names = 1)
colnames(countdata) <- gsub(".bam", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("marked_dups.", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub(".sorted.marked_duplicates", "", colnames(countdata), fixed = T)
countdata <- countdata[ ,c(-1:-5)]
head(countdata)

metadata <- read.delim("metadata.txt", row.names = 1)
metadata

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = metadata,
                                 design = ~Group)

ddsMat <- DESeq(ddsMat)

results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)
summary(results)


mcols(results, use.names = T)


# Human genome database (Select the correct one)
library(org.Hs.eg.db)

# Add gene full name
results$description <- mapIds(x = org.Hs.eg.db,
                              keys = row.names(results),
                              column = "GENENAME",
                              keytype = "ENSEMBL",
                              multiVals = "first")

# Add gene symbol
results$symbol <- row.names(results)

# Add ENTREZ ID
results$entrez <- mapIds(x = org.Hs.eg.db,
                         keys = row.names(results),
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

# Add ENSEMBL
results$ensembl <- mapIds(x = org.Hs.eg.db,
                          keys = row.names(results),
                          column = "ENSEMBL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

# Subset for only significant genes (q < 0.05)
results_sig <- subset(results, padj < 0.05)
head(results_sig)

# Write normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(ddsMat), normalized = T), 
            file = 'normalized_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant normalized gene counts to a .txt file
write.table(x = counts(ddsMat[row.names(results_sig)], normalized = T), 
            file = 'normalized_counts_significant.txt', 
            sep = '\t', 
            quote = F, 
            col.names = NA)

# Write the annotated results table to a .txt file
write.table(x = as.data.frame(results), 
            file = "results_gene_annotated.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(results_sig), 
            file = "results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)


# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Plot PCA by column variable
plotPCA(ddsMat_rlog, intgroup = "Group", ntop = 500) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 5) + # Increase point size
  scale_y_continuous(limits = c(-2, 2.5)) + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 500 most variable genes") +
  theme(
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )


# Gather 30 significant genes and make matrix
mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:40, ]

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$Group), 
  Replicate = factor(colData(ddsMat_rlog)$Replicate),
  row.names = colData(ddsMat_rlog)$sampleid
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c(UHRF1_shRNA = "lightblue", Scramble_shRNA = "darkorange"),
  Replicate = c(Rep3 = "darkred", Rep2 = "forestgreen", Rep1 = "maroon")
)


#################################
col_data <- colData(ddsMat)


colData(ddsMat_rlog)$Group <- factor(colData(ddsMat_rlog)$Group, levels = c("UHRF1_shRNA", "Scramble_shRNA"))
colData(ddsMat_rlog)$Replicate <- factor(colData(ddsMat_rlog)$Replicate, levels = c("Rep1", "Rep2", "Rep3"))

annotation_col <- data.frame(
  Group = colData(ddsMat_rlog)$Group, 
  Replicate = colData(ddsMat_rlog)$Replicate,
  row.names = colData(ddsMat_rlog)$sampleid
)

#################################
# Make Heatmap with pheatmap function.
## See more in documentation for customization
pheatmap(mat = mat, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         annotation_colors = ann_colors,# Change the default colors of the annotations
         fontsize = 6.5, # Make fonts smaller
         cellwidth = 55, # Make the cells wider
         show_colnames = F)

# Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
## - Change pvalues to -log10 (1.3 = 0.05)
data <- data.frame(gene = row.names(results),
                   pval = -log10(results$padj), 
                   lfc = results$log2FoldChange)

# Remove any rows that have NA as an entry
data <- na.omit(data)

# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < 0 & data$pval > 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "nonsignificant"))

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

# Add ggplot2 layers
vol +   
  ggtitle(label = "UHRF1_shRNA vs Scramble_shRNA", subtitle = "gene expression") +
  geom_point(size = 1.5, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "green", Decreased = "red3", nonsignificant = "black")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("UHRF1_shRNA" / "Scramble_shRNA"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p")+
  theme(
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )# Scale yaxis due to large p-values



plotMA(results, ylim = c(-5, 5))

plotDispEsts(ddsMat)


# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Get gene with highest expression
top_gene <- rownames(results)[which.min(results$log2FoldChange)]

# Plot single gene
plotCounts(dds = ddsMat, 
           gene = top_gene, 
           intgroup = "Group", 
           normalized = T, 
           transform = T)


# Remove any genes that do not have any entrez identifiers
results_sig_entrez <- subset(results_sig, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
gene_matrix <- results_sig_entrez$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(gene_matrix) <- results_sig_entrez$entrez

# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(gene_matrix)


kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'human',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)

# Plot results
barplot(kegg_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)

go_enrichBP <- enrichGO(gene = names(gene_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

go_enrichMF <- enrichGO(gene = names(gene_matrix),
                        OrgDb = 'org.Hs.eg.db', 
                        readable = T,
                        ont = "MF",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)

# Plot results
barplot(go_enrichCC, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Cellular Component",
        font.size = 8)

# Load pathview
biocLite("pathview") ; library(pathview)

# Plot specific KEGG pathways (with fold change) 
## pathway.id : KEGG pathway identifier
pathview(gene.data = gene_matrix, 
         pathway.id = "05224", 
         species = "human")
