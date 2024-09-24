library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)
library(apeglm)
library(tibble)
library(dplyr)
library(tidyverse)

#set working directory
setwd("/Users/jkalami/Documents/geneX/featurel")
getwd()


##   STEP:2 Preparing counts data and loading metadata    ##

# Import gene counts table
# - skip first row (general command info)
# - make row names the gene identifiers
leu_count <- read.table("final_count-leu.txt", header = TRUE, skip = 1, row.names = 1)

# Remove ...mapping. + 'Aligned.sortedByCoord.out.bam' from column identifiers
colnames(leu_count) <- gsub("mapping.", "", colnames(leu_count), fixed = T)
colnames(leu_count) <- gsub("Aligned.sortedByCoord.out.bam", "", colnames(leu_count), fixed = T)

# Remove length/char columns
leu_count <- leu_count[ ,c(-1:-5)]

# Make sure ID's are correct
head(leu_count)

# load in sample meta info
leu_meta = read.delim("metadata.csv", sep = ",", stringsAsFactors = TRUE)

#make the Accession row names
rownames(leu_meta) =  leu_meta$Accession
leu_meta
#metadata = metadata[,-c(1)] 

#check if column names of counts matrix is the same as rownames of colData
all(colnames(leu_count) %in% rownames(leu_meta))

#check if they are in the same order
all(colnames(leu_count) == rownames(leu_meta))



##   STEP:3 Creating DESeq object from counts and metadata     ##

# - leu_count : count dataframe
# - leu_meta : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare

#creatind the Deseq dataset
leu_ddsMat <- DESeqDataSetFromMatrix(countData = leu_count,
                                     colData = leu_meta,
                                     design = ~Treatment)


#inspect the levels in the deseq dataset 
leu_ddsMat$Treatment

#prefiltering of deseq dataset
nrow(leu_ddsMat)
leu_keep <- rowSums(counts(leu_ddsMat)) > 10
leu_ddsMat <- leu_ddsMat[leu_keep,]
nrow(leu_ddsMat)


##   STEP:4 Performing differential expression analysis    ##

#call the DESeq function on a DESeqDataSet object
leu_ddsMat <- DESeq(leu_ddsMat)

# Get results from testing with FDR adjust pvalues
leu_res <- results(leu_ddsMat, 
                   contrast=c("Treatment","UHRF1_shRNA","ctrl_shRNA"),
                   pAdjustMethod = "fdr", 
                   alpha = 0.05)


# Generate summary of testing. 
summary(leu_res)


#annotate gene symbols using Human genome database
# Add gene full name
results$description <- mapIds(x = org.Hs.eg.db,
                              keys = row.names(results),
                              column = "GENENAME",
                              keytype = "ENSEMBL",
                              multiVals = "first")

# Add gene symbol
results$symbol <- mapIds(org.Hs.eg.db, 
                         keys = row.names(results), 
                         keytype = "ENSEMBL", 
                         column = "SYMBOL")


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

#writing all important results to text files
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



#PLOT PCA
# Convert all samples to rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# Plot PCA by column variable
plotPCA(ddsMat_rlog, intgroup = "Treatment", ntop = 500) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 7) + # Increase point size
  scale_y_continuous(limits = c(-2.5, 2.5)) + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 500 most variable genes") +
  theme(
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank()  
  )


# Gather 50 significant genes and make matrix
mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:50, ]

# Choose which column variables you want to annotate the columns by.
annotation_col = data.frame(
  Treatment = factor(colData(ddsMat_rlog)$Treatment), 
  Group = factor(colData(ddsMat_rlog)$Group),
  row.names = colData(ddsMat_rlog)$Accession
)

# Specify colors you want to annotate the columns by.
ann_colors = list(
  Group = c(Rep3 = "maroon", Rep2 = "darkorange", Rep1 = "turquoise"),
  Treatment = c(UHRF1_shRNA = "dodgerblue4", SCRAMBLE_shRNA = "#d8b365")
)

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


#Shrinking log2 fold changes to remove the noise associated with fold changes 
#coming from genes with low count levels
resLFC <- lfcShrink(dds = ddsMat, 
                    res = results, 
                    coef="Treatment_UHRF1_shRNA_vs_SCRAMBLE_shRNA", 
                    type="apeglm", 
)


EnhancedVolcano(toptable = resLFC,
                x = "log2FoldChange",
                y = "padj",
                lab = rownames(resLFC),
                xlim = c(-3, +4),
                ylim = c(0,87.5),
                pointSize = 2.0,
                title = "UHRF1_shRNA versus SCRAMBLE_shRNA",
                legendLabels=c(
                  'Not significant',
                  'Log2 fold-change ',
                  'p-value',
                  'p-value & Log2 fold change')
)

