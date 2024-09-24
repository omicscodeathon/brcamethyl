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
library(enrichplot)
library(EnhancedVolcano)

#set working directory
setwd("/Users/jkalami/Documents/geneX/leukemia")
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

#create Kasumi-1 countdata
kas_count <- leu_count[ ,c(-5:-8)]

# Make sure ID's are correct
head(kas_count)

# load in sample meta info
leu_meta = read.delim("metadata.csv", sep = ",", stringsAsFactors = TRUE)

#make the Accession row names
rownames(leu_meta) =  leu_meta$Accession
leu_meta

#create Kasumi-1 coldata
kas_meta <- leu_meta[c(-5:-8),]
kas_meta

#check if column names of counts matrix is the same as rownames of colData
all(colnames(kas_count) %in% rownames(kas_meta))

#check if they are in the same order
all(colnames(kas_count) == rownames(kas_meta))



##   STEP:3 Creating DESeq object from counts and metadata     ##

# - kas_count : count dataframe
# - kas_meta : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare

#creating the Deseq dataset
kas_ddsMat <- DESeqDataSetFromMatrix(countData = kas_count,
                                     colData = kas_meta,
                                     design = ~Treatment)


#inspect the levels in the deseq dataset 
kas_ddsMat$Treatment

#prefiltering of deseq dataset
nrow(kas_ddsMat)
kas_keep <- rowSums(counts(kas_ddsMat)) > 10
kas_ddsMat <- kas_ddsMat[kas_keep,]
nrow(kas_ddsMat)


##   STEP:4 Performing differential expression analysis    ##

#call the DESeq function on a DESeqDataSet object
kas_ddsMat <- DESeq(kas_ddsMat)

# Get results from testing with FDR adjust pvalues
kas_res <- results(kas_ddsMat, 
                   contrast=c("Treatment","UHRF1_shRNA","control_shRNA"),
                   pAdjustMethod = "fdr", 
                   alpha = 0.05)


# Generate summary of testing. 
summary(kas_res)


#annotate gene symbols using Human genome database
# Add gene full name
kas_res$description <- mapIds(x = org.Hs.eg.db,
                              keys = row.names(kas_res),
                              column = "GENENAME",
                              keytype = "ENSEMBL",
                              multiVals = "first")

# Add gene symbol
kas_res$symbol <- mapIds(org.Hs.eg.db, 
                         keys = row.names(kas_res), 
                         keytype = "ENSEMBL", 
                         column = "SYMBOL")


# Add ENTREZ ID
kas_res$entrez <- mapIds(x = org.Hs.eg.db,
                         keys = row.names(kas_res),
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

# Add ENSEMBL
kas_res$ensembl <- mapIds(x = org.Hs.eg.db,
                          keys = row.names(kas_res),
                          column = "ENSEMBL",
                          keytype = "ENSEMBL",
                          multiVals = "first")



# Subset for only significant genes (q < 0.05)
kas_res_sig <- subset(kas_res, padj < 0.05)
head(kas_res_sig)


#writing all important results to text files
# Write normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(kas_ddsMat), normalized = T), 
            file = 'Kasumi-1_normalized_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant normalized gene counts to a .txt file
write.table(x = counts(kas_ddsMat[row.names(kas_res_sig)], normalized = T), 
            file = 'Kasumi-1_normalized_counts_significant.txt', 
            sep = '\t', 
            quote = F, 
            col.names = NA)

# Write the annotated results table to a .txt file
write.table(x = as.data.frame(kas_res), 
            file = "Kasumi-1_results_gene_annotated.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(kas_res_sig), 
            file = "Kasumi-1_results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)



#PLOT PCA
# Convert all samples to rlog
kas_ddsMat_rlog <- rlog(kas_ddsMat, blind = FALSE)

# Plot PCA by column variable
plotPCA(kas_ddsMat_rlog, intgroup = "Treatment", ntop = 500) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 7) + # Increase point size
  scale_y_continuous(limits = c(-7.5, 7.5)) + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Kasumi-1 Leukemia Cell line") +
  theme(
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank()  
)


# Gather 50 significant genes and make matrix
kas_mat <- assay(kas_ddsMat_rlog[row.names(kas_res_sig)])[1:50, ]

# convert to dataframe
kdf <- as.data.frame(kas_mat)

#add gene symbols
kdf$symbol = mapIds(org.Hs.eg.db, 
                    keys = row.names(kdf), 
                    keytype = "ENSEMBL", 
                    column = "SYMBOL")

# removes row with non-unique symbols
nrow(kdf)
kdf <- kdf[!duplicated(kdf$symbol),]
nrow(kdf)

# removes row with NA symbols
nrow(kdf)
kdf <- kdf[!is.na(kdf$symbol),]
nrow(kdf)

# make symbol rownames
rownames(kdf) =  kdf$symbol

# remove symbol column
kdf <- kdf[ ,c(-5)]

# Choose which column variables you want to annotate the columns by.
kas_annotation_col = data.frame(
  Treatment = factor(colData(kas_ddsMat_rlog)$Treatment), 
  Group = factor(colData(kas_ddsMat_rlog)$Group),
  row.names = colData(kas_ddsMat_rlog)$Accession
)

# Specify colors you want to annotate the columns by.
kas_ann_colors = list(
  Group = c(kas1 = "orange", kas2 = "brown"),
  Treatment = c(UHRF1_shRNA = "#66FFFF", control_shRNA = "#66FF00")
)

# Make Heatmap with pheatmap function.
## See more in documentation for customization
pheatmap(mat = kdf, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = kas_annotation_col, # Add multiple annotations to the samples
         annotation_colors = kas_ann_colors,# Change the default colors of the annotations
         fontsize = 6.5, # Make fonts smaller
         cellwidth = 55, # Make the cells wider
         main = "50 significant DEGs (Kasumi-1)",
         show_colnames = F)


#Shrinking log2 fold changes to remove the noise associated with fold changes 
#coming from genes with low count levels
kas_resLFC <- lfcShrink(dds = kas_ddsMat, 
                        res = kas_res, 
                        coef="Treatment_UHRF1_shRNA_vs_control_shRNA", 
                        type="apeglm", 
)

# add gene symbols
kas_resLFC$symbol = mapIds(org.Hs.eg.db, 
                           keys = row.names(kas_resLFC), 
                           keytype = "ENSEMBL", 
                           column = "SYMBOL")

# create an enhanced volcano plot
EnhancedVolcano(kas_resLFC,
                lab = kas_resLFC$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 5),
                title = 'Kasumi-1 UHRF1_shRNA versus control_shRNA',
                pCutoff = 10e-8,
                FCcutoff = 0.6,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'blue', 'red3'),
                colAlpha = 1)

# plot MA 
plotMA(kas_res, 
       ylim = c(-5, 5),
       colNonSig = "#333333",
       colSig = "red",
       colLine = "#FF6666")

# plot Dispersion
plotDispEsts(kas_ddsMat)


# Get gene with highest expression
top_gene <- rownames(kas_res)[which.min(kas_res$log2FoldChange)]

# Plot single gene
plotCounts(dds = kas_ddsMat, 
           gene = top_gene, 
           intgroup = "Treatment", 
           normalized = T, 
           transform = T)


# Remove any genes that do not have any entrez identifiers
kas_res_sig_entrez <- subset(kas_res_sig, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
kas_gene_matrix <- kas_res_sig_entrez$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(kas_gene_matrix) <- kas_res_sig_entrez$ensembl

# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(kas_gene_matrix)

# sort the list in decreasing order (required for clusterProfiler)
kas_gene_matrix = sort(kas_gene_matrix, decreasing = TRUE)



# Gene Set Enrichment
kas_gse <- gseGO(geneList=kas_gene_matrix, 
                 ont ="ALL", 
                 keyType = "ENSEMBL", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Hs.eg.db, 
                 pAdjustMethod = "none")

# plot gse GO dotplot
dotplot(kas_gse, title = "GO Kasumi-1 cell line",showCategory=10, split=".sign") + facet_grid(.~.sign)


# plot enrichment map
cnetplot(kas_gse, categorySize="pvalue", foldChange=kas_gene_matrix, showCategory = 1)


# KEGG PATHWAY
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
original_gene_list <- kas_res_sig$log2FoldChange
ids<-bitr(names(kas_gene_matrix), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = kas_res_sig[kas_res_sig$ensembl %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "hsa",
               #nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, 
        title = "Enriched Pathways (Kasumi-1 cell line)", 
        split=".sign") + facet_grid(.~.sign)

