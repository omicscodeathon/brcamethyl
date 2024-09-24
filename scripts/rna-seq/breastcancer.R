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
setwd("/Users/jkalami/Documents/geneX/brca/")


##   STEP:2 Preparing counts data and loading metadata    ##

# Import gene counts table
# - skip first row (general command info)
# - make row names the gene identifiers
brca_count <- read.table("final_counts.txt", header = TRUE, skip = 1, row.names = 1)

# Remove ...mapping. + 'Aligned.sortedByCoord.out.bam' from column identifiers
colnames(brca_count) <- gsub("...mapping.", "", colnames(brca_count), fixed = T)
colnames(brca_count) <- gsub("Aligned.sortedByCoord.out.bam", "", colnames(brca_count), fixed = T)

# Remove length/char columns
brca_count <- brca_count[ ,c(-1:-5)]

# Make sure ID's are correct
head(brca_count)

# load in sample meta info
brca_meta = read.delim("metadata.csv", sep = ",", stringsAsFactors = TRUE)

#make the Accession row names
rownames(brca_meta) =  brca_meta$Accession
brca_meta

#check if column names of counts matrix is the same as rownames of colData
all(colnames(brca_count) %in% rownames(brca_meta))

#check if they are in the same order
all(colnames(brca_count) == rownames(brca_meta))



##   STEP:3 Creating DESeq object from counts and metadata     ##

# - brca_count : count dataframe
# - brca_meta : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare

#creating the Deseq dataset
brca_ddsMat <- DESeqDataSetFromMatrix(countData = brca_count,
                                 colData = brca_meta,
                                 design = ~Treatment)


#inspect the levels in the deseq dataset 
brca_ddsMat$Treatment

#prefiltering of deseq dataset
nrow(brca_ddsMat)
brca_keep <- rowSums(counts(brca_ddsMat)) > 10
brca_ddsMat <- brca_ddsMat[brca_keep,]
nrow(brca_ddsMat)


##   STEP:4 Performing differential expression analysis    ##

#call the DESeq function on a DESeqDataSet object
brca_ddsMat <- DESeq(brca_ddsMat)

# Get results from testing with FDR adjust pvalues
brca_res <- results(brca_ddsMat, 
                   contrast=c("Treatment","UHRF1_shRNA","control_shRNA"),
                   pAdjustMethod = "fdr", 
                   alpha = 0.05)


# Generate summary of testing. 
summary(brca_res)


#annotate gene symbols using Human genome database
# Add gene full name
brca_res$description <- mapIds(x = org.Hs.eg.db,
                              keys = row.names(brca_res),
                              column = "GENENAME",
                              keytype = "ENSEMBL",
                              multiVals = "first")

# Add gene symbol
brca_res$symbol <- mapIds(org.Hs.eg.db, 
                         keys = row.names(brca_res), 
                         keytype = "ENSEMBL", 
                         column = "SYMBOL")


# Add ENTREZ ID
brca_res$entrez <- mapIds(x = org.Hs.eg.db,
                         keys = row.names(brca_res),
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

# Add ENSEMBL
brca_res$ensembl <- mapIds(x = org.Hs.eg.db,
                          keys = row.names(brca_res),
                          column = "ENSEMBL",
                          keytype = "ENSEMBL",
                          multiVals = "first")



# Subset for only significant genes (q < 0.05)
brca_res_sig <- subset(brca_res, padj < 0.05)
head(brca_res_sig)


#writing all important results to text files
# Write normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(brca_ddsMat), normalized = T), 
            file = 'mcf7_normalized_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant normalized gene counts to a .txt file
write.table(x = counts(brca_ddsMat[row.names(brca_res_sig)], normalized = T), 
            file = 'mcf7_normalized_counts_significant.txt', 
            sep = '\t', 
            quote = F, 
            col.names = NA)

# Write the annotated results table to a .txt file
write.table(x = as.data.frame(brca_res), 
            file = "mcf7_results_gene_annotated.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(brca_res_sig), 
            file = "mcf7_results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)



#PLOT PCA
# Convert all samples to rlog
brca_ddsMat_rlog <- rlog(brca_ddsMat, blind = FALSE)

# Plot PCA by column variable
plotPCA(brca_ddsMat_rlog, intgroup = "Treatment", ntop = 500) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 7) + # Increase point size
  scale_y_continuous(limits = c(-2.5, 2.5)) + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "MCF-7 Breast Cancer Cell line") +
  theme(
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank()  
)


# Gather 50 significant genes and make matrix
brca_mat <- assay(brca_ddsMat_rlog[row.names(brca_res_sig)])[1:50, ]

# convert to dataframe
bdf <- as.data.frame(brca_mat)

#add gene symbols
bdf$symbol = mapIds(org.Hs.eg.db, 
                    keys = row.names(bdf), 
                    keytype = "ENSEMBL", 
                    column = "SYMBOL")

# removes row with non-unique symbols
nrow(bdf)
bdf <- bdf[!duplicated(bdf$symbol),]
nrow(bdf)

# removes row with NA symbols
nrow(bdf)
bdf <- bdf[!is.na(bdf$symbol),]
nrow(bdf)

# make symbol rownames
rownames(bdf) =  bdf$symbol

# remove symbol column
bdf <- bdf[ ,c(-7)]

# Choose which column variables you want to annotate the columns by.
brca_annotation_col = data.frame(
  Treatment = factor(colData(brca_ddsMat_rlog)$Treatment), 
  Group = factor(colData(brca_ddsMat_rlog)$Group),
  row.names = colData(brca_ddsMat_rlog)$Accession
)

# Specify colors you want to annotate the columns by.
brca_ann_colors = list(
  Group = c(Rep3 = "#333CCC", Rep2 = "#660000", Rep1 = "#FFF000"),
  Treatment = c(UHRF1_shRNA = "#66FFFF", control_shRNA = "#66FF00")
)

# Make Heatmap with pheatmap function.
## See more in documentation for customization
pheatmap(mat = bdf, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = brca_annotation_col, # Add multiple annotations to the samples
         annotation_colors = brca_ann_colors,# Change the default colors of the annotations
         fontsize = 6.5, # Make fonts smaller
         cellwidth = 55, # Make the cells wider
         main = "50 significant DEGs (MCF-7)",
         show_colnames = F)


#Shrinking log2 fold changes to remove the noise associated with fold changes 
#coming from genes with low count levels
brca_resLFC <- lfcShrink(dds = brca_ddsMat, 
                    res = brca_res, 
                    coef="Treatment_UHRF1_shRNA_vs_control_shRNA", 
                    type="apeglm", 
)

# add gene symbols
brca_resLFC$symbol = mapIds(org.Hs.eg.db, 
                    keys = row.names(brca_resLFC), 
                    keytype = "ENSEMBL", 
                    column = "SYMBOL")

# create an enhanced volcano plot
EnhancedVolcano(brca_resLFC,
                lab = brca_resLFC$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 5),
                title = 'MCF-7 UHRF1_shRNA versus control_shRNA',
                pCutoff = 10e-8,
                FCcutoff = 0.6,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'blue', 'red3'),
                colAlpha = 1)

# plot MA 
plotMA(brca_res, 
       ylim = c(-5, 5),
       colNonSig = "#333333",
       colSig = "red",
       colLine = "#FF6666",)

# plot Dispersion
plotDispEsts(brca_ddsMat)


# Get gene with highest expression
top_gene <- rownames(brca_res)[which.min(brca_res$log2FoldChange)]

# Plot single gene
plotCounts(dds = brca_ddsMat, 
           gene = top_gene, 
           intgroup = "Treatment", 
           normalized = T, 
           transform = T)

##### GENE SET ENRICHMENT ANALYSIS #####
# Remove any genes that do not have any entrez identifiers
brca_res_sig_entrez <- subset(brca_res_sig, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
brca_gene_matrix <- brca_res_sig_entrez$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(brca_gene_matrix) <- brca_res_sig_entrez$ensembl

# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(brca_gene_matrix)

# sort the list in decreasing order (required for clusterProfiler)
brca_gene_matrix = sort(brca_gene_matrix, decreasing = TRUE)



# Gene Onotology
brca_gse <- gseGO(geneList=brca_gene_matrix, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

# plot gse GO dotplot
dotplot(brca_gse, title = "GO MCF-7 cell line", showCategory=10, split=".sign") + facet_grid(.~.sign)


# plot enrichment map
cnetplot(brca_gse, categorySize="pvalue", foldChange=brca_gene_matrix, showCategory = 2)


# KEGG PATHWAY
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
original_gene_list <- brca_res_sig$log2FoldChange
ids<-bitr(names(brca_gene_matrix), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = brca_res_sig[brca_res_sig$ensembl %in% dedup_ids$ENSEMBL,]

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
        title = "Enriched Pathways (MCF-7 cell line)", 
        split=".sign") + facet_grid(.~.sign)

