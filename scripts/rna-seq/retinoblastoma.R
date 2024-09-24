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
setwd("/Users/jkalami/Documents/geneX/retino/")


##   STEP:2 Preparing counts data and loading metadata    ##

# Import gene counts table
# - skip first row (general command info)
# - make row names the gene identifiers
ret_count <- read.table("final_count-ret.txt", header = TRUE, skip = 1, row.names = 1)

# Remove ...mapping. + 'Aligned.sortedByCoord.out.bam' from column identifiers
colnames(ret_count) <- gsub("...mapping.", "", colnames(ret_count), fixed = T)
colnames(ret_count) <- gsub("Aligned.sortedByCoord.out.bam", "", colnames(ret_count), fixed = T)

# Remove length/char columns
ret_count <- ret_count[ ,c(-1:-5)]

# Make sure ID's are correct
head(ret_count)

# load in sample meta info
ret_meta = read.delim("metadata.csv", sep = ",", stringsAsFactors = TRUE)

#make the Accession row names
rownames(ret_meta) =  ret_meta$Accession
ret_meta

#check if column names of counts matrix is the same as rownames of colData
all(colnames(ret_count) %in% rownames(ret_meta))

#check if they are in the same order
all(colnames(ret_count) == rownames(ret_meta))



##   STEP:3 Creating DESeq object from counts and metadata     ##

# - ret_count : count dataframe
# - ret_meta : sample metadata in the dataframe with row names as sampleID's
# - design : The design of the comparisons to use. 
#            Use (~) before the name of the column variable to compare

#creating the Deseq dataset
ret_ddsMat <- DESeqDataSetFromMatrix(countData = ret_count,
                                      colData = ret_meta,
                                      design = ~Treatment)


#inspect the levels in the deseq dataset 
ret_ddsMat$Treatment

#prefiltering of deseq dataset
nrow(ret_ddsMat)
ret_keep <- rowSums(counts(ret_ddsMat)) > 10
ret_ddsMat <- ret_ddsMat[ret_keep,]
nrow(ret_ddsMat)


##   STEP:4 Performing differential expression analysis    ##

#call the DESeq function on a DESeqDataSet object
ret_ddsMat <- DESeq(ret_ddsMat)

# Get results from testing with FDR adjust pvalues
ret_res <- results(ret_ddsMat, 
                    contrast=c("Treatment","UHRF1_shRNA","control_shRNA"),
                    pAdjustMethod = "fdr", 
                    alpha = 0.05)


# Generate summary of testing. 
summary(ret_res)


#annotate gene symbols using Human genome database
# Add gene full name
ret_res$description <- mapIds(x = org.Hs.eg.db,
                               keys = row.names(ret_res),
                               column = "GENENAME",
                               keytype = "ENSEMBL",
                               multiVals = "first")

# Add gene symbol
ret_res$symbol <- mapIds(org.Hs.eg.db, 
                          keys = row.names(ret_res), 
                          keytype = "ENSEMBL", 
                          column = "SYMBOL")


# Add ENTREZ ID
ret_res$entrez <- mapIds(x = org.Hs.eg.db,
                          keys = row.names(ret_res),
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")

# Add ENSEMBL
ret_res$ensembl <- mapIds(x = org.Hs.eg.db,
                           keys = row.names(ret_res),
                           column = "ENSEMBL",
                           keytype = "ENSEMBL",
                           multiVals = "first")



# Subset for only significant genes (q < 0.05)
ret_res_sig <- subset(ret_res, padj < 0.05)
head(ret_res_sig)


#writing all important results to text files
# Write normalized gene counts to a .txt file
write.table(x = as.data.frame(counts(ret_ddsMat), normalized = T), 
            file = 'y79_normalized_counts.txt', 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant normalized gene counts to a .txt file
write.table(x = counts(ret_ddsMat[row.names(ret_res_sig)], normalized = T), 
            file = 'y79_normalized_counts_significant.txt', 
            sep = '\t', 
            quote = F, 
            col.names = NA)

# Write the annotated results table to a .txt file
write.table(x = as.data.frame(ret_res), 
            file = "y79_results_gene_annotated.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)

# Write significant annotated results table to a .txt file
write.table(x = as.data.frame(ret_res_sig), 
            file = "y79_results_gene_annotated_significant.txt", 
            sep = '\t', 
            quote = F,
            col.names = NA)



#PLOT PCA
# Convert all samples to rlog
ret_ddsMat_rlog <- rlog(ret_ddsMat, blind = FALSE)

# Plot PCA by column variable
plotPCA(ret_ddsMat_rlog, intgroup = "Treatment", ntop = 500) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 7) + # Increase point size
  scale_y_continuous(limits = c(-10, 10)) + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Y79 Retinoblastoma Cell line") +
  theme(
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank()  
  )


# Gather 50 significant genes and make matrix
ret_mat <- assay(ret_ddsMat_rlog[row.names(ret_res_sig)])[1:50, ]

# convert to dataframe
rdf <- as.data.frame(ret_mat)

#add gene symbols
rdf$symbol = mapIds(org.Hs.eg.db, 
                    keys = row.names(rdf), 
                    keytype = "ENSEMBL", 
                    column = "SYMBOL")

# removes row with non-unique symbols
nrow(rdf)
rdf <- rdf[!duplicated(rdf$symbol),]
nrow(rdf)

# removes row with NA symbols
nrow(rdf)
rdf <- rdf[!is.na(rdf$symbol),]
nrow(rdf)

# make symbol rownames
rownames(rdf) =  rdf$symbol

# remove symbol column
rdf <- rdf[ ,c(-5)]

# Choose which column variables you want to annotate the columns by.
ret_annotation_col = data.frame(
  Treatment = factor(colData(ret_ddsMat_rlog)$Treatment), 
  Group = factor(colData(ret_ddsMat_rlog)$Group),
  row.names = colData(ret_ddsMat_rlog)$Accession
)

# Specify colors you want to annotate the columns by.
ret_ann_colors = list(
  Group = c(YT1 = "#333CCC", YT2 = "#660000", YS1 = "orange", YS2 = "brown"),
  Treatment = c(UHRF1_shRNA = "#66FFFF", control_shRNA = "#66FF00")
)

# Make Heatmap with pheatmap function.
## See more in documentation for customization
pheatmap(mat = rdf, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = ret_annotation_col, # Add multiple annotations to the samples
         annotation_colors = ret_ann_colors,# Change the default colors of the annotations
         fontsize = 6.5, # Make fonts smaller
         cellwidth = 55, # Make the cells wider
         main = "50 significant DEGs (Y79)",
         show_colnames = F)


#Shrinking log2 fold changes to remove the noise associated with fold changes 
#coming from genes with low count levels
ret_resLFC <- lfcShrink(dds = ret_ddsMat, 
                         res = ret_res, 
                         coef="Treatment_UHRF1_shRNA_vs_control_shRNA", 
                         type="apeglm", 
)

# add gene symbols
ret_resLFC$symbol = mapIds(org.Hs.eg.db, 
                            keys = row.names(ret_resLFC), 
                            keytype = "ENSEMBL", 
                            column = "SYMBOL")

# create an enhanced volcano plot
EnhancedVolcano(ret_resLFC,
                lab = ret_resLFC$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 5),
                title = 'Y79 UHRF1_shRNA versus control_shRNA',
                pCutoff = 10e-8,
                FCcutoff = 0.6,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'blue', 'red3'),
                colAlpha = 1)

# plot MA 
plotMA(ret_res, 
       ylim = c(-5, 5),
       colNonSig = "#333333",
       colSig = "red",
       colLine = "#FF6666",)

# plot Dispersion
plotDispEsts(ret_ddsMat)


# Get gene with highest expression
top_gene <- rownames(ret_res)[which.min(ret_res$log2FoldChange)]

# Plot single gene
plotCounts(dds = ret_ddsMat, 
           gene = top_gene, 
           intgroup = "Treatment", 
           normalized = T, 
           transform = T)

##### GENE SET ENRICHMENT ANALYSIS #####
# Remove any genes that do not have any entrez identifiers
ret_res_sig_entrez <- subset(ret_res_sig, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
ret_gene_matrix <- ret_res_sig_entrez$log2FoldChange

# Add the entrezID's as names for each logFC entry
names(ret_gene_matrix) <- ret_res_sig_entrez$ensembl

# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(ret_gene_matrix)

# sort the list in decreasing order (required for clusterProfiler)
ret_gene_matrix = sort(ret_gene_matrix, decreasing = TRUE)


# Gene Onotology
ret_gse <- gseGO(geneList=ret_gene_matrix, 
                  ont ="ALL", 
                  keyType = "ENSEMBL", 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "none")

# plot gse GO dotplot
dotplot(ret_gse, title = "GO Y79 cell line", showCategory=10, split=".sign") + facet_grid(.~.sign)


# plot enrichment map
cnetplot(ret_gse, categorySize="pvalue", foldChange=ret_gene_matrix, showCategory = 2)


# KEGG PATHWAY
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
original_gene_list <- ret_res_sig$log2FoldChange
ids<-bitr(names(ret_gene_matrix), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = ret_res_sig[ret_res_sig$ensembl %in% dedup_ids$ENSEMBL,]

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
        title = "Enriched Pathways (Y79 cell line)", 
        split=".sign") + facet_grid(.~.sign)

