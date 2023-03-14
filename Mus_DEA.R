# Differential Gene Expression Analysis
# Input: Counts table from the featureCounts
# Output: DESeq2 downstream analysis results

## -----------------------------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(org.Mm.eg.db)
library(biomaRt)
library(EnhancedVolcano)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(pheatmap)
library(RColorBrewer)
library(viridis)


# Importing the dataset for 13 fastq files after getting a summary from featureCounts as fcount.
# Data cleanup
head(fcount)

rownames(fcount) <- fcount$V1
fcount$V1 <- NULL
fcount$V2 <- NULL
fcount$V3 <- NULL
fcount$V4 <- NULL
fcount$V5 <- NULL
fcount$V6 <- NULL
colnames(fcount) <- c('12_1', '12_2', '12_3', '12_4', '24_1', '24_2', '24_4', '24_5', '24_6', '6_1', '6_2', '6_3', '6_4')
fcount <- fcount[-1,]
fcount$Length <- NULL #cleaning up data
rownames(fcount) <- fcount$gene_id
fcount$gene_id <- NULL

sapply(fcount, class) #class of each column
fcount[] <- lapply(fcount, function(x) as.numeric(as.character(x))) #changing class to numeric

# First visualization
boxplot(fcount) #boxplot of raw data
# The boxplot is not very pretty so, I plotted the log +1 of the counts.
count_mat <- as.matrix(fcount) #conversion to a matrix
count_mat <- log(count_mat+1) #+1 to overcome Inf values
pheatmap(cor(count_mat), color=inferno(10), cluster_rows = TRUE, cluster_cols = TRUE)
boxplot(count_mat, main = 'Boxplots before normalization') #boxplot before normalization


# Forming the metadata for DESeq2 functions
replicate <- factor(c('12m', '12m', '12m', '12m','24m', '24m', '24m', '24m', '24m', '6m', '6m', '6m', '6m')) #creating a factor with three levels
coldata <- data.frame(row.names = colnames(fcount), replicate) #creating metadata for our samples
all(colnames(fcount) == rownames(coldata)) #ensuring the column names from our filtered data are same as rownames from the metadata

# DESeq object and filtering
# featureCounts output is a count matrix; hence DESeqDatasetFromMatrix command is used in this case.
ds <- DESeqDataSetFromMatrix(countData = fcount, colData = coldata, design = ~replicate) #creating a DESeq object with design as the replicate.
ds <- ds[rowSums(counts(ds) >= 5) == 13,] #removing genes with counts less than 5 in all of the samples
dim(ds)

ds <- estimateSizeFactors(ds)

# Normalization seen with boxplot and heatmap
normalizedcounts <- counts(ds, normalized=TRUE) #normalizing the ds with counts function
pheatmap(cor(normalizedcounts), color=mako(10), cluster_rows = T, cluster_cols = T, show_colnames = T, show_rownames = T, main = 'Biological replicate clustering')
normalizedcounts <- log(normalizedcounts+1) #conversion to log +1 for boxplot
boxplot(normalizedcounts, main = 'Boxplots after normalization') #boxplot after normalization of data and filtering

## -----------------------------------------------------------------------------
# DESeq Analysis
ds <- DESeq(ds) #performing DESeq analysis of ds data
head(ds)

# vst for normalization for PCA
vscounts <- vst(ds, blind = FALSE)
plotPCA(vscounts, intgroup = 'replicate') + theme_classic() #PCA plot

# Dispersion estimates
plotDispEsts(ds) #dispersion estimate

# Results 12m/6m
res_6_12 <- results(ds, contrast = c('replicate', '12m', '6m')) #results function to compare between 12m and 6m
# Significant results
sig_6_12 <- na.omit(res_6_12) #omitting NA values
sig_6_12 <- sig_6_12[(sig_6_12$baseMean >= 10) & (abs(sig_6_12$log2FoldChange) >= 0.58) & (sig_6_12$padj <= 0.05),] #significant data from comparing 12m and 6m
plotMA(res_6_12) #MA plot for 12m and 6m
# Creating ensembl, entrez and uniprot ids for significant genes
sig_6_12_df <- as.data.frame(sig_6_12)
sig_6_12_df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(sig_6_12_df), keytype = 'ENSEMBL', column = 'SYMBOL')
sig_6_12_df$entrez_id <- mapIds(org.Mm.eg.db, keys = rownames(sig_6_12_df), keytype = 'ENSEMBL', column = 'ENTREZID')
sig_6_12_df$uniprot <- mapIds(org.Mm.eg.db, keys = rownames(sig_6_12_df), keytype = 'ENSEMBL', column = 'UNIPROT')
# Volcano plot
res_6_12_df <- as.data.frame(res_6_12)
res_6_12_df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(res_6_12_df), keytype = 'ENSEMBL', column = 'SYMBOL')
res_6_12_df$uniprot <- mapIds(org.Mm.eg.db, keys = rownames(res_6_12_df), keytype = 'ENSEMBL', column = 'UNIPROT')
res_6_12_df$entrez_id <- mapIds(org.Mm.eg.db, keys = rownames(res_6_12_df), keytype = 'ENSEMBL', column = 'ENTREZID')
EnhancedVolcano(res_6_12_df, x = 'log2FoldChange', y = 'padj', lab = res_6_12_df$symbol, title = 'DEGs for 12m/6m', gridlines.major = FALSE, gridlines.minor = FALSE)

# Results 24m/12m
res_12_24 <- results(ds, contrast = c('replicate', '24m', '12m')) #results function to compare between 24m and 12m
sig_12_24 <- na.omit(res_12_24) #omitting NA values
sig_12_24 <- sig_12_24[(sig_12_24$baseMean >= 10) & (abs(sig_12_24$log2FoldChange) >= 0.58) & (sig_12_24$padj <= 0.05), ] #significant data from comparing 24m and 12m
plotMA(res_12_24)
# Creating ensembl, entrez and uniprot ids for significant genes
sig_12_24_df <- as.data.frame(sig_12_24)
sig_12_24_df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(sig_12_24_df), keytype = 'ENSEMBL', column = 'SYMBOL')
sig_12_24_df$entrez_id <- mapIds(org.Mm.eg.db, keys = rownames(sig_12_24_df), keytype = 'ENSEMBL', column = 'ENTREZID')
sig_12_24_df$uniprot <- mapIds(org.Mm.eg.db, keys = rownames(sig_12_24_df), keytype = 'ENSEMBL', column = 'UNIPROT')
# Volcano plot
res_12_24_df <- as.data.frame(res_12_24)
res_12_24_df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(res_12_24_df), keytype = 'ENSEMBL', column = 'SYMBOL')
res_12_24_df$uniprot <- mapIds(org.Mm.eg.db, keys = rownames(res_12_24_df), keytype = 'ENSEMBL', column = 'UNIPROT')
res_12_24_df$entrez_id <- mapIds(org.Mm.eg.db, keys = rownames(res_12_24_df), keytype = 'ENSEMBL', column = 'ENTREZID')
EnhancedVolcano(res_12_24_df, x = 'log2FoldChange', y = 'padj', lab = res_12_24_df$symbol, title = 'DEGs for 24m/12m')

# Results 24m/6m
res_6_24 <- results(ds, contrast = c('replicate', '24m', '6m')) ##results function to compare between 24m and 6m
sig_6_24 <- na.omit(res_6_24) #omitting NA values
sig_6_24 <- sig_6_24[(sig_6_24$baseMean >= 10) & (abs(sig_6_24$log2FoldChange) >= 0.58) & (sig_6_24$padj <= 0.05), ]
plotMA(res_6_24) #MA plot for 24m and 6m
# Creating ensembl, entrez and uniprot ids for significant genes
sig_6_24_df <- as.data.frame(sig_6_24)
sig_6_24_df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(sig_6_24_df), keytype = 'ENSEMBL', column = 'SYMBOL')
sig_6_24_df$entrez_id <- mapIds(org.Mm.eg.db, keys = rownames(sig_6_24_df), keytype = 'ENSEMBL', column = 'ENTREZID')
sig_6_24_df$uniprot <- mapIds(org.Mm.eg.db, keys = rownames(sig_6_24_df), keytype = 'ENSEMBL', column = 'UNIPROT')
# Volcano plot
res_6_24_df <- as.data.frame(res_6_24)
res_6_24_df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(res_6_24_df), keytype = 'ENSEMBL', column = 'SYMBOL')
res_6_24_df$uniprot <- mapIds(org.Mm.eg.db, keys = rownames(res_6_24_df), keytype = 'ENSEMBL', column = 'UNIPROT')
res_6_24_df$entrez_id <- mapIds(org.Mm.eg.db, keys = rownames(res_6_24_df), keytype = 'ENSEMBL', column = 'ENTREZID')
EnhancedVolcano(res_6_24_df, x = 'log2FoldChange', y = 'padj', lab = res_6_24_df$symbol, title = 'DEGs for 24m/6m')
