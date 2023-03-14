# Differential Transcript Expression Analysis
# Input: Quantification output from Salmon
# Output: DESeq2 downstream analysis results

## ---------------------------------------------------------------------------------------------------------------------------------------
library(tximport)
library(dplyr)
library(DESeq2)
library(readr)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(biomaRt)
library(EnhancedVolcano)


# Metadata
rownames(coldata.2.1) <- coldata.2.1$samples
coldata.2.1

# Transcript files
sample_files <- paste0(pull(coldata.2.1, 'samples'), '/quant.sf')
names(sample_files) <- pull(coldata.2.1, 'samples')
sample_files

# Importing the dataset using tximport
txi_data <- tximport(files = sample_files,
         type = 'salmon',
         txOut = TRUE,
         ignoreTxVersion = TRUE,
         )

head(txi_data$counts)
all(colnames(txi_data$counts) == rownames(coldata.2.1))

# DESeq object and filtering
txi <- DESeqDataSetFromTximport(txi = txi_data,
                                colData = coldata.2.1,
                                design = ~replicate.r)
txi <- txi[rowSums(counts(txi) >= 5) == 13,]
dim(txi)

# DESeq analysis
txi <- DESeq(txi)
head(txi)

# PCA plot
vsTcounts <- vst(txi, blind = FALSE)
plotPCA(vsTcounts, intgroup = 'replicate')

# DEA
# Results 12m/6m
tx_6_12 <- results(txi, contrast = c('replicate.r', '12m', '6m')) #results function to compare between 12m and 6m
tx_6_12

tx_6_12_df <- as.data.frame(tx_6_12)
rownames(tx_6_12_df) <- sub("\\.\\d+", "", rownames(tx_6_12_df)) #removing the version ID from transcripts
tx_6_12_df$ensmusTid <- rownames(tx_6_12_df)
tx_6_12_df <- merge(tx_6_12_df, gene_map) #merging the transcript IDs with their mapped gene IDs (done in command line)
tx_6_12_df <- tx_6_12_df[!duplicated(tx_6_12_df$ensmusTid),]
rownames(tx_6_12_df) <- tx_6_12_df$ensmusTid
tx_6_12_df$ensmusTid <- NULL
tx_6_12_df$symbol <- mapIds(org.Mm.eg.db, keys = tx_6_12_df$ensmusGid, keytype = 'ENSEMBL', column = 'SYMBOL')
EnhancedVolcano(tx_6_12_df, x = 'log2FoldChange', y = 'padj', lab = tx_6_12_df$symbol, title = 'DETs for 12m/6m')

sigtx_6_12_df <- na.omit(tx_6_12_df)
sigtx_6_12_df <- sigtx_6_12_df[(sigtx_6_12_df$baseMean >= 10) & (abs(sigtx_6_12_df$log2FoldChange) >= 0.58) & (sigtx_6_12_df$padj <= 0.05),] #significant data from comparing 12m and 6m
dim(sigtx_6_12_df)

# Results 24m/12m
tx_12_24 <- results(txi, contrast = c('replicate.r', '24m', '12m'))
tx_12_24_df <- as.data.frame(tx_12_24)
rownames(tx_12_24_df) <- sub("\\.\\d+", "", rownames(tx_12_24_df))
tx_12_24_df$ensmusTid <- rownames(tx_12_24_df)
tx_12_24_df <- merge(tx_12_24_df, gene_map)
tx_12_24_df <- tx_12_24_df[!duplicated(tx_12_24_df$ensmusTid),]
rownames(tx_12_24_df) <- tx_12_24_df$ensmusTid
tx_12_24_df$ensmusTid <- NULL
tx_12_24_df$symbol <- mapIds(org.Mm.eg.db, keys = tx_12_24_df$ensmusGid, keytype = 'ENSEMBL', column = 'SYMBOL')
EnhancedVolcano(tx_12_24_df, x = 'log2FoldChange', y = 'padj', lab = tx_12_24_df$symbol, title = 'DETs for 24m/12m')

sigtx_12_24_df <- na.omit(tx_12_24_df)
sigtx_12_24_df <- sigtx_12_24_df[(sigtx_12_24_df$baseMean >= 10) & (abs(sigtx_12_24_df$log2FoldChange) >= 0.58) & (sigtx_12_24_df$padj <= 0.05),]
dim(sigtx_12_24_df)

# Results 24m/6m
tx_6_24 <- results(txi, contrast = c('replicate.r', '24m', '6m'))
tx_6_24_df <- as.data.frame(tx_6_24)
rownames(tx_6_24_df) <- sub("\\.\\d+", "", rownames(tx_6_24_df))
tx_6_24_df$ensmusTid <- rownames(tx_6_24_df)
tx_6_24_df <- merge(tx_6_24_df, gene_map)
tx_6_24_df <- tx_6_24_df[!duplicated(tx_6_24_df$ensmusTid),]
rownames(tx_6_24_df) <- tx_6_24_df$ensmusTid
tx_6_24_df$ensmusTid <- NULL
tx_6_24_df$symbol <- mapIds(org.Mm.eg.db, keys = tx_6_24_df$ensmusGid, keytype = 'ENSEMBL', column = 'SYMBOL')
EnhancedVolcano(tx_6_24_df, x = 'log2FoldChange', y = 'padj', lab = tx_6_24_df$symbol, title = 'DETs for 24m/6m')
head(tx_6_24_df)

sigtx_6_24_df <- na.omit(tx_6_24_df)
sigtx_6_24_df <- sigtx_6_24_df[(sigtx_6_24_df$baseMean >= 10) & (abs(sigtx_6_24_df$log2FoldChange) >= 0.58) & (sigtx_6_24_df$padj <= 0.05),]
dim(sigtx_6_24_df)
