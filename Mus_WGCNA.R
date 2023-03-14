# Weighted Gene Correlation Network Analysis
# Input: Counts table from featureCounts
# Output: Modules from WGCNA and their analysis

## ---------------------------------------------------------------------------------------------------------------------------------------
library(WGCNA)
library(DESeq2)
library(CorLevelPlot)
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)
library(reshape)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GO.db)
library(rstatix)
library(dplyr)
library(broman)
library(clusterProfiler)
library(enrichplot)

# Initial dataset (fcount) from featureCounts which was cleaned during DEA
# Normalization
count_mat_fil <- as.matrix(fcount) #converting fcount to matrix
filtered.counts <- fcount %>%  
  filter_all(all_vars(. > 5)) #filtering
fcount.fil <- as.data.frame(filtered.counts) #converting to a dataframe
quan_counts <- normalize(filtered.counts) #quantile normalization with broman package
# Clean up of quantile normalized data
colnames(quan_counts) <- c('12_1', '12_2', '12_3', '12_4', '24_1', '24_2', '24_4', '24_5', '24_6', '6_1', '6_2', '6_3', '6_4')
quan_counts <- as.data.frame(quan_counts) 
rownames(quan_counts) <- rownames(fcount.fil)
# This was uploaded to excel and row-wise median normalization was performed - taking median of rows and dividing each column by the median
QUANT_RMN <- as.data.frame(new_qr) #uploading final normalized data to R
# Clean up
rownames(QUANT_RMN) <- QUANT_RMN$Column1
QUANT_RMN$Column1 <- NULL

# Log transformation of above counts
QUANT_RMN[] <- lapply(QUANT_RMN, function(x) as.numeric(as.character(x)))
fcount.rwm <- log(QUANT_RMN+1)
dim(fcount.rwm) #dimensions of our final data
boxplot(fcount.rwm, main = 'Boxplots after normalization')

# Checking for outliers
g <- goodSamplesGenes(t(fcount.rwm)) #transposing data for this command
summary(g) #checking which genes and/or samples contain outliers
g$allOK #TRUE or FALSE report for a quick check

table(g$goodGenes) #found in genes
table(g$goodSamples) #all samples are good to go

# Identifying the outlier genes
if (!g$allOK)
{
  if (sum(!g$goodGenes)>0)
    print(paste("Removing genes:", paste(rownames(fcount.rwm)[!g$goodGenes], collapse = ", ")));
  if (sum(!g$goodSamples)>0)
    print(paste("Removing samples:", paste(names(fcount.rwm)[!g$goodSamples], collapse = ", ")));
}

# Removing outlier genes
fcount.g <- fcount.rwm[g$goodGenes == TRUE,]
dim(fcount.g)

# Dataframe for final counts for later use
fcount.g.df <- as.data.frame(fcount.g)
fcount.g.df$entrez_id <- mapIds(org.Mm.eg.db, keys = rownames(fcount.g.df), keytype = 'ENSEMBL', column = 'ENTREZID')
head(fcount.g.df)

# Soft threshold power
temp_cor <- cor
cor <- WGCNA::cor
power <- c(c(1:10), seq(from = 12, to = 20, by = 2))
fcount.g <- t(fcount.g) #transposing the data for WGCNA command
sft.r <- pickSoftThreshold(fcount.g, powerVector = power, networkType = 'signed', verbose = 5)

sft_data.r <- sft.r$fitIndices #matrix for all the powers. Ideally we want a power with max R sq and min mean connectivity. That should give us a scale free topology.
head(sft_data.r) #sft data at a glance

# Plot for visualization
plot(sft.r$fitIndices[,1], -sign(sft.r$fitIndices[,3])*sft.r$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))

text(sft.r$fitIndices[,1], -sign(sft.r$fitIndices[,3]) * sft.r$fitIndices[,2],
     labels=power,cex=0.5,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft.r$fitIndices[,1], sft.r$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft.r$fitIndices[,1], sft.r$fitIndices[,5], labels=power, cex=0.5,col="red")

# Choosing soft power
fcount.g[] <- sapply(fcount.g, as.numeric) #conversion of normalized counts to numeric
soft_power.rw <- 14

# Adjacency matrix
adjacency.r <- adjacency(fcount.g, power = soft_power.rw)
TOMadj.r <- TOMsimilarity(adjacency.r) #turn adjacency into topological overlap matrix

# Dissimilarity matrix and heirarchical clustering of genes
disTOMadj.r <- 1- TOMadj.r #dissimilarity matrix
hgenetree.r <- hclust(as.dist(disTOMadj.r), method = 'average')
plot(hgenetree.r, xlab = '', sub = '', main = 'Gene clustering on TOM based dissimilarity', labels = FALSE, hang = 0.04)

# Creating a data set with various cutoff heights and their resultant modules to evaluate and choose a final cutoff
No.Genes = nrow(t(fcount.g))
ConsCutoffs = seq(from = 0.80, to = 0.98, by = 0.01)
No.Cutoffs = length(ConsCutoffs)
ConsensusColor = array(dim = c(No.Genes, No.Cutoffs))

for (cut in 1:No.Cutoffs)
{
    ConsensusColor[, cut] = cutreeDynamic(dendro = hgenetree.r , distM = disTOMadj.r ,
                                          deepSplit = 2, pamRespectsDendro = TRUE,
                                          minClusterSize = 200, cutHeight =  ConsCutoffs[cut])
}


plotDendroAndColors(hgenetree.r, ConsensusColor, c(ConsCutoffs), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# Choosing cutoff height
minModuleSize <- 200 #to make the modules larger
cutoff <- 0.93
# Module identification using dynamic tree cut
dynamic_modules.r <- cutreeDynamic(dendro = hgenetree.r, distM = disTOMadj.r, cutHeight = cutoff, deepSplit = 2, pamRespectsDendro = TRUE, minClusterSize = minModuleSize)

# Dynamic modules based on cut off heights
table(dynamic_modules.r)
dynamic_colors.r <- labels2colors(dynamic_modules.r) #changin numeric labels to colors
table(dynamic_colors.r)
plotDendroAndColors(hgenetree.r, dynamic_colors.r, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

# GO analysis of the modules using entrez IDs from the initial dataset and module colors
GOenr <- GOenrichmentAnalysis(dynamic_colors.r, fcount.g.df$entrez, organism = 'mouse', backgroundType = 'allGiven', nBestP = 2); #choosing to see the top 2 terms
tab = GOenr$bestPTerms[[4]]$enrichment
# Displaying the results from GO analysis
keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab1 = tab[, keepCols];
# Round the numeric columns to 2 decimal places
numCols = c(3, 4);
screenTab1[, numCols] = signif(apply(screenTab1[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab1[, 7] = substring(screenTab1[, 7], 1, 40)
# Shorten the column names
colnames(screenTab1) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab1) = NULL;
# Set the width of R’s output
options(width=95)
# Finally, display the enrichment table
screenTab1

# Significance tests and temporal patterns of the modules
# Combined gene data with normalized counts (non log counts) and module colors
QUANT_RMN['ENSMUSG00000119584',] #outlier gene identified initially
QUANT_RMN <- QUANT_RMN[!(row.names(QUANT_RMN) %in% c('ENSMUSG00000119584')), ] #removing the one outlier gene
fcount.g.df2 <- cbind(QUANT_RMN, dynamic_colors.r) #binding normalized counts with module colors
geneInfo <- as.data.frame(fcount.g.df2)
head(geneInfo)

# New metadata (6m, 12m, 24m)
coldata.r$replicate <- c('six','six','six','six','twelve','twelve','twelve','twelve','twentyfour','twentyfour','twentyfour','twentyfour','twentyfour')
coldata.r$samples <- rownames(coldata.r)
coldata.r

# Melted dataset
geneInfo_melted <-melt(geneInfo, id = 'dynamic_colors.r') %>% 
  dplyr::inner_join(coldata.r %>%
                      dplyr::select(samples, replicate),
                    by = c('variable' = 'samples')
                    )

# Significance
gene_separate <- geneInfo_melted[geneInfo_melted$dynamic_colors.r=='black',]
tukey <- tukey_hsd(gene_separate, value ~ replicate)
tukey <- as.data.frame(tukey)
tukey$color <- c('black','black','black')
tukey
write.csv(tukey, 'csv/black.csv') #test was done for all modules and all csv results were combined in excel

# Boxplots for all modules
ggplot(geneInfo_melted, aes(x = replicate, y = as.numeric(value), color = replicate)) + 
  ggtitle('Dynamic Modules') +
  ylim(0, 1.5) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~dynamic_colors.r) +
  theme_classic()


# Calculating module eigenegenes
dynamic_MElist.r <- moduleEigengenes(fcount.g, colors = dynamic_colors.r) #calculating eigengenes
dynamic_ME.r <- dynamic_MElist.r$eigengenes
dis_ME.r <- 1-cor(dynamic_ME.r)
ME_tree.r <- hclust(as.dist(dis_ME.r))
dynamic_ME.r

# Dendrogram of MEs
plot(ME_tree.r, main = 'Dynamic clustering of module eigengenes', xlab = '', sub = '')
ME_threshold <- 0.08 #to merge similar modules
abline(h = ME_threshold, col = 'red')

# Merging close modules depending on GO terms, significance and boxplots
merge_MED <- mergeCloseModules(fcount.g, dynamic_colors.r, cutHeight = ME_threshold, verbose = 3)

# Merged colors
merged_colors <- merge_MED$colors
table(merged_colors)

# Combined new dataset
fcount.merged.df <- cbind(QUANT_RMN, merged_colors) #binding normalized counts with merged module colors
geneInfo.merged <- as.data.frame(fcount.merged.df)
head(geneInfo.merged)

# Melted merged data-set
geneInfo_meltmerge <-melt(geneInfo.merged, id = 'merged_colors') %>% 
  dplyr::inner_join(coldata.r %>%
                      dplyr::select(samples, replicate),
                    by = c('variable' = 'samples')
                    )
# Boxplots for merged colours
ggplot(geneInfo_meltmerge, aes(x = replicate, y = value, color = replicate)) + 
  ggtitle('Merged Modules') +
  ylim(0, 1.5) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~merged_colors) +
  theme_classic()

# Significance
gene_separate2 <- geneInfo_meltmerge[geneInfo_meltmerge$merged_colors=='red',]
tukey2 <- tukey_hsd(gene_separate2, value ~ replicate)
tukey2 <- as.data.frame(tukey2)
tukey2$color <- c('red','red','red')
tukey2
write.csv(tukey2, 'csv/red.csv')

# GO analysis of merged modules
GOenr.m <- GOenrichmentAnalysis(merged_colors, fcount.g.df$entrez, organism = 'mouse', backgroundType = 'allGiven', nBestP = 2)
tab.m = GOenr.m$bestPTerms[[4]]$enrichment
write.table(tab.m, file = "GOEnrichmentTable_merged14-200.csv", sep = ",", quote = TRUE, row.names = FALSE) #to open in excel

# ORA of interesting merged modules
# Creating datasets with l2fc, padj and module colours
all_6_12 <- merge(geneInfo.merged, sig_6_12_df, by="row.names", all=TRUE) #sig_6_12_df is from our DEA analysis
all_12_24 <- merge(geneInfo.merged, sig_12_24_df, by="row.names", all=TRUE)
all_6_24 <- merge(geneInfo.merged, sig_6_24_df, by="row.names", all=TRUE)

# Subsetting according to module colours
small.merge <- subset(all_6_24, merged_colors == 'blue')
sigsmall_OE <- (small.merge$Row.names)[(small.merge$log2FoldChange) >= 0.58]
sigsmall_UE <- (small.merge$Row.names)[(small.merge$log2FoldChange) <= -0.58]
smallOE6_12 <- enrichGO(gene = sigsmall_OE,
                     keyType = 'ENSEMBL',
                     universe = rownames(res_6_24_df),
                     OrgDb = 'org.Mm.eg.db',
                     ont = 'ALL',
                     readable = TRUE)

smallUE6_12 <- enrichGO(gene = sigsmall_UE,
                     keyType = 'ENSEMBL',
                     universe = rownames(res_6_24_df),
                     OrgDb = 'org.Mm.eg.db',
                     ont = 'ALL',
                     readable = TRUE)

dotplot(smallOE6_12, showCategory = 5, split = 'ONTOLOGY', label_format = 50, title = 'Blue Module - Up regulated terms') + facet_grid(ONTOLOGY~., scale="free")
dotplot(smallUE6_12, showCategory = 5, split = 'ONTOLOGY', label_format = 50, title = 'Blue Module - Down regulated terms') + facet_grid(ONTOLOGY~., scale="free")
