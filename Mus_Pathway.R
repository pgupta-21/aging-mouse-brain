# Pathway Analysis
# Input: Differentially expressed genes from DEA
# Output: Results from clusterProfiler for enriched and over-represented pathways

## ---------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(org.Mm.eg.db)
library(biomaRt)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(DOSE)
library(pathview)

# 12m/6m genelist
# Creating genelist for pathway analysis using Uniprot IDs
uni_6_12 <- sig_6_12_df$log2FoldChange #extracting the log2foldchange values from our dataset above
names(uni_6_12) <- sig_6_12_df$uniprot #naming the numbers with uniprot ids
uni_6_12 <- na.omit(uni_6_12) #omitting na
uni_6_12 = sort(uni_6_12, decreasing = TRUE) #sorting into descending order
uni_6_12 <- uni_6_12[!duplicated(names(uni_6_12))] #removing duplicates
head(uni_6_12)

# gseKEGG function- GSEA of KEGG pathways
k6_12 <- gseKEGG(geneList = uni_6_12, organism = 'mmu', keyType = 'uniprot', pvalueCutoff = 0.05, verbose = TRUE, eps = 1e-300)
pathway_6_12 <- as.data.frame(k6_12)
dotplot(k6_12, x = 'enrichmentScore' ,showCategory = 30)
ridgeplot(k6_12, showCategory = 20) + labs(x = "enrichmentScore") + ggtitle("Enriched pathways for 12m/6m") + theme_classic()

# enrichKEGG function- ORA of KEGG pathways
sig_p_OE_6_12 <- sig_6_12_df$uniprot[(sig_6_12_df$log2FoldChange) > 0.58] #input for ORA
sig_p_UE_6_12 <- sig_6_12_df$uniprot[(sig_6_12_df$log2FoldChange) <= -0.58]
kp_OE_6_12 <- enrichKEGG(gene = sig_p_OE_6_12,
                 universe = res_6_12_df$uniprot,
                 organism = 'mmu',
                 keyType = 'uniprot',
                 pvalueCutoff = 0.05)
kp_UE_6_12 <- enrichKEGG(gene = sig_p_UE_6_12,
                 universe = res_6_12_df$uniprot,
                 organism = 'mmu',
                 keyType = 'uniprot',
                 pvalueCutoff = 0.05)

enrichplot::dotplot(kp_OE_6_12)
enrichplot::dotplot(kp_UE_6_12)

# Visualizing few pathways from 12m/6m
# Creating a new dataset for pathview function for successfully mapped uniprot ids
uni_6_12_df <- as.data.frame(uni_6_12) #dataframe of successfully mapped uniprot ids
uni_6_12_df$entrez_id <- mapIds(org.Mm.eg.db, keys = rownames(uni_6_12_df), keytype = 'UNIPROT', column = 'ENTREZID') #conversion to entrez ids for pathview function
rownames(uni_6_12_df) <- NULL
rownames(uni_6_12_df) <- uni_6_12_df$entrez_id #renaming the rows
uni_6_12_df$entrez_id <- NULL
pathview(gene.data = uni_6_12_df,
         pathway.id = 'mmu00190',
         species = 'mmu',
         kegg.native = T,
         out.suffix = 'oxiphosp',
         limit = list(gene = 2, cpd = 1)) #example pathway from 12m/6m for oxidative phosphorylation

## ---------------------------------------------------------------------------------------------------------------------------------------
# 24m/12m
uni_12_24 <- sig_12_24_df$log2FoldChange #extracting the log2foldchange values from our dataset
names(uni_12_24) <- sig_12_24_df$uniprot #naming the numbers with uniprot ids
uni_12_24 <- na.omit(uni_12_24) #omitting na
uni_12_24 = sort(uni_12_24, decreasing = TRUE) #sorting into descending order
uni_12_24 <- uni_12_24[!duplicated(names(uni_12_24))] #removing duplicates

k12_24 <- gseKEGG(geneList = uni_12_24,
                  organism = 'mmu',
                  keyType = 'uniprot',
                  pvalueCutoff = 0.5,
                  verbose = TRUE) #pvaluecutoff increased to 1 because there are no significnatly enriched pathways between 12 and 24 samples

ridgeplot(k12_24, showCategory = 20) + labs(x = "enrichmentScore") + ggtitle("Enriched pathways for 24m/12m") + theme_classic()

sig_p_OE_12_24 <- sig_12_24_df$uniprot[(sig_12_24_df$log2FoldChange) > 0.58]
sig_p_UE_12_24 <- sig_12_24_df$uniprot[(sig_12_24_df$log2FoldChange) <= 0.58]
kp_OE_12_24 <- enrichKEGG(gene = sig_p_OE_12_24,
                 universe = res_12_24_df$uniprot,
                 organism = 'mmu',
                 keyType = 'uniprot',
                 pvalueCutoff = 0.05)
kp_UE_12_24 <- enrichKEGG(gene = sig_p_UE_12_24,
                 universe = res_12_24_df$uniprot,
                 organism = 'mmu',
                 keyType = 'uniprot',
                 pvalueCutoff = 0.05)

enrichplot::dotplot(kp_OE_12_24)
enrichplot::dotplot(kp_UE_12_24)

## ---------------------------------------------------------------------------------------------------------------------------------------
# 24m/6m
uni_6_24 <- sig_6_24_df$log2FoldChange #extracting the log2foldchange values from our dataset above
names(uni_6_24) <- sig_6_24_df$uniprot #naming the numbers with uniprot ids
uni_6_24 <- na.omit(uni_6_24) #omitting na
uni_6_24 = sort(uni_6_24, decreasing = TRUE) #sorting into descending order
uni_6_24 <- uni_6_24[!duplicated(names(uni_6_24))] #removing duplicates

k6_24 <- gseKEGG(geneList = uni_6_24, organism = 'mmu', keyType = 'uniprot', pvalueCutoff = 0.05, verbose = TRUE, eps = 1e-300)
pathway_6_24 <- as.data.frame(k6_24)
ridgeplot(k6_24, showCategory = 20) + labs(x = "enrichmentScore") + ggtitle("Enriched pathways for 24m/6m") + theme_classic()

sig_p_OE_6_24 <- sig_6_24_df$uniprot[(sig_6_24_df$log2FoldChange) > 0.58]
sig_p_UE_6_24 <- sig_6_24_df$uniprot[(sig_6_24_df$log2FoldChange) <= 0.58]
kp_OE_6_24 <- enrichKEGG(gene = sig_p_OE_6_24,
                         universe = res_6_24_df$uniprot,
                         organism = 'mmu',
                         keyType = 'uniprot',
                         pvalueCutoff = 0.05)
kp_UE_6_24 <- enrichKEGG(gene = sig_p_UE_6_24,
                         universe = res_6_24_df$uniprot,
                         organism = 'mmu',
                         keyType = 'uniprot',
                         pvalueCutoff = 0.05)

enrichplot::dotplot(kp_OE_6_24) #no enriched term!
enrichplot::dotplot(kp_UE_6_24)

uni_6_24_df <- as.data.frame(uni_6_24) #dataframe of successfully mapped uniprot ids
uni_6_24_df$entrez_id <- mapIds(org.Mm.eg.db, keys = rownames(uni_6_24_df), keytype = 'UNIPROT', column = 'ENTREZID') #conversion to entrez ids for pathview function
rownames(uni_6_24_df) <- NULL
rownames(uni_6_24_df) <- uni_6_24_df$entrez_id #renaming the rows
uni_6_24_df$entrez_id <- NULL
pathview(gene.data = uni_6_24_df,
         pathway.id = 'mmu00190',
         species = 'mmu',
         kegg.native = T,
         out.suffix = 'oxyphosp.2',
         limit = list(gene = 2, cpd = 1)) #oxidative phosphorylation pathway for 24m/6m
