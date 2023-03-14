# Gene Set Enrichment Analysis
# Input: Differentially expressed genes from DEA
# Output: Results from clusterProfiler for enriched terms

## ---------------------------------------------------------------------------------------------------------------------------------------
library(org.Mm.eg.db)
library(biomaRt)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(DOSE)

# Geneset 12m/6m
sig_6_12 <- sig_6_12[order(-sig_6_12$log2FoldChange),] #ordering of input list in descending order by log2FC
genes_6_12 <- sig_6_12$log2FoldChange
names(genes_6_12) <- rownames(sig_6_12) #naming the log2fc values
# GSEA command
gsea_6_12 <- gseGO(genes_6_12, keyType = 'ENSEMBL', ont = 'ALL', OrgDb = 'org.Mm.eg.db', eps = 1e-300) #maxing out the eps so that p values are precise
# Visualization of the first term via gseaplot
gseaplot(gsea_6_12, geneSetID = 1)
# Dotplot 12m/6m
dotplot(gsea_6_12, showCategory=9, split='ONTOLOGY', x = 'NES', title = 'GSEA for 12m/6m')

# Geneset 24m/12m
sig_12_24 <- sig_12_24[order(-sig_12_24$log2FoldChange),]
genes_12_24 <- sig_12_24$log2FoldChange
names(genes_12_24) <- rownames(sig_12_24)
gsea_12_24 <- gseGO(genes_12_24, keyType = 'ENSEMBL', OrgDb = 'org.Mm.eg.db', eps = 1e-300, ont = 'ALL')
gseaplot(gsea_12_24, geneSetID = 1)
dotplot(gsea_12_24, showCategory=10, split='ONTOLOGY', x = 'NES', title = 'GSEA for 24m/12m') #+ facet_grid(.~.sign)

# Geneset 24m/6m
sig_6_24 <- sig_6_24[order(-sig_6_24$log2FoldChange),]
genes_6_24 <- sig_6_24$log2FoldChange
names(genes_6_24) <- rownames(sig_6_24)
gsea_6_24 <- gseGO(genes_6_24, keyType = 'ENSEMBL', OrgDb = 'org.Mm.eg.db', eps = 1e-300, ont = 'ALL')
gseaplot(gsea_6_24, geneSetID = 1)
dotplot(gsea_6_24, showCategory=10, split='ONTOLOGY',  x = 'NES', title = 'GSEA for 24m/6m')
