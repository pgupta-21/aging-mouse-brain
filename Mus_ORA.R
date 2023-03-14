# Over Representation Analysis
# Input: Differentially expressed genes from DEA
# Output: Results from clusterProfiler for over-represented terms

## ---------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(DOSE)
library(patchwork)

# Extracting significant genes based on over regulation or under regulation for 12m/6m
sigOE_6_12 <- rownames(sig_6_12)[(sig_6_12$log2FoldChange) >= 0.58] #only names are taken as the input for ORA
sigUE_6_12 <- rownames(sig_6_12)[(sig_6_12$log2FoldChange) <= -0.58]
dim(sig_6_12)

# enrich GO command
ora_O_6_12 <- enrichGO(gene = sigOE_6_12,
                     keyType = 'ENSEMBL',
                     universe = rownames(res_6_12_df),
                     OrgDb = 'org.Mm.eg.db',
                     ont = 'ALL',
                     readable = TRUE)

ora_U_6_12 <- enrichGO(gene = sigUE_6_12,
                     keyType = 'ENSEMBL',
                     universe = rownames(res_6_12_df),
                     OrgDb = 'org.Mm.eg.db',
                     ont = 'ALL',
                     readable = TRUE)


# Calculated enrichment ratio in excel by dividing GeneRatio by BgRatio and re-uploaded the resulting dataframe
# Barplots for visualizing results of 12m/6m for over and under regulated functional terms in our dataset
b1 <- ggplot(head(O6_12, 10), aes(x=Description, y=ES, label=sprintf('%0.2f', round(ES, digits = 2)), fill=p.adjust)) +
  ylab('Enrichment ratio') +
  geom_bar(position = position_dodge() ,stat = 'identity') +
  coord_flip() +
  ggtitle('Up-regulated GO terms for 12m/6m') +
  theme_classic() +
  theme(axis.text = element_text(face="bold", size = 12))

b2 <- ggplot(head(U6_12, 10), aes(x=Description, y=ES, label=sprintf('%0.2f', round(ES, digits = 2)), fill=p.adjust )) +
  ylab('Enrichment ratio') +
  geom_bar(position = position_dodge() ,stat = 'identity') +
  coord_flip() +
  ggtitle('Down-regulated GO terms for 12m/6m') +
  theme_classic() +
  theme(axis.text = element_text(face="bold", size = 12))

b1 + b2

# 24m/12m
sigOE_12_24 <- rownames(sig_12_24)[(sig_12_24$log2FoldChange) >= 0.58]
sigUE_12_24 <- rownames(sig_12_24)[(sig_12_24$log2FoldChange) <= -0.58]
ora_O_12_24 <- enrichGO(gene = sigOE_12_24,
                     keyType = 'ENSEMBL',
                     universe = rownames(res_12_24_df),
                     OrgDb = 'org.Mm.eg.db',
                     ont = 'ALL',
                     readable = TRUE)
ora_U_12_24 <- enrichGO(gene = sigUE_12_24,
                     keyType = 'ENSEMBL',
                     universe = rownames(res_12_24_df),
                     OrgDb = 'org.Mm.eg.db',
                     ont = 'ALL',
                     readable = TRUE)

b3 <- ggplot(head(O12_24, 10), aes(x=Description, y=ES, label=sprintf('%0.2f', round(ES, digits = 2)), fill=p.adjust)) +
  geom_bar(position = position_dodge() ,stat = 'identity') +
  ylab('Enrichment ratio') +
  coord_flip() +
  ggtitle('Up-regulated GO terms for 12m/24m') +
  theme_classic() +
  theme(axis.text = element_text(face="bold", size = 12))

b4 <- ggplot(head(U12_24, 10), aes(x=Description, y=ES, label=sprintf('%0.2f', round(ES, digits = 2)), fill=p.adjust )) +
  geom_bar(position = position_dodge() ,stat = 'identity') +
  ylab('Enrichment ratio') +
  coord_flip() +
  ggtitle('Down-regulated GO terms for 12m/24m') +
  theme_classic() +
  theme(axis.text = element_text(face="bold", size = 12))

b3 + b4

# 24m/6m
sigOE_6_24 <- rownames(sig_6_24)[(sig_6_24$log2FoldChange) >= 0.58]
sigUE_6_24 <- rownames(sig_6_24)[(sig_6_24$log2FoldChange) <= -0.58]
ora_O_6_24 <- enrichGO(gene = sigOE_6_24,
                     keyType = 'ENSEMBL',
                     universe = rownames(res_6_24_df),
                     OrgDb = 'org.Mm.eg.db',
                     ont = 'ALL',
                     readable = TRUE)
ora_U_6_24 <- enrichGO(gene = sigUE_6_24,
                     keyType = 'ENSEMBL',
                     universe = rownames(res_6_24_df),
                     OrgDb = 'org.Mm.eg.db',
                     ont = 'ALL',
                     readable = TRUE)

b5 <- ggplot(head(O6_24, 10), aes(x=Description, y=ES, label=sprintf('%0.2f', round(ES, digits = 2)), fill=p.adjust )) +
  geom_bar(position = position_dodge() ,stat = 'identity') +
  coord_flip() +
  ylab('Enrichment ratio') +
  ggtitle('Up-regulated GO terms for 24m/6m') +
  theme_classic() +
  theme(axis.text = element_text(face="bold", size = 12))

b6 <- ggplot(head(U6_24, 10), aes(x=Description, y=ES, label=sprintf('%0.2f', round(ES, digits = 2)), fill=p.adjust )) +
  geom_bar(position = position_dodge() ,stat = 'identity') +
  coord_flip() +
  ylab('Enrichment ratio') +
  ggtitle('Down-regulated GO terms for 24m/6m') +
  theme_classic() +
  theme(axis.text = element_text(face="bold", size = 12))

b5 + b6
