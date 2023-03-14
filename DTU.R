# Differential Transcript Usage Analysis
# Input: Quantification output from Salmon
# Output: DEXSeq results from downstream analysis

## ---------------------------------------------------------------------------------------------------------------------------------------
library(IsoformSwitchAnalyzeR)

# Importing salmon results
salmonQuant <- importIsoformExpression(
  parentDir = "files/"
)

# Abundance matrix
head(salmonQuant$abundance)

# Count matrix
head(salmonQuant$counts)

# Design formula
myDesign <- data.frame(sampleID = colnames(salmonQuant$abundance)[-1],
                       condition = c("twelve","twelve","twelve","twelve","twentyFour","twentyFour","twentyFour","twentyFour","twentyFour","six","six","six","six"))
myDesign

# Creating switchAnalyseRList
sSwitchList <- importRdata(
  isoformCountMatrix = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix = myDesign,
  isoformExonAnnoation = "gencode.vM31.annotation.gtf.gz",
  isoformNtFasta = "gencode.vM31.transcripts.fa.gz",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = FALSE
)

summary(sSwitchList)

# Filtering
sSwitchList.fil <- preFilter(
  switchAnalyzeRlist = sSwitchList,
  geneExpressionCutoff = 5,
  isoformExpressionCutoff = 5,
  removeSingleIsoformGenes = TRUE
)

head(sSwitchList.fil$isoformFeatures)

# Identifying isoform switches via DEXSeq
aSwitchList <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = sSwitchList.fil,
  reduceToSwitchingGenes = TRUE
)

# Re-written switches
extractSwitchSummary(aSwitchList)

# Predicting alternative splicing
asSwitchList <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = aSwitchList,
  quiet = TRUE
)

extractSplicingSummary(
  asSwitchList,
  asFractionTotal = FALSE,
  plotGenes = FALSE
)

# Summarizing uneven usage
splicingEnrichment <- extractSplicingEnrichment(
  asSwitchList,
  splicingToAnalyze = "all",
  returnResult = TRUE,
  returnSummary = TRUE
)

# Individual splice types
# How does the isoform usage of all isoforms utilizing a particular splicing type change?
extractSplicingGenomeWide(
  asSwitchList,
  featureToExtract = "all",
  splicingToAnalyze = c("ATSS","A3"),
  plot = TRUE,
  returnResult = FALSE
)

# AS sites dataset
asSwitchList$AlternativeSplicingAnalysis

# Isoform switch dataset
asSwitchList$isoformSwitchAnalysis

# Merging the two data-sets above
splicingAnalysis <- merge(asSwitchList$isoformSwitchAnalysis, asSwitchList$AlternativeSplicingAnalysis, by="isoform_id")
dim(splicingAnalysis)

# filtering by pvalue cutoff
splicingAnalysis <- subset(splicingAnalysis, padj <= 0.05)
splicingAnalysis$isoform_id <- sub("\\.\\d+", "", splicingAnalysis$isoform_id)
dim(splicingAnalysis)

# For 6m and 12m 
as_switch_6_12 <- subset(splicingAnalysis, (condition_1=="six" & condition_2=="twelve"))
# Combining with significant data from DTE analysis of those genes with two or more transcripts
tu.6.12 <- merge(as_switch_6_12, duplicates_6_12, by="isoform_id")

# For 12m and 24m
as_switch_12_24 <- subset(splicingAnalysis, (condition_1=="twelve" & condition_2=="twentyFour"))
tu.12.24 <- merge(as_switch_12_24, duplicates_12_24, by="isoform_id")

# For 6m and 24m
as_switch_6_24 <- subset(splicingAnalysis, (condition_1=="six" & condition_2=="twentyFour"))
tu.6.24 <- merge(as_switch_6_24, duplicates_6_24, by="isoform_id")
