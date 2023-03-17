# aging-mouse-brain
RNA Seq Analysis

To address the gaps in transcriptomic studies of brain aging, we decided to study the whole brain transcriptome
of mice at three different ages: six, twelve and twenty-four months. My aim is to conduct an exploratory analysis
of the aging mouse brain for identification of changes at the transcriptome level. The above mentioned aim was 
achieved by the following specific steps-

  - Gathering of RNA-Seq data from the three ages and pre-processing the reads. Performing thorough quality checks
    and then aligning the FASTQ files to their reference genome with HISAT2. Identifying the length of the reads 
    and quantifying the reads with Subread.
  - Performing Differential Expression Analysis using DESeq2 and identifying various differentially expressed genes
    for all three conditions and three comparisons (twelve versus six month, twenty-four versus twelve month and 
    twenty-four versus six month of age).
  - In depth analysis of the DEGs to identify up-regulated and down-regulated pathways from results of each comparison
    using Over-Representation Analysis, Gene Set Enrichment Analysis and Pathway Analysis using clusterProfiler and 
    DOSE.
  - Weighted Gene Correlation Network Analysis with WGCNA for identifying modules with co-expressed genes and 
    executing functional annotation- GSEA and ORA, for interesting modules.
  - From transcripts quantified with Salmon, performing downstream Differential Transcript Expression Analysis using 
    tximport and DESeq2. Further distinguishing isoform switches and alternative splicing patterns in transcripts as
    a part of Differential Transcript Usage Analysis using IsoformSwitchAnalzeR and lastly, comparison of results
    obtained.
    
The following repository contains Bash and R code for the above mentioned steps.
