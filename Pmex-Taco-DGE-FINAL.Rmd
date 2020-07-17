---
title: "Poecilia mexicana RNA-seq analysis"
author: "Kerry McGowan"
output:
  html_notebook:
    fig_caption: yes
    toc: yes
    toc_depth: 3
    toc_float: yes
  pdf_document:
    toc: yes
    toc_depth: '3'
  html_document:
    keep_md: TRUE
---

# Initate the project

Install libraries
```{r, results='hide'}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

install.packages("xlsx", "gplots", "data.table", "dplyr", "ggplot2", "d3heatmap")
install.packages ("d3heatmap")
```

Load libraries
```{r, results='hide'}
library(DESeq2)
library(xlsx)
library(gplots)
library(data.table)
library(dplyr)
library(pcaExplorer)
library(ggplot2)
options(stringsAsFactors = FALSE);
```

Load gene counts matrix and sample sheet
```{r}
cts = read.csv("gene_count_matrix-2019-10-01_no_STRG.csv", header = T, row.names = 1)
# Subset gene counts matrix for Tacotalpa samples
taco_cts = cts[,88:99]
write.csv(taco_cts, "taco_cts.csv")

key = read.xlsx("PmexRNA_samples_2010.xlsx", sheetName = "RawReadData")
# Subset sample sheet for Tacotalpa samples
taco_key = key[24:35,]
```

# Differential expression analysis with DEseq
##  Tacotalpa sulfidic vs non-sulfidic gills

Create DESeqDataSet object
```{r}
dds = DESeqDataSetFromMatrix(countData = taco_cts, colData = taco_key, design = ~ Sulfur.1)
```

Specify the reference level (e.g., the control group) for factors (default is alphabetical)
```{r}
dds$Sulfur.1 = relevel(dds$Sulfur.1, ref = "fresh")
```

Estimate size factors, dispersion, fit Negative Binomial GLM and calculate Wald statistics.
Values that have NA for p-value are flagged as being the Cook's distance cutoff. (See http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#approach-to-count-outliers "Approach to outliers").
```{r}
# parametric = default
dds = DESeq(dds, fitType = "parametric")
res = results(dds)
res
# sulfur vs fresh means that the estimates are log2(sulfur/fresh)

# local
dds_local = DESeq(dds, fitType = "local")
res_local = results(dds_local)
res_local
# sulfur vs fresh means that the estimates are log2(sulfur/fresh)

# mean
dds_mean = DESeq(dds, fitType = "mean")
res_mean = results(dds_mean)
res_mean
# sulfur vs fresh means that the estimates are log2(sulfur/fresh)
```

Plot padj parametric
```{r}
pdf("ParametricFit.pdf")
hist(res$padj[res$baseMean >1], main="Parametric", ylim = c(0,3500), xlab = "Adjusted p-values")
while (!is.null(dev.list()))
dev.off()
```

Plot padj local. Note how this one remains the flattest as the x-axis increases (this is ideal)
```{r}
pdf("LocalFit.pdf")
hist(res_local$padj[res_local$baseMean >1], main="Local", ylim = c(0,3500), xlab = "Adjusted p-values")
while (!is.null(dev.list()))
dev.off()
```

Plot padj mean
```{r}
pdf("MeanFit.pdf")
hist(res_mean$padj[res_mean$baseMean >1], main="Mean", ylim = c(0,3500), xlab = "Adjusted p-values")
while (!is.null(dev.list()))
dev.off()
```

3-panel figure of parametric, local, and mean
```{r}
pdf("3PanelFit.pdf")
par(mfrow=c(3,1))
hist(res$padj[res$baseMean >1], main="Parametric", ylim = c(0,3500), xlab = "Adjusted p-values")
hist(res$padj[res$baseMean >1], main="Local", ylim = c(0,3500), xlab = "Adjusted p-values")
hist(res_mean$padj[res_mean$baseMean >1], main="Mean", ylim = c(0,3500), xlab = "Adjusted p-values")
while (!is.null(dev.list()))
dev.off()
```

Local fit looks the best, so all following analyses use this fit.

Plot PCA (Principal Components Analysis)
```{r}
# rlog transformed data
rld <- rlog(dds_local)
pdf("Sulfur_vs_Fresh-localFit-PCA.pdf")
plotPCA(rld, intgroup = "Sulfur.1")
while (!is.null(dev.list()))
dev.off()
```

Filtering by adjusted p-value
```{r}
resOrdered <- res_local[order(res_local$padj),]

resSig <- subset(resOrdered, padj < 0.01)
resSig
```

MA-plot for all results, red points have adjusted p-value < 0.01
```{r}
plotMA(res_local, ylim=c(-2,2))
```

MA-plot for significant results, red points have adjusted p-value < 0.01
```{r}
plotMA(resSig, ylim=c(-2,2))
```

Plot counts
```{r}
# smallest padj
plotCounts(dds, gene=which.min(res_local$padj), intgroup="Sulfur.1")
```

```{r}
# largest padj
plotCounts(dds, gene=which.max(res_local$padj), intgroup="Sulfur.1")
```

Save differential expression results as a CSV file
```{r}
write.csv(as.data.frame(resOrdered), file="Sulfur_vs_Fresh-localFit.csv")
```

Merge differential expression results with gene IDs 
```{r}
resOrderedDF <- as.data.frame(resOrdered)

# Move gene IDs to column 1
resOrderedDF <- setDT(resOrderedDF, keep.rownames = TRUE)[]

# Rename first column
names(resOrderedDF)[1] <- "geneID"

# Read in PmexGeneNameMatching from Passow
PmexGeneNameMatching <- read.csv(file = "PmexGeneNameMatching.csv")

# Make row 1 the headers for the columns
names(PmexGeneNameMatching) <- as.matrix(PmexGeneNameMatching[1, ])
PmexGeneNameMatching <- PmexGeneNameMatching[-1, ]
PmexGeneNameMatching[] <- lapply(PmexGeneNameMatching, function(x) type.convert(as.character(x)))

# Merge by gene ID
merged <- merge(x = resOrderedDF, y = PmexGeneNameMatching, by.x = "geneID", by.y = "gene.ID", all.x = TRUE)
merged <- select(merged, geneID,gene.name,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj)
```

Merge differential expression results with annotations
```{r}
# Read in annotations file from Passow et al. Mol Ecol 2017 "The roles of plasticity and evolutionary change in shaping gene expression variation in natural populations of extremophile fish"
anno <- read.csv(file = "Table_S2-Poecilia_mexicana_Annotations.csv")

# Merge merged_S_vs_NS and anno by gene ID
merged2 <- merge(x = merged, y = anno, by.x = "geneID", by.y = "gene.ID", all.x = TRUE)
# note not all genes have annotations, which is why 2 geneID columns appear, gene.name.x is more complete

# Subset data frame to only contain relevant columns
merged_subset <- select(merged2, geneID,gene.name.x,Subject.sequence.ID,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,Protein.annotations)
dim(merged_subset)

# Remove duplicate rows
merged_unique = unique(merged_subset)
dim(merged_unique)

# Sort by padj
sorted = merged_unique[order(padj),]

# Write CSV of results
write.csv(x=sorted, file = "Sulfur_vs_Fresh-localFit_WITH_ANNOTATIONS.csv")
```

Merge differential expression results with Ensembl IDs. Ensembl IDs were downloaded (6/18/2020) from Ensembl's BioMart [database = Ensembl Genes 100, dataset = Shortfin molly genes (P_mexicana-1.0)] in a table containing attributes: (column 1) GENE: Ensembl: Gene stable ID; (column 2) EXTERNAL: External References: NCBI gene (formerly Entrezgene) accession.
```{r}
ensembl = read.csv(file = "mart_export.txt")

final = merge(x = sorted, y = ensembl, by.x = "gene.name.x", by.y = "NCBI.gene..formerly.Entrezgene..accession", all.x = TRUE)

# Reorder columns
final = select(final, geneID,gene.name.x,Gene.stable.ID,Subject.sequence.ID,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,Protein.annotations)

# Sort by padj
final = final[order(padj),]

write.csv(x=final, file = "Sulfur_vs_Fresh-localFit_WITH_ANNOTATIONS_AND_ENSEMBL.csv")
```
