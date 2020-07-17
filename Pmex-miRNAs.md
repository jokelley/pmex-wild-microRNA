---
title: "Poecilia mexicana microRNA analysis"
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
### Load libraries
```{r, results='hide'}
library(DESeq2)
library(reshape2)
library(gplots)
options(stringsAsFactors = FALSE);
```

### Load counts file and sample sheet
```{r}
mirnas = read.csv("Pme_prost_counts_191128_td.csv", header = T, row.names = 1)
key = read.csv("sample_key.csv", header = T)
```

## MDS for entire dataset
```{r}
cds = DESeqDataSetFromMatrix( countData = mirnas, colData = key, design = ~ Population )
dds <- DESeq(cds)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
head(assay(vsd), 3)
plotPCA(vsd, intgroup=c("Population", "Sex"))
```

## Remove outlier
### Replot MDS
```{r}
mirnas = mirnas[,-which(colnames(mirnas) %in% c("sample_591"))]
key = key[-which(key$Sample.ID == "sample_591"),]

cds = DESeqDataSetFromMatrix( countData = mirnas, colData = key, design = ~ Population )
dds <- DESeq(cds)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
head(assay(vsd), 3)
plotPCA(vsd, intgroup=c("Population", "Sex"))
```

### Extract normalized counts for heatmap
```{r}
norm_counts = counts(dds, normalized = TRUE)
write.csv(norm_counts, file = "NormalizedCounts.csv", quote = F)
```

# Differential expression analysis with DEseq
##  Tacotalpa sulfidic vs non-sulfidic gills

## Multifactor analysis 
Calculate normalization factors
The calcNormFactors function normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes.

```{r}
cds = DESeqDataSetFromMatrix( countData = mirnas, colData = key, design = ~ Sex + Population )
dds <- DESeq(cds)

ddsMF <- dds
design(ddsMF) <- formula(~ Sex + Population)
ddsMF <- DESeq(ddsMF)

resMF <- results(ddsMF)
head(resMF)

hist(resMF$padj[resMF$baseMean > 1])

resOrdered <- resMF[order(resMF$pvalue),]
summary(resMF, alpha = 0.01)
sum(resMF$padj < 0.01, na.rm=TRUE)

resSig <- subset(resOrdered, padj < 0.01)
resSig
rownames(resSig)

write.csv(as.data.frame(resOrdered), file="Population_PSO_vs_Bonita-MFdesign.csv")

plotMA(resMF, alpha = 0.01, ylim=c(-2,2), cex=0.75)

# looking at SEX
resMFSex <- results(ddsMF, contrast=c("Sex", "F", "M"))
head(resMFSex)
resOrderedSex <- resMFSex[order(resMFSex$pvalue),]
summary(resMFSex, alpha = 0.01)
sum(resMFSex$padj < 0.01, na.rm=TRUE)
#pme-miR-8160a-5p
plotCounts(cds, gene="pme-miR-8160a-5p", intgroup="Sex")
```

## Combining datasets

```{r}
targets = read.csv(file = "pmex-miranda-3prime-predictedTargets.csv", header = T)
genes = read.csv(file = "TableS4-mRNA Expression.csv", header = T)
mirna.de = read.csv(file = "TableS3-miRNA Expression.csv", header = T)

head(mirna.de)
head(genes)
head(targets)

mm = merge(targets,mirna.de,by.x = "Seq1", by.y = "miRNA")
dim(mm)
head(mm)
mm2 = merge(mm,genes,by.x = "Ensembl_Gene", by.y = "Gene.stable.ID")
dim(mm2)
head(mm2)

mm3 = mm2[which(mm2$padj.y < 0.01),]
head(mm3)
dim(mm3)

miup.genedown = mm3[which(mm3$log2FoldChange.x > 0 & mm3$log2FoldChange.y < 0),]
dim(miup.genedown)
midown.geneup = mm3[which(mm3$log2FoldChange.x < 0 & mm3$log2FoldChange.y > 0),]
dim(midown.geneup)
```


```{r}
utr.annotated = read.csv(file = "GenesWithAnnontated3UTR.csv", header = F)
head(utr.annotated)

combo = merge(genes,utr.annotated,by.x = "Gene.stable.ID", by.y = "V1")
dim(combo)
head(combo)

write.csv(combo,file="genes_with_annotated3UTR.csv")

```
