# pmex-wild-microRNA
microRNA expression variation in a natural system sheds light into the molecular targets of hydrogen sulfide

## Run prost


## miRanda

Manual at: http://cbio.mskcc.org/microrna_data/manual.html

```
cd /data/kelley/projects/programs
ls -l miRanda-aug2010.tar.gz 
tar -xzf miRanda-aug2010.tar.gz 
cd miRanda-3.3a/
vi INSTALL 
./configure 
make
make install
cd src/
./miranda 
```

Download 5' and 3' UTR from Ensembl BioMart
Extract all genes that do not have annotated UTR with split.pl script (from https://www.biostars.org/p/127842/)
```
perl script.pl 3primeUTR.fa > 3primeUTR-sequencesOnly.fa
```

Extract all DE miRNAs from fasta of all mature miRNAs
DE-microRNA-p01.txt is from the R code or TableS3-miRNA Expression.csv (see Pmex-miRNAs.md)
```
seqtk subseq AdditionalFile2_pme_matures.fa DE-microRNA-p01.txt > DE-microRNA-p01.fa
miranda DE-microRNA-p01.fa 3primeUTR-sequencesOnly.fa -out pmex-miranda-3primeoutput-MFdesign-DE-en20.txt -en -20
grep ">>" pmex-miranda-3primeoutput-MFdesign-DE-en20.txt > pmex-miranda-3primeoutput-MFdesign-DE-hitsOnly.txt

```
Turn pmex-miranda-3primeoutput-MFdesign-DE-hitsOnly.txt into pmex-miranda-3prime-predictedTargets.csv for ease of use in R

See Pmex-miRNAs.md for R analyses

## Supplementary files for manuscript 
 - AdditionalFile1_pme_hairpins.fa
 - AdditionalFile2_pme_matures.fa
 - TableS2-Pmex miRNA Annotation.xlsx
 - TableS3-miRNA Expression.csv
 - TableS4-mRNA Expression.csv
 - TableS5-GO-Enrichment.xlsx

## Hand-edited gff for mitochondrial annotation compatibility with stringtie
 - GCF_001443325.1_P_mexicana-1.0_genomic_with_mito_sequence-edited.gff.gz


## mRNA analysis
 - Pmex-Taco-DGE-FINAL.Rmd
 - PmexGeneNameMatching.csv
 - PmexRNA_samples_2010.xlsx
 - gene_count_matrix-2019-10-01_no_STRG.csv
 
## miRNA analysis
 - Pmex-miRNAs.md
 - GenesWithAnnontated3UTR.csv
 - pmex-miranda-3prime-predictedTargets.csv
 - split.pl
