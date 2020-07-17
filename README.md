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
```
seqtk subseq pmex-microRNAs.fa DE-microRNA-p01.txt > DE-microRNA-p01.fa
miranda DE-microRNA-p01.fa 3primeUTR-sequencesOnly.fa -out pmex-miranda-3primeoutput-MFdesign-DE-en20.txt -en -20
grep ">>" pmex-miranda-3primeoutput-MFdesign-DE-en20.txt > pmex-miranda-3primeoutput-MFdesign-DE-hitsOnly.txt

```
Turn pmex-miranda-3primeoutput-MFdesign-DE-hitsOnly.txt into pmex-miranda-3prime-predictedTargets.csv for ease of use in R
