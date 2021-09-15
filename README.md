# Sequence Data Similarity Detection / Nucleotide Fingerprinter
## Summary

This tools counts the number of specific k-mers within sequence data. The counts can then be compare to other counts to determine to compute the probability that sample are of the same origin to discover incongruent samples or sample swaps.

## Dependencies

* Python (Tested on 3.8.5)
* pyfaidx python module
* GCC (Tested on 9.3.0)
* zlibdev
* Autotools (if directly cloning from repo)

## Installation

If cloning directly from the repository make sure you get the required submodules:
```bash
git submodule update --init
```
If cloning directly from the repository run:
```bash
./autogen.sh
```
Compiling should be as easy as:
```bash
./configure && make
```
To install in a specified directory:
```bash
./configure --prefix=/PATH && make install
```

## Usage

##### Generating k-mers from fasta file:

Given a VCF file and a reference genome you can produce fasta files with k-mers that one can use to create a fingerprinting. We have provided a set for human data based on similar criterion found in SNP microarrays.

Example:

```bash
python scripts/extractSNPsfromVCF.py -p prefix_31 -v snps.vcf -f reference.fa -k 31
```

Creates 2 fasta files ( with the 31-mers specified in the VCF file. All non C/G <-> A/T conversions are ignored.

##### Counting the k-mers:

Using these set of k-mers we can then count all of these k-mers within a fastq file. Files may be gziped and multiple threads can be used.

Example:

```bash
ntfp -k 19 -t 2 -r prefix_31_AT.fa -a prefix_31_GC.fa <(pigz -cd sample_part1.fq) <(pigz -cd sample_part2.fq) > K19_31_counts.txt
```

Creates count file using 2 threads. A sliding window using 19-mers is used in this case and the highest count in the window is recorded.

Outout Example:

```
rs16824588	26	0
rs1181883	80	12
rs200458	27	29
rs228648	16	14
rs848209	33	0
rs11203366	7	16
rs2240335	25	0
rs1541185	21	18
rs2254358	0	35
rs6678540	12	15
rs2294228	0	29
```


##### Evaluating the samples:

Example:

```bash
ntfpEval K19_31_counts1.txt K19_31_counts2.txt K19_31_counts3.txt > summary.tsv
```

A tsv file is produced that contains all combinations of geometric mean of Fisher's exact tests on all k-mer sites (p-values).

Output Example:
```
K19_31_HG002_CCS1.txt	K19_31_HG002_CCS2.txt	0.541319
K19_31_HG002_CCS1.txt	K19_31_HG002_HC1.txt	0.162735
K19_31_HG002_CCS1.txt	K19_31_HG002_HC2.txt	0.342228
K19_31_HG002_CCS1.txt	K19_31_HG002_HC3.txt	0.288665
K19_31_HG002_CCS1.txt	K19_31_HG002_HC4.txt	0.473151
K19_31_HG002_CCS1.txt	K19_31_HG002_NP.txt	0.134269
K19_31_HG002_CCS1.txt	K19_31_HG002_SR1.txt	0.503195
K19_31_HG002_CCS1.txt	K19_31_HG002_SR2.txt	0.478716
K19_31_HG002_CCS1.txt	K19_31_NA12878_CCS1.txt	0.000272751
K19_31_HG002_CCS1.txt	K19_31_NA12878_HC1.txt	9.68865e-05
K19_31_HG002_CCS1.txt	K19_31_NA12878_SR1.txt	5.47665e-05
K19_31_HG002_CCS1.txt	K19_31_NA12878_SR2.txt	0.000114591
K19_31_HG002_CCS1.txt	K19_31_NA12878_SR3.txt	0.000190844
K19_31_HG002_CCS1.txt	K19_31_NA12878_SR4.txt	0.000222317
K19_31_HG002_CCS1.txt	K19_31_NA12878_SR5.txt	0.000468057
K19_31_HG002_CCS1.txt	K19_31_NA12878_SR6.txt	0.000474764
K19_31_HG002_CCS1.txt	K19_31_NA12878_SR7.txt	0.000582513
K19_31_HG002_CCS1.txt	K19_31_NA12878_SR8.txt	0.000174436
```


