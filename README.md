# Sequence Data Similarity Detection / Nucleotide Fingerprinter
##Summary

This tools counts the number of specific k-mers within sequence data. The counts can then be compare to other counts to determine to compute the probability that sample are of the same origin to discover incongruent samples or sample swaps.

##Dependencies

* Python (Tested on 3.8.5)
* pyfaidx python module
* GCC (Tested on 9.3.0)
* zlibdev
* Autotools (if directly cloning from repo)

##Installation

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

##Usage

#####Generating k-mers from fasta file:

Given a VCF file and a reference genome you can produce fasta files with k-mers that one can use to create a fingerprinting. We have provided a set for human data based on similar criterion found in SNP microarrays.

Example:

```bash
python python scripts/extractSNPsfromVCF.py -p prefix_31 -v snps.vcf -f reference.fa -k 31
```

Creates 2 fasta files ( with the 31-mers specified in the VCF file. All non C/G <-> A/T conversions are ignored.

#####Counting the k-mers:

Using these set of k-mers we can then count all of these k-mers within a fastq file. Files may be gziped and multiple threads can be used.

Example:

```bash
ntfp -k 19 -t 2 -r prefix_31_AT.fa -a prefix_31_GC.fa sample_part1.fq sample_part2.fq > K19_31_counts.txt
```

Creates count file using 2 threads. A sliding window using 19-mers is used in this case and the highest count in the window is recorded.

You can specify the counts of each k-mer within

#####Evaluating the samples:

Example:

```bash
ntfpEval K19_31_counts1.txt K19_31_counts2.txt K19_31_counts3.txt > summary.tsv
```

A tsv file is produced that contains all combinations and p-values denoting the chance that the samples come from a different origin.


