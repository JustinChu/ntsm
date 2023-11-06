# ntsm - Nucleotide Sequence/Sample Matcher
## Summary

This tools counts the number of specific k-mers within sequence data. The counts can then be compare to other counts to determine to compute the probability that sample are of the same origin to discover incongruent samples or sample swaps.

## Dependencies

General:
* GCC (Tested on 9.3.0)
* zlibdev
* Autotools (if directly cloning from repo)

For generating site fasta files given a VCF file
* Python (Tested with 3.8.5)
* pyfaidx python module
* Perl (Tested with v5.26.2)
* bwa (Tested with 0.7.17)

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
scripts/generateSites name=prefix ref=reference.fa vcf=snps.vcf
```

Creates a fasta files VCF file. All non C/G <-> A/T conversions are ignored.

Parameters:

```
w=31 #window size to consider sequences in this region
k=19 #kmer size used in window region
t=4 #threads for any subprocess or tools
n=0 #number of sub k-mers to allow
```

If you do not wish to select your own sites, we currently include `human_sites_n10.fa` a fasta selected sites with 96287 sites adequate for sample swap detection for human samples in the `data` folder.

##### Counting the k-mers:

Using these set of k-mers we can then count all of these k-mers within a fastq file. Files may be gziped and multiple threads can be used.

Example:

```bash
ntsmCount -t 2 -s sites.fa sample_part1.fq sample_part2.fq > counts.txt
```

Creates count file using 2 threads. A sliding window using 19-mers is used in this case and the highest count in the window is recorded.

Outout Example:

```
#@TK	119443488624
#@KS	19
#locusID	countAT	countCG	sumAT	sumCG	distinctAT	distinctCG
rs1741692	23	22	68	190	3	9
rs6419870	0	43	0	86	1	2
rs3171927	19	20	91	20	5	1
rs12057128	16	0	31	0	2	5
rs11976368	43	0	43	0	1	10
rs4545798	17	13	34	65	2	5
rs10888802	0	37	0	355	9	10
...
```

Header lines (`#@`) help in error rate estimation.

##### Evaluating the samples:

Example:

```bash
ntsmEval sampleA_counts.txt sampleB_counts.txt sampleC_counts.txt > summary.tsv
```

Output Example:

```
sample1	sample2	relate	ibs0	ibs2	homConcord	hets1	hets2	sharedHets	hom1	hom2	sharedHom	n	score	same	cov1	cov2	error_rate1	error_rate2
a1	a2	0.997193	0	209139	0.998521	72327	72463	72124	137354	137218	137015	209681	0.079416	1	37.275289	45.248168	0.004990	0.005695
a1	b	0.473059	902	136180	0.717727	72327	71862	35799	137346	137811	100381	209673	1.690037	0	37.275289	45.040483	0.004990	0.007278
a1	c	0.499938	57	136556	0.734409	72327	73271	36273	137338	136394	100283	209665	1.643607	0	37.275289	44.389736	0.004990	0.005253
a2	b	0.473880	899	136153	0.717900	72463	71862	35852	137210	137811	100301	209673	1.790874	0	45.248168	45.040483	0.005695	0.007278
a2	c	0.499814	56	136535	0.733852	72463	73271	36330	137202	136394	100205	209665	1.741786	0	45.248168	44.389736	0.005695	0.005253
b	c	-0.017423	14756	106307	0.355821	71860	73269	28260	137812	136403	78047	209672	3.430927	0	45.040483	44.389736	0.007278	0.005253
...
```

Column explainations:

* relate: Relatedness determined via shared heterozygous sites
* ibs0: Number of sites with alleles not shared between two samples
* ibs2: Number of sites with the same genotype between two samples
* homConcord: Homozygous concordance determined via shared homozygous sites
* hetsX: Number of heterzygous sites for sample X
* sharedHets: Number of shared hetero
* homX: Number of homozygous sites for sample X
* sharedHom: Number of shared homozygous sites
* n: number of unfiltered sites used in comparison
* score: log-likelihood based score to determine in samples are the same or differ
* same: 1 means the tool thinks the sample is the same and 0 is if the tools thinks they differ
* covX: coverage of sample X
* error_rateX: error rate of sample X. May underestimate error caused by long indels.


