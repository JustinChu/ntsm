# ntsm - Nucleotide Sequence/Sample Matcher
## Summary

This tools counts the number of specific k-mers within sequence data. The counts can then be compare to other counts to determine to compute the probability that sample are of the same origin to discover incongruent samples or sample swaps.

Intended to be run at before any analysis and can provide some additional QC information like sequencing error rate.

Manuscript preprint: https://doi.org/10.1101/2023.11.01.565041

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

If you do not wish to select your own sites, we currently include `data/human_sites_n10.fa` a fasta selected sites with 96287 sites adequate for sample swap detection for human samples in the `data` folder.

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

Example command:

```bash
ntsmEval HG002_rep1_counts.txt HG002_rep2_counts.txt HG003_counts.txt HG004_counts.txt > summary.tsv
```

or if you wish optionally to speed up the analysis using a PCA rotation matrix:

```bash
ntsmEval -a -t 16 -n data/human_sites_center.txt -p data/human_sites_rotationMat.tsv HG002_rep1_counts.txt HG002_rep2_counts.txt HG003_counts.txt HG004_counts.txt > summary.tsv
```

Output Example (with -a option):
```
sample1	sample2	score	same	dist	relate	ibs0	ibs2	homConcord	het1	het2	sharedHet	hom1	hom2	sharedHom	n	cov1	cov2	errorRate1	errorRate2	miss1	miss2	allHom1	allHom2	allHet1	allHet2
HG002_rep1_counts.txt	HG002_rep2_counts.txt	0.07988	1	0.004839	0.996827	0	95971	0.998287	33720	33787	33613	62532	62465	62358	96252	37.416162	45.260554	0.003493	0.004301	35	35	62532	62465	33720	33787
HG003_counts.txt	HG004_counts.txt	3.430842	0	7.512569	-0.003973	6649	48672	0.355549	33473	33781	13165	62772	62464	35507	96245	44.931787	44.068285	0.005968	0.004208	34	38	62779	62466	33474	33783
HG002_rep1_counts.txt	HG003_counts.txt	1.660803	0	4.732675	0.498775	24	62518	0.731288	33720	33474	16744	62528	62774	45774	96248	37.416162	44.931787	0.003493	0.005968	35	34	62532	62779	33720	33474
HG002_rep1_counts.txt	HG004_counts.txt	1.653478	0	2.872071	0.500089	19	62525	0.72982	33720	33783	16901	62525	62462	45624	96245	37.416162	44.068285	0.003493	0.004208	35	38	62532	62466	33720	33783
HG002_rep2_counts.txt	HG003_counts.txt	1.760081	0	4.707002	0.499821	24	62521	0.73156	33787	33474	16779	62461	62774	45742	96248	45.260554	44.931787	0.004301	0.005968	35	34	62465	62779	33787	33474
HG002_rep2_counts.txt	HG004_counts.txt	1.74858	0	2.78644	0.4996	19	62488	0.729034	33787	33783	16916	62458	62462	45572	96245	45.260554	44.068285	0.004301	0.004208	35	38	62465	62466	33787	33783
...
```

Column explainations:
* sampleX: Filename for sample X
* score: Log-likelihood based score to determine in samples are the same or differ
* same: 1 means the tool thinks the sample is the same and 0 is if the tools thinks they differ
* dist: Distance of sample in PCA space
* relate: Relatedness determined via shared heterozygous sites
* relate: Relatedness determined via shared heterozygous sites
* ibs0: Number of sites with alleles not shared between two samples
* ibs2: Number of sites with the same genotype between two samples
* homConcord: Homozygous concordance determined via shared homozygous sites
* hetX: Number of heterzygous sites for sample X for all sites considered in comparison
* sharedHets: Number of shared heterozygous sites
* homX: Number of homozygous sites for sample X for all sites considered in comparison
* sharedHom: Number of shared homozygous sites
* n: number of unfiltered sites used in comparison
* covX: Coverage of sample X
* errorRateX: Error rate of sample X. May underestimate error caused by long indels.
* missX: Total number of missing sites in sample X
* allHomX: Total number of homozygous sites in sample X
* allHetX: Total number of heterozygous sites in sample X

If run on a single counts file the output can look like this:
```
sample	cov	errorRate	miss	hom	het	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10	PC11	PC12	PC13	PC14	PC15	PC16	PC17	PC18	PC19	PC20
HG002_rep1_counts.txt	37.416162	0.003493	35	62532	33720	-14.254352	-23.285693	-0.373179	-7.315873	-1.187992	-5.494577	-0.434657	-0.334589	-0.832297	-1.160507	0.286102	-0.114464	1.013333	0.252766	-0.204102	0.465836	0.694361	0.099620	0.019345	-1.195279
```

This provides generic QC information useful without other samples (e.g. error rate) and if inclined a means of plotting the sample relative to others on a PCA plot. Number of columns is variable dependant on the number of principle components used.

