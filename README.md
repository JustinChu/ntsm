# ntsm - Nucleotide Sequence/Sample Matcher
## Summary
Publication: [https://doi.org/10.1101/2023.11.01.565041](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giae024/7687245?login=false)

This tool counts the number of specific k-mers within sequence data. The counts can then be compared to other counts to determine and compute the probability that samples are of the same origin to discover incongruent samples or sample swaps. It is intended to be run before any analysis and can provide additional QC information like sequencing error rate.

By default, this tool will only return pairs of samples with the same origin. This tool can theoretically also be used for relatedness inference using the `-a` parameter, however in this case the PCA-based heuristic should not be used (see below).

## Dependencies

General:
* GCC (Tested on 9.3.0)
* zlibdev
* Autotools (if directly cloning from repo)

For generating site fasta files given a VCF file
* Python (Tested with 3.8.5)
* pyfaidx python module
* scikit-learn python module
* Perl (Tested with v5.26.2)
* bwa (Tested with 0.7.17)

## Installation

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

Given a VCF file and a reference genome, you can produce fasta files with k-mers that one can use to create fingerprinting. We have provided a set of human data based on similar criteria found in SNP microarrays.

The VCF file in this stage can be a single sample VCF, it just needs the variants.

Example:

```bash
ntsmSiteGen generate-sites name=prefix ref=reference.fa vcf=snps.vcf
```

Creates site fasta files referred to as `prefix_n{min missing sub k-mers}.fa` which is the set of k-mers used by the tool. By default all non C/G <-> A/T conversions are ignored.

Parameters:

```
w=31 #window size to consider sequences in this region
k=19 #kmer size used in the window region
t=4 #threads for any subprocess or tools
```

The sites fasta file generated is a set of k-mers for each site. `n` for each fasta file generated refers to the number of missing k-mers allowed after filtering out repetitive k-mers (e.g. `n=0` would mean you need each site to have all possible k-mers). The max value of `n` is `w - k` (so with default parameters 12). Small values of `n` will help you handle sequencing errors, but you need to have enough sites for our statistical analysis to work. In our testing retaining around ~10^5 sites seemed to work well, though if you end up with slightly less or a lot more ntsm should still function well.

If you do not wish to select your own sites, we currently include `data/human_sites_n10.fa` a fasta file with selected 96287 sites adequate for sample swap detection for human samples in the `data` folder.

##### Generating PCA from multiVCF file:

Once fasta files for sites have been created, it is possible to create a PCA rotation matrix for speeding up the analysis. To do so you must supply a multiVCF file from which the PCA will be built. This multi-sample VCF ideally should not contain the same samples as the VCF used in the sample swap detection process. It should be a set of reliable samples on which a PCA and rotational matrix would be based on. We note that the use of a rotational matrix is optional.

Example (can be run in same directory as `ntsmSiteGen generate-sites` after picking a sites file to use:

```{bash}
ntsmSiteGen generate-pca-rot-mat name=prefix ref=reference.fa multivcf=/mnt/1886F90E86F8ED5EmultiVCF.vcf sites=prefix_n10.fa
```

Again if you are working with human samples and do not wish to generate your own, we currently include `data/human_sites_rotationMat.tsv` and `human_sites_center.txt` to use in our PCA-based heurstic. We based our PCA and rotation matrix on 3202 samples from the 1000 Genomes Project.

##### Counting the k-mers:

Using this set of k-mers we can then count all of these k-mers within a fastq file. Files may be gzipped and multiple threads can be used. Each sample needs a separate run of this command and its own count files.

Example:

```bash
ntsmCount -t 2 -s sites.fa sample_part1.fq sample_part2.fq > counts.txt
```

Creates count file using 2 threads. A sliding window using 19-mers is used in this case and the highest count in the window is recorded.

If your files are unsorted and have massive coverage, you may also intentionally run less reads using the `-m` parameter:

```bash
ntsmCount -t 2 -m 10 -s sites.fa sample_part1.fq sample_part2.fq > counts.txt
```

This will run the file until the average site coverage reaches 10x, which should be adequate for most sequencing data types. Lower coverage is possible if the read error rate is low enough.

Output Example:

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

Header lines (`#@`) helps in error rate estimation.

##### Evaluating the samples:

Example command:

```bash
ntsmEval HG002_rep1_counts.txt HG002_rep2_counts.txt HG003_counts.txt HG004_counts.txt > summary.tsv
```

or if you wish optionally to speed up the analysis using a PCA rotation matrix:

```bash
ntsmEval -a -t 16 -n data/human_sites_center.txt -p data/human_sites_rotationMat.tsv HG002_rep1_counts.txt HG002_rep2_counts.txt HG003_counts.txt HG004_counts.txt > summary.tsv
```

Output Example (with `-a` option):
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
Note that the `-a` parameter will not output all pairwise comparisons if the PCA heuristic is used and will only consider comparing samples within a certain distance in PCA space. Though related samples are more likely to exist closer in PCA space, if using ntsm for relatedness inference you may miss some related pairs if you use this heuristic. You would not likely want to use `-a` with the PCA heuristic in most cases and the example above is for more illustrative purposes of the tool's use and output.

Column explanations:
* sampleX: Filename for sample X
* score: Log-likelihood-based score to determine if samples are the same or differ
* same: 1 means the tool thinks the sample is the same and 0 is if the tool thinks they differ
* dist: Distance of sample in PCA space
* relate: Relatedness determined via shared heterozygous sites
* relate: Relatedness determined via shared heterozygous sites
* ibs0: Number of sites with alleles not shared between two samples
* ibs2: Number of sites with the same genotype between two samples
* homConcord: Homozygous concordance determined via shared homozygous sites
* hetX: Number of heterozygous sites for sample X for all sites considered in comparison
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

This provides generic QC information useful without other samples (e.g. error rate) and if inclined a means of plotting the sample relative to others on a PCA plot. The number of columns is dependant on the number of principal components used.

