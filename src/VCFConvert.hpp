/*
 * VCFConvert.hpp
 *
 *  Created on: Nov 8, 2022
 *      Author: cjustin
 *
 *  Original intent: Take in a multiVCF file, reference and outputs sequence to be counted or
 *  Rounds up coverage to value specified, even in the case of multiple counts
 *  If rsIDs exist convert to counts
 */

#ifndef SRC_VCFCONVERT_HPP_
#define SRC_VCFCONVERT_HPP_

#include <vector>
#include <string>
#include <zlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <sstream>

#include "MultiCount.hpp"

#include "vendor/tsl/robin_map.h"
#include "vendor/tsl/robin_set.h"
#include "vendor/kseq.h"
#include "vendor/kseq_util.h"
#include "vendor/KseqHashIterator.hpp"

#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class VCFConvert {
public:
	VCFConvert() {
		//load in reference
		if(opt::verbose > 1){
			cerr << "Loading Reference " << opt::ref << endl;
		}
		gzFile fp;
		fp = gzopen(opt::ref.c_str(), "r");
		//read in seq
		kseq_t *seq = kseq_init(fp);
		int l = kseq_read(seq);
		while (l >= 0) {
			m_chrIDs[seq->name.s] = m_ref.size();
			m_ref.push_back(kseq_t());
			cpy_kseq(&m_ref.back(), seq);
			l = kseq_read(seq);
		}
		kseq_destroy(seq);
		gzclose(fp);
	}

	//loads in reference file
	void count(string filename) {

		if(opt::verbose > 1){
			cerr << "Reading VCF file: " << filename << endl;
		}

		//load in vcf file
		ifstream fh(filename.c_str());
		{
			string line;
			while (getline(fh, line)) {
				if (line.at(0) == '#') {
					std::stringstream ss;
					ss.str(line);
					string item;
					getline(ss, item, '\t');
					//parse out header and IDs
					//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	...
					if (item.compare("#CHROM") == 0) {
						//skip first 9 entries
						for (unsigned i = 0; i < 8; ++i) {
							getline(ss, item, '\t');
						}
						while (getline(ss, item, '\t')) {
							m_sampleIDs.push_back(item);
						}
						break;
					}
				}
			}
		}

		if(opt::verbose > 1){
			cerr << "Starting multicount of each rsID for " << m_sampleIDs.size() << " samples." << endl;
		}

		MultiCount counts(m_sampleIDs);

		enum Genotype {hom1, het, hom2};

#pragma omp parallel
		while (fh.good()) {
			std::stringstream ss;
			string line;
#pragma omp critical(good)
			getline(fh, line);

			if (fh.good()) {
				ss.str(line);
				string item;
				getline(ss, item, '\t');
				string chr = item;
				getline(ss, item, '\t');
				size_t loc = stoi(item);
				getline(ss, item, '\t');
				string rsID = item;
				if (opt::verbose > 2) {
#pragma omp critical
					cerr << "Processing site: " << rsID << endl;
				}
				getline(ss, item, '\t'); //skip
				if (item == ".") { //if vcf entry isn't a snp
					continue;
				}
				getline(ss, item, '\t');
				if (item.size() != 1) { //if vcf entry isn't a snp
					continue;
				}
				char var = item.at(0);
				pair<string, string> seqs = getSeqFromSite(chr, loc, var);

				//skip 4 columns
				for (unsigned i = 0; i < 4; ++i) {
					getline(ss, item, '\t');
				}

				unsigned sampleIndex = 0;
				vector<Genotype> genotypes(m_sampleIDs.size());
				while (getline(ss, item, '\t')) {
					if (item == "0|0") {
						genotypes[sampleIndex] = hom1;
					} else if (item == "0|1" || item == "1|0") {
						genotypes[sampleIndex] = het;
					} else if (item == "1|1") {
						genotypes[sampleIndex] = hom2;
					}
					sampleIndex++;
				}
				assert(sampleIndex == m_sampleIDs.size());

				for (KseqHashIterator itr(seqs.first.c_str(), seqs.first.size(),
						opt::k); itr != itr.end(); ++itr) {
					for (unsigned i = 0; i < sampleIndex; ++i) {
						if (genotypes[i] == hom1) {
							counts.insertCount(i, *itr, opt::multi * 2);
						} else if (genotypes[i] == het) {
							counts.insertCount(i, *itr, opt::multi);
						}
					}
				}
				for (KseqHashIterator itr(seqs.second.c_str(),
						seqs.second.size(), opt::k); itr != itr.end(); ++itr) {
					for (unsigned i = 0; i < sampleIndex; ++i) {
						if (genotypes[i] == hom2) {
							counts.insertCount(i, *itr, opt::multi * 2);
						} else if (genotypes[i] == het) {
							counts.insertCount(i, *itr, opt::multi);
						}
					}
				}
			}
		}
		fh.close();

		if(opt::verbose > 1){
			cerr << "Outputting counts" << endl;
		}

#pragma omp parallel for
		for(unsigned i = 0; i < m_sampleIDs.size(); ++i){
			if(opt::verbose > 1){
#pragma omp critical(cerr)
				cerr << "Outputting counts for " << m_sampleIDs.at(i) << endl;
			}
			ofstream out(m_sampleIDs.at(i) + ".counts.txt");
			counts.printCountsMax(i, out);
			out.close();
		}
		if (!opt::pca.empty()) {
			if(opt::verbose > 1){
				cerr << "Outputting matrix and normalization values for PCA" << endl;
			}
			ofstream out(opt::pca + "_matrix.tsv");
			ofstream centerFile(opt::pca + "_center.txt");
			counts.printNormMatrix(out, centerFile);
		}
	}

	~VCFConvert() {
		// TODO Auto-generated destructor stub
	}
private:
	vector<kseq_t> m_ref;
	vector<string> m_sampleIDs;
	tsl::robin_map<string,unsigned> m_chrIDs;

	pair<string,string> getSeqFromSite(string chr, size_t pos, char var){
		unsigned chrIndex = m_chrIDs.at(chr);
		char refStr[opt::window + 1];
		char varStr[opt::window + 1];
        size_t offset = pos - opt::window/2 - 1;
		strncpy(refStr, &(m_ref.at(chrIndex).seq.s[offset]), opt::window );
		strncpy(varStr, &(m_ref.at(chrIndex).seq.s[offset]), opt::window );
		varStr[opt::window/2] = var;
		varStr[opt::window] = '\0';
		refStr[opt::window] = '\0';
		return(std::make_pair(string(refStr), string(varStr)));
	}
};
#endif /* SRC_VCFCONVERT_HPP_ */
