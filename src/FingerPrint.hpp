/*
 * FingerPrint.hpp
 *
 *  Created on: Jul 6, 2021
 *      Author: cjustin
 */

#ifndef SRC_FINGERPRINT_HPP_
#define SRC_FINGERPRINT_HPP_
#include <vector>
#include <string>
#include <omp.h>
#include <stdio.h>
#include <zlib.h>
#include <iostream>

#include "Options.h"

#include "vendor/tsl/robin_map.h"
#include "vendor/ntHash/ntHashIterator.hpp"
#include "vendor/kseq.h"
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

//TODO current implementation may be incorrect due to hash collisions

using namespace std;

class FingerPrint {
public:
	FingerPrint(const vector<string> &filenames) :
			m_filenames(filenames), m_totalCounts(0){
		//read in fasta files
		//generate hash table
		initCountsHash();
	}

	void computeCounts(){
#pragma omp parallel for
		for (unsigned i = 0; i < m_filenames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_filenames[i].c_str(), "r");
			if (fp == Z_NULL) {
				std::cerr << "file " << m_filenames[i] << " cannot be opened"
						<< std::endl;
				exit(1);
			} else if (opt::verbose) {
				std::cerr << "Opening " << m_filenames[i] << std::endl;
			}
			//read in seq
			kseq_t *seq = kseq_init(fp);
			int l = kseq_read(seq);
			unsigned index = 0;
			while (l >= 0) {
				//k-merize and insert
				for (ntHashIterator itr(seq->seq.s, 1, opt::k);
						itr != itr.end(); ++itr) {
					//if kmer exists inside set
					uint64_t hv = (*itr)[0];
					if(m_counts.find(hv) != m_counts.end()){
#pragma omp atomic update
						++m_counts[hv];
					}
#pragma omp atomic update
					++m_totalCounts;
				}
				l = kseq_read(seq);
				index++;
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
	}

	/*
	 * 0 = missing alleles
	 * 1 = homozygous wt
	 * 2 = homozygous var
	 * 3 = heterzygous wt/var
	 */
	void computeFingerPrint(){
		double expectedCoverage = opt::genomeSize/m_totalCounts;
		for (vector<pair<uint64_t, uint64_t>>::iterator itr =
				m_allelePairs.begin(); itr != m_allelePairs.end(); ++itr) {
			double freqAlle1 = m_counts[itr->first];
			double freqAlle2 = m_counts[itr->second];
			bool wt, var = false;
			if(freqAlle1 >= expectedCoverage * opt::minProp){
				wt = true;
			}
			if (freqAlle2 >= expectedCoverage * opt::minProp) {
				var = true;
			}
			if(!wt && !var){
				cout << "0";
			}
			else if(wt && var){
				cout << "3";
			}
			else if(wt){
				cout << "1";
			}
			else if(var){
				cout << "2";
			}
		}
		cout << endl;
	}

	virtual ~FingerPrint() {
	}
private:
	const vector<string> &m_filenames;
	size_t m_totalCounts;
	tsl::robin_map<uint64_t, size_t> m_counts;
	vector<pair<uint64_t, uint64_t>> m_allelePairs;

	void initCountsHash(){
		gzFile fp1, fp2;
		fp1 = gzopen(opt::ref.c_str(), "r");
		fp2 = gzopen(opt::var.c_str(), "r");
		validateFile(fp1);
		validateFile(fp2);
		kseq_t *seq1 = kseq_init(fp1);
		kseq_t *seq2 = kseq_init(fp2);
		int l1 = kseq_read(seq1);
		int l2 = kseq_read(seq2);
		unsigned index = 0;
		while (l1 >= 0 && l2) {
			ntHashIterator itr1(seq1->seq.s, 1, opt::k);
			ntHashIterator itr2(seq1->seq.s, 1, opt::k);
			//k-merize and insert
			//TODO Add some more file and sanity checks
			for (;itr1 != itr1.end() && itr2 != itr2.end(); ++itr1, ++itr2) {
				uint64_t hv1 = (*itr1)[0];
				uint64_t hv2 = (*itr2)[0];
				m_allelePairs.emplace_back(make_pair(hv1, hv2));
				assert(m_counts.find(hv1) == m_counts.end());
				assert(m_counts.find(hv2) == m_counts.end());
				m_counts[hv1] = 0;
				m_counts[hv2] = 0;
			}
			l1 = kseq_read(seq1);
			l2 = kseq_read(seq2);
			index++;
		}
		kseq_destroy(seq1);
		kseq_destroy(seq2);
		gzclose(fp1);
		gzclose(fp2);
	}

	void validateFile(gzFile &fp){
		if (fp == Z_NULL) {
			std::cerr << "file " << opt::ref.c_str() << " cannot be opened"
					<< std::endl;
			exit(1);
		} else if (opt::verbose) {
			std::cerr << "Opening " << opt::ref.c_str() << std::endl;
		}
	}

};
#endif /* SRC_FINGERPRINT_HPP_ */
