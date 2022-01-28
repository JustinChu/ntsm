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
#include "vendor/KseqHashIterator.hpp"
#include "vendor/kseq.h"
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

//TODO current implementation may be incorrect due to hash collisions

using namespace std;

class FingerPrint {
public:

	typedef uint16_t alleleID;

	FingerPrint(const vector<string> &filenames) :
			m_filenames(filenames), m_totalCounts(0), m_maxCounts(0), m_totalBases(0) {
		//read in fasta files
		//generate hash table
		initCountsHash();
		if(opt::covThresh != 0){
			m_maxCounts = (m_counts.size() * opt::covThresh)/2;
		}
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
				for (KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k);
						itr != itr.end(); ++itr) {
					if(m_counts.find(*itr) != m_counts.end()){
#pragma omp atomic update
						++m_counts[*itr];
#pragma omp atomic update
						++m_totalCounts;
					}
				}
#pragma omp atomic update
				m_totalBases += seq->seq.l;
				l = kseq_read(seq);
				index++;
				//terminate early
				if (m_maxCounts != 0 && m_totalCounts > m_maxCounts) {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
	}

	void printCounts(){
		string tempStr = "";
		for (size_t i = 0; i < m_alleleIDs.size() ; ++i) {
			const vector<pair<uint64_t, uint64_t>> &allelePair = *m_alleleIDToKmer[i];
			size_t maxCountWT = 0;
			size_t maxCountVAR = 0;
			for(size_t j = 0; j < allelePair.size() ; ++j) {
				double freqAlle1 = m_counts[allelePair.at(j).first];
				double freqAlle2 = m_counts[allelePair.at(j).second];
				if(maxCountWT < freqAlle1){
					maxCountWT = freqAlle1;
				}
				if(maxCountVAR < freqAlle2){
					maxCountVAR = freqAlle2;
				}
			}
			cout << m_alleleIDs.at(i) << "\t" << maxCountWT << "\t" << maxCountVAR << endl;
		}
		cout << endl;
	}

	uint64_t getTotalKmerCounts(){
		return m_totalCounts;
	}

	uint64_t getTotalCounts(){
		return m_totalBases;
	}

//	/*
//	 * 0 = missing alleles
//	 * 1 = homozygous wt
//	 * 2 = homozygous var
//	 * 3 = heterzygous wt/var
//	 */
//	void computeFingerPrint(){
//		//assume count is not propotional to k-mer count
//		double expectedCoverage = m_totalCounts;
//		for (vector<pair<uint64_t, uint64_t>>::iterator itr =
//				m_allelePairs.begin(); itr != m_allelePairs.end(); ++itr) {
//			double freqAlle1 = m_counts[itr->first];
//			double freqAlle2 = m_counts[itr->second];
//			bool wt, var = false;
//
//			if(freqAlle1 >= expectedCoverage * opt::minProp){
//				wt = true;
//			}
//			if (freqAlle2 >= expectedCoverage * opt::minProp) {
//				var = true;
//			}
//			if(!wt && !var){
//				cout << "0";
//			}
//			else if(wt && var){
//				cout << "3";
//			}
//			else if(wt){
//				cout << "1";
//			}
//			else if(var){
//				cout << "2";
//			}
//		}
//		cout << endl;
//	}

	virtual ~FingerPrint() {
	}
private:
	const vector<string> &m_filenames;
	size_t m_totalCounts;
	size_t m_maxCounts;
	tsl::robin_map<uint64_t, size_t> m_counts;
	tsl::robin_map<alleleID, shared_ptr<vector<pair<uint64_t, uint64_t>>>> m_alleleIDToKmer;
//	vector<pair<uint64_t, uint64_t>> m_allelePairs;
	vector<string> m_alleleIDs;
	uint64_t m_totalBases;

	static const unsigned interval = 65536; //interval to be used when considering

	void initCountsHash(){
		gzFile fp1, fp2;
		fp1 = gzopen(opt::ref.c_str(), "r");
		fp2 = gzopen(opt::var.c_str(), "r");
		if (fp1 == Z_NULL) {
			std::cerr << "file " << opt::ref.c_str() << " cannot be opened"
					<< std::endl;
			exit(1);
		} else if (opt::verbose) {
			std::cerr << "Opening " << opt::ref.c_str() << std::endl;
		}
		if (fp2 == Z_NULL) {
			std::cerr << "file " << opt::var.c_str() << " cannot be opened"
					<< std::endl;
			exit(1);
		} else if (opt::verbose) {
			std::cerr << "Opening " << opt::var.c_str() << std::endl;
		}
		kseq_t *seq1 = kseq_init(fp1);
		kseq_t *seq2 = kseq_init(fp2);
		int l1 = kseq_read(seq1);
		int l2 = kseq_read(seq2);
		unsigned index = 0;
		while (l1 >= 0 && l2) {
			KseqHashIterator itr1(seq1->seq.s, seq1->seq.l, opt::k);
			KseqHashIterator itr2(seq2->seq.s, seq2->seq.l, opt::k);
			m_alleleIDToKmer[m_alleleIDs.size()] = shared_ptr<
					vector<pair<uint64_t, uint64_t>>>(
					new vector<pair<uint64_t, uint64_t>>());
			m_alleleIDToKmer[m_alleleIDs.size()]->reserve(seq1->seq.l - opt::k + 1);

			//k-merize and insert
			//TODO Add some more file and sanity checks
			for (;itr1 != itr1.end() && itr2 != itr2.end(); ++itr1, ++itr2) {
				uint64_t hv1 = *itr1;
				uint64_t hv2 = *itr2;
				//check for duplicates
				if (m_counts.find(hv1) != m_counts.end()
						|| m_counts.find(hv2) != m_counts.end()) {
					cerr << seq1->name.s << " has a k-mer collision" << endl;
				} else {
					m_alleleIDToKmer[m_alleleIDs.size()]->emplace_back(
							make_pair(hv1, hv2));
					m_counts[hv1] = 0;
					m_counts[hv2] = 0;
				}
			}
			m_alleleIDs.emplace_back(seq1->name.s);
			l1 = kseq_read(seq1);
			l2 = kseq_read(seq2);
			index++;
		}
		kseq_destroy(seq1);
		kseq_destroy(seq2);
		gzclose(fp1);
		gzclose(fp2);
	}
};
#endif /* SRC_FINGERPRINT_HPP_ */
