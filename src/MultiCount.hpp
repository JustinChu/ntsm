/*
 * MultiCount.hpp
 *
 *  Created on: Nov 22, 2022
 *      Author: cjustin
 *
 *      TODO: Shares function to FingerPrint.hpp, possible to refactor and integrate
 */

#ifndef SRC_MULTICOUNT_HPP_
#define SRC_MULTICOUNT_HPP_
#include <vector>
#include <string>
#include <omp.h>
#include <stdio.h>
#include <zlib.h>
#include <iostream>

#include "Options.h"

#include "vendor/tsl/robin_map.h"
#include "vendor/tsl/robin_set.h"
#include "vendor/KseqHashIterator.hpp"
#include "vendor/kseq_util.h"

#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "vendor/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class MultiCount {
public:

	typedef uint16_t SampleID;
	typedef uint64_t HashedKmer;
	static constexpr double UNDEF = numeric_limits<double>::max();

	MultiCount(const vector<string> &sampleIDs) :
			m_sampleIDs(sampleIDs) {
//	MultiCount(size_t size){
		//read in fasta files
		//generate hash table
		initCountsHash(sampleIDs.size());
	}

	void insertCount(unsigned index, const char *seqs, uint64_t seql, unsigned multi = 1) {
		for (KseqHashIterator itr(seqs, seql, opt::k); itr != itr.end();
				++itr) {
			if (m_kmerToHash.find(*itr) != m_kmerToHash.end()) {
//				if ((*m_counts[m_kmerToHash.at(*itr)]).at(index) > 0) {
//					cerr << "Warning double count detected: " << index << " "
//							<< to_string(
//									(*m_counts[m_kmerToHash.at(*itr)]).at(
//											index)) << " " << itr.getPos()
//							<< " " << seqs << endl;
//				}
#pragma omp atomic update
				(*m_counts[m_kmerToHash.at(*itr)])[index] += multi;
			}
		}
	}

	void printCountsMax(unsigned index, ostream &out = cout) const{
		out << "\n#locusID\tcountAT\tcountCG\tsumAT\tsumCG\tdistinctAT\tdistinctCG\n";
		string outStr = "";
		for (size_t i = 0; i < m_alleleIDs.size() ; ++i) {
			outStr.clear();
			const vector<uint64_t> &allele1 = *m_alleleIDToKmerRef.at(i);
			const vector<uint64_t> &allele2 = *m_alleleIDToKmerVar.at(i);
			unsigned maxCountREF = 0;
			unsigned maxCountVAR = 0;
			unsigned countSumAT = 0;
			unsigned countSumCG = 0;
			for(size_t j = 0; j < allele1.size() ; ++j) {
				unsigned freqAlle = m_counts.at(m_kmerToHash.at(allele1.at(j)))->at(index);
				if(maxCountREF < freqAlle){
					maxCountREF = freqAlle;
				}
				countSumAT += freqAlle;
			}
			for(size_t j = 0; j < allele2.size() ; ++j) {
				unsigned freqAlle = m_counts.at(m_kmerToHash.at(allele2.at(j)))->at(index);
				if(maxCountVAR < freqAlle){
					maxCountVAR = freqAlle;
				}
				countSumCG += freqAlle;
			}
			outStr += m_alleleIDs.at(i);
			outStr += "\t";
			outStr += std::to_string(maxCountREF);
			outStr += "\t";
			outStr += std::to_string(maxCountVAR);
			outStr += "\t";
			outStr += std::to_string(countSumAT);
			outStr += "\t";
			outStr += std::to_string(countSumCG);
			outStr += "\t";
			outStr += std::to_string(allele1.size());
			outStr += "\t";
			outStr += std::to_string(allele2.size());
			outStr += "\n";
			out << outStr;
		}
	}

	/*
	 * Prints centered and normalized matrix intended for PCA analysis
	 * Converts missing values to average value
	 * Missing values are not use in calculation
	 */
	void printNormMatrix(ostream &out, ostream &centerFile) {
		out << "alleleID";
		//print allele ids
		for (unsigned i = 0; i < m_sampleIDs.size(); ++i) {
			out << "\t" << m_sampleIDs[i];
		}
		out << endl;

		for (unsigned i = 0; i < m_alleleIDs.size(); ++i) {
			const vector<uint64_t> &allele1 = *m_alleleIDToKmerRef.at(i);
			const vector<uint64_t> &allele2 = *m_alleleIDToKmerVar.at(i);
			vector<double> values(m_sampleIDs.size());
			double sum = 0;
			unsigned size = 0;
			for (size_t j = 0; j < m_sampleIDs.size(); ++j) {
				unsigned maxCountREF = 0;
				unsigned maxCountVAR = 0;
				for (size_t k = 0; k < allele1.size(); ++k) {
					unsigned freqAlle = m_counts.at(
							m_kmerToHash.at(allele1.at(k)))->at(j);
					if (maxCountREF < freqAlle) {
						maxCountREF = freqAlle;
					}
				}
				for (size_t k = 0; k < allele2.size(); ++k) {
					unsigned freqAlle = m_counts.at(
							m_kmerToHash.at(allele2.at(k)))->at(j);
					if (maxCountVAR < freqAlle) {
						maxCountVAR = freqAlle;
					}
				}
				unsigned denom = maxCountREF + maxCountVAR;
				if (denom == 0) {
					values[j] = UNDEF;
				} else {
					values[j] = double(maxCountREF) / double(denom);
					++size;
					sum += values.at(j);
				}
			}
			out << m_alleleIDs.at(i);
			double center = sum / double(size);
			for (unsigned j = 0; j < values.size(); ++j) {
				if (values.at(j) == UNDEF) {
					out << "\t0";
				} else {
					out << "\t" << to_string(values.at(j) - center);
				}
			}
//			centerFile << m_alleleIDs[i];
			centerFile << center << endl;
			out << endl;
		}
	}

private:
	tsl::robin_map<HashedKmer, SampleID> m_kmerToHash; //hashvalue->sample
	vector<unique_ptr<vector<uint8_t>>> m_counts; //sample->counts
	vector<unique_ptr<vector<HashedKmer>>> m_alleleIDToKmerRef; //alleleID_ref->hashvalue
	vector<unique_ptr<vector<HashedKmer>>> m_alleleIDToKmerVar; //alleleID_var->hashvalue
	vector<string> m_alleleIDs;
	const vector<string> &m_sampleIDs;

	void initCountsHash(size_t size){
		gzFile fp = gzopen(opt::snp.c_str(), "r");
		tsl::robin_set<uint64_t> dupes;
		if (fp == Z_NULL) {
			std::cerr << "file " << opt::snp.c_str() << " cannot be opened"
					<< std::endl;
			exit(1);
		} else if (opt::verbose) {
			std::cerr << "Opening " << opt::snp.c_str() << std::endl;
		}
		{
			kseq_t *seq = kseq_init(fp);
			int l = kseq_read(seq);
			unsigned entryNum = 0;
			unsigned totalCount = 0;
			while (l >= 0) {
				if (entryNum % 2 == 0) {
					unsigned index = entryNum / 2;
					assert(index == m_alleleIDToKmerRef.size());
					m_alleleIDToKmerRef.emplace_back(
							unique_ptr<vector<uint64_t>>(
									new vector<uint64_t>()));
					//k-merize and
					for (KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k);
							itr != itr.end(); ++itr) {
						uint64_t hv = *itr;
						//check for duplicates
						if (m_kmerToHash.find(hv) != m_kmerToHash.end()) {
							cerr << "Warning: " << seq->name.s
									<< " of REF file has a k-mer collision at pos: "
									<< itr.getPos() << endl;
							dupes.insert(hv);
						} else {
							m_alleleIDToKmerRef[index]->emplace_back(hv);
							m_kmerToHash[hv] = totalCount;
							m_counts.emplace_back(
									unique_ptr<vector<uint8_t>>(
											new vector<uint8_t>(size, 0)));
							totalCount++;
						}
					}
					m_alleleIDs.emplace_back(seq->name.s);
				} else {
					unsigned index = entryNum / 2;
					assert(index == m_alleleIDToKmerVar.size());
					m_alleleIDToKmerVar.emplace_back(
							unique_ptr<vector<uint64_t>>(
									new vector<uint64_t>()));
					//k-merize and insert
					for (KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k);
							itr != itr.end(); ++itr) {
						uint64_t hv = *itr;
						//check for duplicates
						if (m_kmerToHash.find(hv) != m_kmerToHash.end()) {
							cerr << "Warning: " << seq->name.s
									<< " of VAR file has a k-mer collision at pos: "
									<< itr.getPos() << endl;
							dupes.insert(hv);
						} else {
							m_alleleIDToKmerVar[index]->emplace_back(hv);
							m_kmerToHash[hv] = totalCount;
							m_counts.emplace_back(
									unique_ptr<vector<uint8_t>>(
											new vector<uint8_t>(size, 0)));
							totalCount++;
						}
					}
				}
				l = kseq_read(seq);
				entryNum++;
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
		//TODO remove total count indexes?
		if (!opt::dupes) {
			//remove dupes
			for (tsl::robin_set<uint64_t>::iterator itr = dupes.begin();
					itr != dupes.end(); ++itr) {
				m_kmerToHash.erase(*itr);
			}
		}
	}
};

#endif /* SRC_MULTICOUNT_HPP_ */
