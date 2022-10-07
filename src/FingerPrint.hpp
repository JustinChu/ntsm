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
#include "vendor/tsl/robin_set.h"
#include "vendor/KseqHashIterator.hpp"
#include "vendor/kseq.h"
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class FingerPrint {
public:

	typedef uint16_t AlleleID;
	typedef uint64_t HashedKmer;

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
#pragma omp atomic update
					++m_totalKmers;
				}
#pragma omp atomic update
				m_totalBases += seq->seq.l;
				l = kseq_read(seq);
				//terminate early
				if (m_maxCounts != 0 && m_totalCounts > m_maxCounts) {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
	}

	/*
	 * Prints header for summary info
	 */
	void printHeader() const{
		string outStr = "";
		outStr += "#@TK\t";
		outStr += std::to_string(m_totalKmers);
		outStr += "\n#@KS\t";
		outStr += std::to_string(opt::k);
		outStr += "\n#locusID\tcountAT\tcountCG\tsumAT\tsumCG\tdistinctAT\tdistinctCG\n";
		cout << outStr;
	}

	void printCountsMax() const{
		printHeader();
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
				unsigned freqAlle = m_counts.at(allele1.at(j));
				if(maxCountREF < freqAlle){
					maxCountREF = freqAlle;
				}
				countSumAT += freqAlle;
			}
			for(size_t j = 0; j < allele2.size() ; ++j) {
				unsigned freqAlle = m_counts.at(allele2.at(j));
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
			cout << outStr;
		}
	}

	string printInfoSummary(){
		unsigned siteCoverage = getSitesCoveredInSample();
		string outStr = "";
		outStr += "Total Bases Considered: ";
		outStr += std::to_string(getTotalCounts());
		outStr += "\n";
		outStr += "Total k-mers Considered: ";
		outStr += std::to_string(m_totalKmers);
		outStr += "\n";
		outStr += "Total k-mers Recorded: ";
		outStr += std::to_string(getTotalKmerCounts());
		outStr += "\n";
		outStr += "Distinct k-mers in initial set: ";
		outStr += std::to_string(m_counts.size());
		outStr += "\n";
		outStr += "Total Sites: ";
		outStr += std::to_string(m_alleleIDToKmerRef.size());
		outStr += "\n";
		outStr += "Sites Covered by at least one k-mer: ";
		outStr += std::to_string(siteCoverage);
		outStr += "\n";
		if(!opt::summary.empty()){
			ofstream fh;
			fh.open(opt::summary);
			fh << outStr;
			fh.close();
		}
		double covPer = double(siteCoverage)
				/ double(m_alleleIDToKmerRef.size());
		if (covPer < opt::siteCovThreshold) {
			cerr << "Warning: site coverage is : " << covPer
					<< "(<75%). Data may be sorted or sparse along the genome. Any PCA projection may be inaccurate."
					<< endl;
		}

		return(outStr);
	}

//	void printCountsAllCounts(){
//		string tempStr;
//		for (size_t i = 0; i < m_alleleIDs.size() ; ++i) {
//			tempStr.clear();
//			const vector<uint64_t> &allele1 = *m_alleleIDToKmerRef[i];
//			const vector<uint64_t> &allele2 = *m_alleleIDToKmerVar[i];
//			tempStr += m_alleleIDs.at(i);
//			tempStr += "\tR\t";
//			for (size_t j = 0; j < allele1.size(); ++j) {
//				if(m_counts.find(allele1.at(j)) != m_counts.end()){
//					tempStr += std::to_string(m_counts[allele1.at(j)]);
//					tempStr += ",";
//				}
//			}
//			tempStr.pop_back();
//			tempStr += "\n";
//			tempStr += m_alleleIDs.at(i);
//			tempStr += "\tV\t";
//			for (size_t j = 0; j < allele2.size(); ++j) {
//				if (m_counts.find(allele2.at(j)) != m_counts.end()) {
//					tempStr += std::to_string(m_counts[allele2.at(j)]);
//					tempStr += ",";
//				}
//			}
//			tempStr.pop_back();
//			tempStr += "\n";
//			cout << tempStr;
//		}
//	}

	uint64_t getTotalKmerCounts(){
		return m_totalCounts;
	}

	uint64_t getTotalCounts(){
		return m_totalBases;
	}

	unsigned getSitesCoveredInSample(){
		unsigned count = 0;
		for (size_t i = 0; i < m_alleleIDs.size() ; ++i) {
			const vector<uint64_t> &allele1 = *m_alleleIDToKmerRef[i];
			const vector<uint64_t> &allele2 = *m_alleleIDToKmerVar[i];
			size_t maxCountREF = 0;
			size_t maxCountVAR = 0;
			for(size_t j = 0; j < allele1.size() ; ++j) {
				double freqAlle = m_counts[allele1.at(j)];
				if(maxCountREF < freqAlle){
					maxCountREF = freqAlle;
				}
			}
			for(size_t j = 0; j < allele2.size() ; ++j) {
				double freqAlle = m_counts[allele2.at(j)];
				if(maxCountVAR < freqAlle){
					maxCountVAR = freqAlle;
				}
			}
			if(maxCountREF > 0 || maxCountVAR > 0){
				++count;
			}
		}
		return count;
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

private:
	const vector<string> &m_filenames;
	uint64_t m_totalCounts;
	uint64_t m_totalKmers;
	uint64_t m_maxCounts;
	tsl::robin_map<uint64_t, size_t> m_counts; //k-mer to count
	vector<shared_ptr<vector<HashedKmer>>> m_alleleIDToKmerRef;
	vector<shared_ptr<vector<HashedKmer>>> m_alleleIDToKmerVar;
	vector<string> m_alleleIDs;
	uint64_t m_totalBases;
//	unsigned m_maxSiteSize;
//	static const unsigned interval = 65536;

	void initCountsHash(){
		gzFile fp1, fp2;
		fp1 = gzopen(opt::ref.c_str(), "r");
		fp2 = gzopen(opt::var.c_str(), "r");
		tsl::robin_set<uint64_t> dupes;
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
		{
			kseq_t *seq = kseq_init(fp1);
			int l = kseq_read(seq);
			unsigned index = 0;
			while (l >= 0) {
				assert(index == m_alleleIDToKmerRef.size());
				m_alleleIDToKmerRef.emplace_back(
						shared_ptr<vector<uint64_t>>(new vector<uint64_t>()));
//				unsigned count = 0;
				//k-merize and
				for (KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k);
						itr != itr.end(); ++itr) {
					uint64_t hv = *itr;
//					count++;
					//check for duplicates
					if (m_counts.find(hv)!= m_counts.end()) {
						cerr << "Warning: " << seq->name.s
								<< " of REF file has a k-mer collision at pos: "
								<< itr.getPos() << endl;
						dupes.insert(hv);
					} else {
						m_alleleIDToKmerRef[index]->emplace_back(
								hv);
						m_counts[hv] = 0;
					}
				}
//				if(m_maxSiteSize < count){
//					m_maxSiteSize = count;
//				}
				m_alleleIDs.emplace_back(seq->name.s);
				l = kseq_read(seq);
				index++;
			}
			kseq_destroy(seq);
			gzclose(fp1);
		}

		kseq_t *seq = kseq_init(fp2);
		int l = kseq_read(seq);
		unsigned index = 0;
		while (l >= 0) {
			assert(index == m_alleleIDToKmerVar.size());
			m_alleleIDToKmerVar.emplace_back(
					shared_ptr<vector<uint64_t>>(new vector<uint64_t>()));
			//k-merize and insert
//			unsigned count = 0;
			for (KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k);
					itr != itr.end(); ++itr) {
				uint64_t hv = *itr;
//				count++;
				//check for duplicates
				if (m_counts.find(hv) != m_counts.end()) {
					cerr << "Warning: " << seq->name.s
							<< " of VAR file has a k-mer collision at pos: "
							<< itr.getPos() << endl;
					dupes.insert(hv);
				} else {
					m_alleleIDToKmerVar[index]->emplace_back(hv);
					m_counts[hv] = 0;
				}
			}
//			if(m_maxSiteSize < count){
//				m_maxSiteSize = count;
//			}
			l = kseq_read(seq);
			index++;
		}
		kseq_destroy(seq);
		gzclose(fp2);

		if (!opt::dupes) {
			//remove dupes
			for (tsl::robin_set<uint64_t>::iterator itr = dupes.begin();
					itr != dupes.end(); ++itr) {
				m_counts.erase(*itr);
			}
		}
	}
};
#endif /* SRC_FINGERPRINT_HPP_ */
