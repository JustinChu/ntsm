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
#include "vendor/concurrentqueue.h"
#include "vendor/kseq_util.h"
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "vendor/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class FingerPrint {
public:

	typedef uint16_t AlleleID;
	typedef uint64_t HashedKmer;

	FingerPrint() : m_totalCounts(0), m_maxCounts(0), m_totalBases(0), m_totalReads(0) {
		//read in fasta files
		//generate hash table
		initCountsHash();
		if(opt::covThresh != 0){
			m_maxCounts = (m_counts.size() * opt::covThresh)/2;
		}
	}

	void computeCounts(const vector<string> &filenames){
#pragma omp parallel for
		for (unsigned i = 0; i < filenames.size(); ++i) {
			gzFile fp;
			fp = gzopen(filenames[i].c_str(), "r");
			if (fp == Z_NULL) {
				std::cerr << "file " << filenames[i] << " cannot be opened"
						<< std::endl;
				exit(1);
			} else if (opt::verbose) {
				std::cerr << "Opening " << filenames[i] << std::endl;
			}
			//read in seq
			kseq_t *seq = kseq_init(fp);
			int l = kseq_read(seq);
			while (l >= 0 && !m_earlyTerm) {
				processSingleRead(seq);
				l = kseq_read(seq);
				if (opt::verbose > 2 && (m_totalReads % 100000) == 0) {
					cerr << "Current Total: " << m_totalKmers << " reads, "
							<< m_totalCounts << " total counts, and "
							<< m_totalBases << " total bases " << endl;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
	}

	void insertCount(const char *seqs, uint64_t seql, unsigned multiplier = 1) {
		for (KseqHashIterator itr(seqs, seql, opt::k); itr != itr.end();
				++itr) {
			if (m_counts.find(*itr) != m_counts.end()) {
#pragma omp atomic update
				m_counts[*itr] += multiplier;
#pragma omp atomic update
				m_totalCounts += multiplier;
			}
#pragma omp atomic update
			++m_totalKmers;
		}
#pragma omp atomic update
		m_totalBases += seql;
	}

	//use only if threads > number of files
	void computeCountsProducerConsumer(const vector<string> &filenames) {
		if (opt::threads <= filenames.size()) {
			//not enough threads to saturate
			computeCounts(filenames);
		} else {
			uint64_t numReads = 0, processedCount = 0;

			moodycamel::ConcurrentQueue<kseq_t> workQueue(
					opt::threads * s_bulkSize);
			moodycamel::ConcurrentQueue<kseq_t> recycleQueue(
					opt::threads * s_bulkSize * 2);
			bool good = true;
			typedef std::vector<kseq_t>::iterator iter_t;

			//fill recycleQueue with empty objects
			{
				std::vector<kseq_t> buffer(opt::threads * s_bulkSize * 2,
						kseq_t());
				recycleQueue.enqueue_bulk(
						std::move_iterator<iter_t>(buffer.begin()),
						buffer.size());
			}

#pragma omp parallel
			{
				std::vector<kseq_t> readBuffer(s_bulkSize);
				string outBuffer;
				if (unsigned(omp_get_thread_num()) < filenames.size()) {
					//file reading init
					gzFile fp;
					fp = gzopen(filenames.at(omp_get_thread_num()).c_str(), "r");
					kseq_t *seq = kseq_init(fp);

					//per thread token
					moodycamel::ProducerToken ptok(workQueue);

					//tokens for recycle queue
					moodycamel::ConsumerToken rctok(recycleQueue);
					moodycamel::ProducerToken rptok(recycleQueue);

					unsigned dequeueSize = recycleQueue.try_dequeue_bulk(rctok,
							std::move_iterator<iter_t>(readBuffer.begin()),
							s_bulkSize);
					while (dequeueSize == 0) {
						dequeueSize = recycleQueue.try_dequeue_bulk(rctok,
								std::move_iterator<iter_t>(readBuffer.begin()),
								s_bulkSize);
					}

					unsigned size = 0;
					while (kseq_read(seq) >= 0 && !m_earlyTerm) {
						cpy_kseq(&readBuffer[size++], seq);
						if (dequeueSize == size) {
							//try to insert, if cannot queue is full
							while (!workQueue.try_enqueue_bulk(ptok,
									std::move_iterator<iter_t>(
											readBuffer.begin()), size)) {
								//try to work
								if (kseq_read(seq) >= 0) {
									//------------------------WORK CODE START---------------------------------------
									processSingleRead(seq);
									//------------------------WORK CODE END-----------------------------------------
								} else {
									goto fileEmpty;
								}
							}
							//reset buffer
							dequeueSize = recycleQueue.try_dequeue_bulk(rctok,
									std::move_iterator<iter_t>(
											readBuffer.begin()), s_bulkSize);
							while (dequeueSize == 0) {
								//try to work
								if (kseq_read(seq) >= 0) {
									//------------------------WORK CODE START---------------------------------------
									processSingleRead(seq);
									//------------------------WORK CODE END-----------------------------------------
								} else {
									goto fileEmpty;
								}
								dequeueSize = recycleQueue.try_dequeue_bulk(
										rctok,
										std::move_iterator<iter_t>(
												readBuffer.begin()),
										s_bulkSize);
							}
							size = 0;
						}
					}
					fileEmpty:
					//finish off remaining work
					for (unsigned i = 0; i < size; ++i) {
						//------------------------WORK CODE START---------------------------------------
						processSingleRead(seq);
						//------------------------WORK CODE END-----------------------------------------
					}
					assert(
							recycleQueue.enqueue_bulk(rptok,
									std::move_iterator<iter_t>(
											readBuffer.begin()), size));
					if (processedCount < numReads) {
						moodycamel::ConsumerToken ctok(workQueue);
						//join in if others are still not finished
						while (processedCount < numReads) {
							size_t num = workQueue.try_dequeue_bulk(ctok,
									std::move_iterator<iter_t>(
											readBuffer.begin()), s_bulkSize);
							if (num) {
								for (unsigned i = 0; i < num; ++i) {
									//------------------------WORK CODE START---------------------------------------
									processSingleRead(seq);
									//------------------------WORK CODE END-----------------------------------------
								}
								assert(
										recycleQueue.enqueue_bulk(rptok,
												std::move_iterator<iter_t>(
														readBuffer.begin()),
												num));
							}
						}
					}
#pragma omp atomic update
					good &= false;
					kseq_destroy(seq);
					gzclose(fp);
				} else {
					moodycamel::ConsumerToken ctok(workQueue);
					moodycamel::ProducerToken rptok(recycleQueue);
					while (good) {
						if (workQueue.size_approx() >= s_bulkSize) {
							size_t num = workQueue.try_dequeue_bulk(ctok,
									std::move_iterator<iter_t>(
											readBuffer.begin()), s_bulkSize);
							if (num) {
								for (unsigned i = 0; i < num; ++i) {
									//------------------------WORK CODE START---------------------------------------
									processSingleRead(&readBuffer[i]);
									//------------------------WORK CODE END-----------------------------------------
								}
								assert(
										recycleQueue.enqueue_bulk(rptok,
												std::move_iterator<iter_t>(
														readBuffer.begin()),
												num));
							}
						}
					}
				}
			}
		}
	}


	/*
	 * Prints summary info needed for computing error rate
	 */
	void printOptionalHeader() const{
		string outStr = "";
		outStr += "#@TK\t";
		outStr += std::to_string(m_totalKmers);
		outStr += "\n#@KS\t";
		outStr += std::to_string(opt::k);
		cout << outStr;
	}

	void printCountsMax(ostream &out = cout) const{
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
			out << outStr;
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
	const static size_t s_bulkSize = 1024;
	uint64_t m_totalCounts;
	uint64_t m_totalKmers;
	uint64_t m_maxCounts;
	tsl::robin_map<uint64_t, size_t> m_counts; //k-mer to count
	vector<shared_ptr<vector<HashedKmer>>> m_alleleIDToKmerRef;
	vector<shared_ptr<vector<HashedKmer>>> m_alleleIDToKmerVar;
	vector<string> m_alleleIDs;
	uint64_t m_totalBases;
	uint64_t m_totalReads; //for debugging purposes
	bool m_earlyTerm;
//	unsigned m_maxSiteSize;
//	static const unsigned interval = 65536;

	void processSingleRead(kseq_t *seq){
#pragma omp atomic
		m_totalReads++;
		//k-merize and insert
		insertCount(seq->seq.s, seq->seq.l);
		if (m_maxCounts != 0 && m_totalCounts > m_maxCounts) {
			if (opt::verbose > 2) {
				cerr << "max count reached at " << m_totalKmers << " reads, "
						<< m_totalCounts << " total counts, and "
						<< m_totalBases << " total bases " << endl;
			}
			m_earlyTerm = true;
		}
	}

	void initCountsHash(){
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
			size_t entryNum = 0;
			while (l >= 0) {
				if (entryNum % 2 == 0) {
					unsigned index = entryNum / 2;
					assert(index == m_alleleIDToKmerRef.size());
					m_alleleIDToKmerRef.emplace_back(
							shared_ptr<vector<uint64_t>>(
									new vector<uint64_t>()));
					//k-merize and
					for (KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k);
							itr != itr.end(); ++itr) {
						uint64_t hv = *itr;
						//check for duplicates
						if (m_counts.find(hv) != m_counts.end()) {
							cerr << "Warning: " << seq->name.s
									<< " of REF file has a k-mer collision at pos: "
									<< itr.getPos() << endl;
							dupes.insert(hv);
						} else {
							m_alleleIDToKmerRef[index]->emplace_back(hv);
							m_counts[hv] = 0;
						}
					}
					m_alleleIDs.emplace_back(seq->name.s);
				} else {
					unsigned index = entryNum / 2;
					assert(index == m_alleleIDToKmerVar.size());
					m_alleleIDToKmerVar.emplace_back(
							shared_ptr<vector<uint64_t>>(
									new vector<uint64_t>()));
					//k-merize and insert
					for (KseqHashIterator itr(seq->seq.s, seq->seq.l, opt::k);
							itr != itr.end(); ++itr) {
						uint64_t hv = *itr;
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
				}
				l = kseq_read(seq);
				entryNum++;
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
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
