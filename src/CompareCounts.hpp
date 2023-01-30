/*
 * CompareCounts.hpp
 *
 *  Created on: Aug 31, 2021
 *      Author: cjustin
 */

#ifndef SRC_COMPARECOUNTS_HPP_
#define SRC_COMPARECOUNTS_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <math.h>
#include <numeric>

#include "vendor/kfunc.c"
#include "vendor/tsl/robin_map.h"

using namespace std;

class CompareCounts {
public:
	CompareCounts(const vector<string> &filenames) :
			m_filenames(filenames), m_sumlogPSingle(
					vector<double>(filenames.size())), m_rawTotalCounts(
					vector<uint64_t>(filenames.size(), 0)), m_kmerSize(
					vector<unsigned>(filenames.size(), 0)), m_totalCounts(
					vector<uint64_t>(filenames.size(), 0)) {
		//read first file twice to init vectors
		{
			ifstream fh(m_filenames.at(0));
			string line;
			if (fh.is_open()) {
				unsigned count = 0;
				while (getline(fh, line)) {
					if (line.length() > 0) {
						stringstream ss;
						ss.str(line);
						string item;
						getline(ss, item, '\t');
						if (line.at(0) != '#') {
							string locusID = item;
							sampleIDs[locusID] = count++;
							getline(ss, item, '\t');
							getline(ss, item, '\t');

							getline(ss, item, '\t');
							getline(ss, item, '\t');

							getline(ss, item, '\t');
							m_distinct.push_back(loadPair(ss, item));
						}
					}
				}
			}
		}
		m_counts = PairedCount(filenames.size(), vector<pair<unsigned, unsigned>>(m_distinct.size()));
		m_sum = PairedCount(filenames.size(), vector<pair<unsigned, unsigned>>(m_distinct.size()));

#pragma omp parallel for
		for (unsigned i = 0; i < m_filenames.size(); ++i) {
			ifstream fh(m_filenames.at(i));
			string line;
			if (fh.is_open()) {
				while (getline(fh, line)) {
					if (line.length() > 0) {
						stringstream ss;
						ss.str(line);
						string item;
						getline(ss, item, '\t');
						if(line.at(0) == '#'){
							if (item == "#@TK") {
								getline(ss, item, '\t');
								m_rawTotalCounts[i] = std::stoull(item);
							} else if (item == "#@KS") {
								getline(ss, item, '\t');
								m_kmerSize[i] = std::stoull(item);
							}
						} else {
							string locusID = item;
							//locusID\tcountAT\tcountCG\tsumAT\\tsumCG\tdistinctAT\tdistinctCG\n
							getline(ss, item, '\t');
							m_counts[i][sampleIDs.at(locusID)] = loadPair(ss,
									item);
							m_totalCounts[i] += m_counts[i][sampleIDs.at(
									locusID)].first
									+ m_counts[i][sampleIDs.at(locusID)].second;
							m_sum[i][sampleIDs.at(locusID)] = loadPair(ss,
									item);
						}
					}
				}
			}
		}
	}

//	void runLogLikelihood() {
//		initLogPSum();
//		string temp = "";
//		for (unsigned i = 0; i < m_counts.size(); ++i) {
//			for (unsigned j = i + 1; j < m_counts.size(); ++j) {
//				double score = computeLogLikelihood(i, j)/m_counts[0]->size();
//				temp += m_filenames[i];
//				temp += "\t";
//				temp += m_filenames[j];
//				temp += "\t";
//				temp += to_string(score);
//				if (score < opt::scoreThresh) {
//					temp += "\tY\t";
//				} else {
//					temp += "\tN\t";
//				}
//				temp += to_string(double(m_totalCounts[i])/double(m_counts[0]->size()));
//				temp += "\t";
//				temp += to_string(double(m_totalCounts[j])/double(m_counts[0]->size()));
//				temp += "\n";
//				cout << temp;
//				temp.clear();
//			}
//		}
//	}

//	void runLogLikelihoodRemove() {
//		vector<unsigned> validIndexes = gatherValidEntries();
//		cerr << "Retained " << validIndexes.size()
//				<< " SNP sites for evaluation" << endl;
//		initLogPSum(validIndexes);
//		string temp = "";
//		for (unsigned i = 0; i < m_counts.size(); ++i) {
//			for (unsigned j = i + 1; j < m_counts.size(); ++j) {
//				double score = computeLogLikelihood(i, j, validIndexes)
//						/ validIndexes.size();
//				temp += m_filenames[i];
//				temp += "\t";
//				temp += m_filenames[j];
//				temp += "\t";
//				temp += to_string(score);
//				if (score < opt::scoreThresh) {
//					temp += "\tY\t";
//				} else {
//					temp += "\tN\t";
//				}
//				temp += to_string(
//						double(m_totalCounts[i]) / double(m_counts[0]->size()));
//				temp += "\t";
//				temp += to_string(
//						double(m_totalCounts[j]) / double(m_counts[0]->size()));
//				temp += "\n";
//				cout << temp;
//				temp.clear();
//			}
//		}
//	}

	void projectPCs() {
		//load in normalization values
		vector<long double> normVals;
		{
			ifstream fh(opt::norm);
			string line;
			if (fh.is_open()) {
				while (getline(fh, line)) {
					stringstream ss(line);
					long double value;
					ss >> value;
					normVals.emplace_back(value);
				}
			}
		}

		//load in rotational components
		//log number of components found
		unsigned compNum = 0;
		ifstream fh(opt::pca);
		string line;
		if (fh.is_open()) {
			//skip first line
			getline(fh, line);
			stringstream ss(line);
			string val;
			//skip first line
			ss >> val;
			ss >> val;
			//count number of components
			while(val){
				ss >> val;
				++compNum;
			}
		}
		vector<vector<long double>> rotVals(compNum, vector<long double>(normVals.size(),0));
		{
			ifstream fh(opt::pca);
			string line;
			if (fh.is_open()) {
				//skip first line
				getline(fh, line);
				unsigned index = 0;
				while (getline(fh, line)) {
					stringstream ss(line);
					//skip rsid
					string rsID;
					ss >> rsID;
					for (unsigned i = 0; 0 < compNum; ++i) {
						ss >> rotVals[i][index];
					}
					++index;
				}
			}
		}
		//for each sample
		for (unsigned i = 0; i < m_counts.size(); ++i) {
			vector<double> pcs(compNum,0);
			//normalize
			vector<double> normVals(m_counts.at(i).size(), 0.0);
			for (unsigned j = 0; j < normVals.size(); ++j) {
				const pair<unsigned, unsigned> &tempCounts = m_counts.at(i).at(
						j);
				normVals[j] = (double(tempCounts.first)
						/ double(tempCounts.first + tempCounts.second))
						- normVals[j];
			}
			//compute dot product for each PC
			for (unsigned j = 0; 0 < compNum; ++j) {
				pcs[j] = inner_product(normVals.begin(), normVals.end(), rotVals.at(j).begin(), 0);
			}
			//load into datastructure

		}
	}

	void computeScore() {
		vector<double> cov(m_totalCounts.size(),0);
		vector<double> errorRate(m_totalCounts.size(), 0);
		for (unsigned i = 0; i < m_totalCounts.size(); ++i) {
			cov[i] = double(m_totalCounts[i]) / double(m_distinct.size());
			errorRate[i] = computeErrorRate(i);
		}
		string temp = "sample1\tsample2\trelate\tibs0\tibs2\thomConcord\thets1"
				"\thets2\tsharedHets\thom1\thom2\tsharedHom\tn\tscore\tsame"
				"\tcov1\tcov2\terror_rate1\terror_rate2\n";
		cout << temp;
#pragma omp parallel for private(temp)
		for (unsigned i = 0; i < m_counts.size(); ++i) {
			for (unsigned j = i + 1; j < m_counts.size(); ++j) {
				unsigned indexesUsed = 0;
				vector<unsigned> validIndexes = gatherValidEntries(i, j);
				double score = skew(
						computeLogLikelihood(i, j, indexesUsed, validIndexes),
						cov[i], cov[j]);
				score /= double(indexesUsed);
				if (opt::all || score < opt::scoreThresh) {
					temp.clear();
					Relate info = calcRelatedness(i, j, validIndexes);
					temp += m_filenames[i];
					temp += "\t";
					temp += m_filenames[j];
					temp += "\t";
					temp += to_string(info.relatedness);
					temp += "\t";
					temp += to_string(info.ibs0);
					temp += "\t";
					temp += to_string(info.ibs2);
					temp += "\t";
					temp += to_string(info.homConcord);
					temp += "\t";
					temp += to_string(info.hets1);
					temp += "\t";
					temp += to_string(info.hets2);
					temp += "\t";
					temp += to_string(info.sharedHets);
					temp += "\t";
					temp += to_string(info.homs1);
					temp += "\t";
					temp += to_string(info.homs2);
					temp += "\t";
					temp += to_string(info.sharedHoms);
					temp += "\t";
					temp += to_string(indexesUsed);
					temp += "\t";
					temp += to_string(score);
					if (opt::all) {
						if (score < opt::scoreThresh) {
							temp += "\t1\t";
						} else {
							temp += "\t0\t";
						}
					} else {
						temp += "\t1\t";
					}
					temp += to_string(cov[i]);
					temp += "\t";
					temp += to_string(cov[j]);
					temp += "\t";
					temp += to_string(errorRate[i]);
					temp += "\t";
					temp += to_string(errorRate[j]);
					temp += "\n";
#pragma omp critical(cout)
					{
						cout << temp;
					}
				}
			}
		}
	}

//	void runFET() {
//		string temp = "";
//		for (unsigned i = 0; i < m_counts.size(); ++i) {
//			for (unsigned j = i + 1; j < m_counts.size(); ++j) {
//				double pVal = runCombinedPval(i, j);
//				temp += m_filenames[i];
//				temp += "\t";
//				temp += m_filenames[j];
//				temp += "\t";
//				temp += to_string(pVal);
//				if (pVal < opt::scoreThresh) {
//					temp += "\tY\t";
//				} else {
//					temp += "\tN\t";
//				}
//				temp += to_string(double(m_totalCounts[i])/double(m_counts[0]->size()));
//				temp += "\t";
//				temp += to_string(double(m_totalCounts[j])/double(m_counts[0]->size()));
//
//				temp += "\n";
//				cout << temp;
//				temp.clear();
//			}
//		}
//	}

	~CompareCounts() {
		// TODO Auto-generated destructor stub
	}

private:
	struct Relate{
		double relatedness = 0;
		unsigned ibs0 = 0;
		unsigned ibs2 = 0;
		double homConcord = 0;
		unsigned sharedHoms = 0;
		unsigned sharedHets = 0;
		unsigned hets1 = 0;
		unsigned homs1 = 0;
		unsigned hets2 = 0;
		unsigned homs2 = 0;
	};

	const vector<string> &m_filenames;
	typedef vector<vector<pair<unsigned, unsigned>>> PairedCount;
	vector<double> m_sumlogPSingle;
	PairedCount m_counts;
	PairedCount m_sum;
	vector<pair<unsigned, unsigned>> m_distinct;
	vector<uint64_t> m_rawTotalCounts;
	vector<unsigned> m_kmerSize;
	vector<uint64_t> m_totalCounts;
	tsl::robin_map<string,unsigned> sampleIDs;

	pair<unsigned, unsigned> loadPair(stringstream &ss, string &item) {
		unsigned count1 = std::stoul(item);
		getline(ss, item, '\t');
		unsigned count2 = std::stoul(item);
		getline(ss, item, '\t');
		return (std::make_pair(count1, count2));
	}

	//new method for calculating similarity using log likelihood
//	double computeSumLogPSingle(unsigned index) const{
//		double sumLogP = 0;
//		for (unsigned i = 0; i < m_counts.at(index).size(); ++i) {
//			double freqAT = 0;
//			double freqCG = 0;
//			if(m_counts.at(index).at(i).first > opt::minCov)
//			{
//				freqAT = double(m_counts.at(index).at(i).first)
//						/ double(
//								m_counts.at(index).at(i).first
//										+ m_counts.at(index).at(i).second);
//			}
//			if(m_counts.at(index).at(i).second > opt::minCov)
//			{
//				freqCG = double(m_counts.at(index).at(i).second)
//						/ double(
//								m_counts.at(index).at(i).first
//										+ m_counts.at(index).at(i).second);
//			}
//			sumLogP += m_counts.at(index).at(i).first * freqAT
//					+ m_counts.at(index).at(i).second * freqCG;
//		}
//		return(sumLogP);
//	}

	double computeSumLogPSingle(unsigned index, const vector<unsigned> &pos) const{
		double sumLogP = 0;
		for (vector<unsigned>::const_iterator i = pos.begin(); i != pos.end(); ++i) {
			double freqAT = 0;
			double freqCG = 0;
			if(m_counts.at(index).at(*i).first > opt::minCov)
			{
				freqAT = double(m_counts.at(index).at(*i).first)
						/ double(
								m_counts.at(index).at(*i).first
										+ m_counts.at(index).at(*i).second);
			}
			if(m_counts.at(index).at(*i).second > opt::minCov)
			{
				freqCG = double(m_counts.at(index).at(*i).second)
						/ double(
								m_counts.at(index).at(*i).first
										+ m_counts.at(index).at(*i).second);
			}
			sumLogP += m_counts.at(index).at(*i).first * freqAT
					+ m_counts.at(index).at(*i).second * freqCG;
		}
		return(sumLogP);
	}

//	double computeSumLogPJoint(unsigned index1, unsigned index2) const{
//		double sumLogP = 0;
//		for (unsigned i = 0; i < m_counts.at(index1).size(); ++i) {
//			double freqAT = 0;
//			double freqCG = 0;
//			unsigned countAT = m_counts.at(index1).at(i).first
//					+ m_counts.at(index2).at(i).first;
//			unsigned countCG = m_counts.at(index1).at(i).second
//					+ m_counts.at(index2).at(i).second;
//			if (countAT > opt::minCov) {
//				freqAT = double(countAT) / double(countAT + countCG);
//			}
//			if (countCG > opt::minCov) {
//				freqCG = double(countCG) / double(countAT + countCG);
//			}
//			sumLogP += countAT * freqAT + countCG * freqCG;
//		}
//		return(sumLogP);
//	}

	double computeSumLogPJoint(unsigned index1, unsigned index2,
			const vector<unsigned> &pos, unsigned covThresh = opt::minCov) const {
		double sumLogP = 0;
		for (vector<unsigned>::const_iterator i = pos.begin(); i != pos.end();
				++i) {
			double freqAT = 0;
			double freqCG = 0;
			unsigned countAT = m_counts.at(index1).at(*i).first
					+ m_counts.at(index2).at(*i).first;
			unsigned countCG = m_counts.at(index1).at(*i).second
					+ m_counts.at(index2).at(*i).second;
			if (countAT > covThresh) {
				freqAT = double(countAT) / double(countAT + countCG);
			}
			if (countCG > covThresh) {
				freqCG = double(countCG) / double(countAT + countCG);
			}
			sumLogP += countAT * freqAT + countCG * freqCG;
		}
		return (sumLogP);
	}

//	vector<unsigned> gatherValidEntries() {
//		vector<unsigned> valid;
//		vector<bool> binValid(m_counts[0]->size(), true);
//		unsigned count = 0;
//		for (unsigned i = 0; i != m_counts.size(); ++i) {
//			for (unsigned j = 0; j < m_counts[i]->size(); ++j) {
//				if (m_counts[i]->at(j).first <= opt::minCov
//						&& m_counts[i]->at(j).second <= opt::minCov) {
//					count++;
//					binValid[j] = false;
//				}
//			}
//		}
//		for (unsigned i = 0; i < binValid.size(); ++i) {
//			if (binValid[i]) {
//				valid.push_back(i);
//			}
//		}
//		return (valid);
//	}

	//TODO: Easily parelizable
	vector<unsigned> gatherValidEntries(unsigned index1, unsigned index2) {
		vector<unsigned> valid;
		vector<bool> binValid(m_distinct.size(), true);
		for (unsigned j = 0; j < m_distinct.size(); ++j) {
			if (m_counts.at(index1).at(j).first <= opt::minCov
					&& m_counts.at(index1).at(j).second <= opt::minCov) {
				binValid[j] = false;
			}
		}
		for (unsigned j = 0; j < m_distinct.size(); ++j) {
			if (m_counts.at(index2).at(j).first <= opt::minCov
					&& m_counts.at(index2).at(j).second <= opt::minCov) {
				binValid[j] = false;
			}
		}
		for (unsigned i = 0; i < binValid.size(); ++i) {
			if (binValid[i]) {
				valid.push_back(i);
			}
		}
		return(valid);
	}

	//skews the score by coverage
	double skew(double score, double cov1, double cov2) const{
		return (score/pow(cov1*cov2, opt::covSkew));
	}

	//standard computation of Likelihood
//	double computeLogLikelihood(unsigned index1, unsigned index2) const {
//		return -2
//				* (computeSumLogPJoint(index1, index2)
//						- (m_sumlogPSingle[index1] + m_sumlogPSingle[index2]));
//	}

	//compute only sites that aren't missing
	double computeLogLikelihood(unsigned index1, unsigned index2,
			unsigned &numRetained, const vector<unsigned> &validIndexes) {
		numRetained = validIndexes.size();
		return -2.0
				* (computeSumLogPJoint(index1, index2, validIndexes)
						- (computeSumLogPSingle(index1, validIndexes)
								+ computeSumLogPSingle(index2, validIndexes)));
	}

//	//computes sites that aren't missing across all datasets
//	double computeLogLikelihood(unsigned index1, unsigned index2,
//			const vector<unsigned> &pos) const {
//		return -2
//				* (computeSumLogPJoint(index1, index2, pos)
//						- (m_sumlogPSingle[index1] + m_sumlogPSingle[index2]));
//	}
//
//	void initLogPSum(){
//		for (unsigned i = 0; i < m_counts.size(); ++i) {
//			m_sumlogPSingle[i] = computeSumLogPSingle(i);
//		}
//	}

	void initLogPSum(const vector<unsigned> &pos){
		for (unsigned i = 0; i < m_counts.size(); ++i) {
			m_sumlogPSingle[i] = computeSumLogPSingle(i, pos);
		}
	}


	//old method for calculating similarity using FET
//	double runCombinedPval(unsigned index1, unsigned index2) {
//		double sumLogPVal = 0.0;
//		unsigned totalCount = 0;
//		for (unsigned i = 0; i < m_counts.at(index1).size(); ++i) {
//			if ((m_counts.at(index1).at(i).first + m_counts.at(index1).at(i).second
//					>= opt::minCov)
//					&& (m_counts.at(index2).at(i).first
//							+ m_counts.at(index2).at(i).second >= opt::minCov)) {
//				double fisher_left_p, fisher_right_p, fisher_twosided_p;
//				kt_fisher_exact(m_counts.at(index1).at(i).first,
//						m_counts.at(index2).at(i).first,
//						m_counts.at(index1).at(i).second,
//						m_counts.at(index2).at(i).second, &fisher_left_p,
//						&fisher_right_p, &fisher_twosided_p);
//				sumLogPVal += log(fisher_twosided_p);
//				totalCount++;
//			}
//		}
//		return(exp(sumLogPVal/totalCount));
//	}

	Relate calcRelatedness(unsigned index1, unsigned index2,
			const vector<unsigned> &validIndexes) {
		Relate info = {};
		enum AlleleType {
			HET, HOM_AT, HOM_CG, UNKNOWN
		};

		for (vector<unsigned>::const_iterator i = validIndexes.begin();
				i != validIndexes.end(); ++i) {
			AlleleType type1 = UNKNOWN, type2 = UNKNOWN;
			if (m_counts.at(index1).at(*i).first > opt::minCov) {
				if (m_counts.at(index1).at(*i).second > opt::minCov) {
					type1 = HET;
					++info.hets1;
				} else {
					type1 = HOM_AT;
					++info.homs1;
				}
			} else if (m_counts.at(index1).at(*i).second > opt::minCov) {
				type1 = HOM_CG;
				++info.homs1;
			}
			if (m_counts.at(index2).at(*i).first > opt::minCov) {
				if (m_counts.at(index2).at(*i).second > opt::minCov) {
					type2 = HET;
					++info.hets2;
				} else {
					type2 = HOM_AT;
					++info.homs2;
				}
			} else if (m_counts.at(index2).at(*i).second > opt::minCov) {
				type2 = HOM_CG;
				++info.homs2;
			}
			if (type1 == HET && type2 == HET) {
				++info.sharedHets;
				++info.ibs2;
			} else if ((type1 == HOM_AT && type2 == HOM_AT)
					|| (type1 == HOM_CG && type2 == HOM_CG)) {
				++info.sharedHoms;
				++info.ibs2;
			} else if ((type1 == HOM_CG && type2 == HOM_AT)
					|| (type1 == HOM_AT && type2 == HOM_CG)){
				++info.ibs0;
			}
		}
//		(shared-homozygous-alts(i,j)-2∗ibs0(i,j))/min (homozygous-alts(i),homozygous-alts(j))
		info.homConcord = (double(info.sharedHoms) - 2.0 * double(info.ibs0))
						/ double(info.homs1 < info.homs2 ? info.homs1 : info.homs2);
		info.relatedness = (double(info.sharedHets) - 2.0 * double(info.ibs0))
						/ double(info.hets1 < info.hets2 ? info.hets1 : info.hets2);
		return info;
	}

	double computeErrorRate(unsigned index) const{
		//1-(dat$recordedKmers*2/(dat$totalKmers*(dat$distinctKmers/genomeSize)))^(1/kmerSize)
		if (m_rawTotalCounts.at(index) > 0 && m_kmerSize.at(index) > 0) {
			uint64_t sum = 0;
			uint64_t distinctKmers = 0;
			for (unsigned i = 0; i < m_distinct.size(); ++i) {
				sum += m_sum.at(index).at(i).first
						+ m_sum.at(index).at(i).second;
				distinctKmers += m_distinct.at(i).first
						+ m_distinct.at(i).second;
			}
			double expected = double(m_rawTotalCounts.at(index))
					* double(distinctKmers) / double(opt::genomeSize);
			return (1.0
					- pow(double(sum) / expected,
							1.0 / double(m_kmerSize.at(index))));
		} else {
			return 0.0;
		}
	}

//	double calcRelatedness(unsigned index1, unsigned index2,
//			const vector<unsigned> &validIndexes ) {
//		unsigned sharedHets = 0, hets1 = 0, hets2 = 0, ibs0 = 0, sharedHom = 0, homs1 = 0, homs2 = 0;
//		enum AlleleType {HET, HOM_AT, HOM_CG, UNKNOWN};
//
//		for (vector<unsigned>::const_iterator i = validIndexes.begin(); i != validIndexes.end();
//				++i) {
//			AlleleType type1 = UNKNOWN, type2 = UNKNOWN;
//			if (m_counts.at(index1).at(*i).first > opt::minCov) {
//				if(m_counts.at(index1).at(*i).second > opt::minCov){
//					type1 = HET;
//					++hets1;
//				}
//				else{
//					type1 = HOM_AT;
//
//				}
//			}
//			else if(m_counts.at(index1).at(*i).second > opt::minCov){
//				type1 = HOM_CG;
//			}
//			if (m_counts.at(index2).at(*i).first > opt::minCov) {
//				if(m_counts.at(index2).at(*i).second > opt::minCov){
//					type2 = HET;
//					++hets2;
//				}
//				else{
//					type2 = HOM_AT;
//				}
//			}
//			else if(m_counts.at(index2).at(*i).second > opt::minCov){
//				type2 = HOM_CG;
//			}
//			if(type1 == HET && type2 == HET){
//				++sharedHets;
//			}
//			else if ((type1 == HOM_AT && type2 == HOM_CG) || (type1 == HOM_CG && type2 == HOM_AT)){
//				++ibs0;
//			}
//		}
////		cout << sharedHets << " " << ibs0 << " "  << (hets1 < hets2 ? hets1 : hets2) << endl;
//		return (double(sharedHets)-2.0*double(ibs0))/double(hets1 < hets2 ? hets1 : hets2);
//	}


//	double calcRelatedness(unsigned index1, unsigned index2,
//			const vector<unsigned> &validIndexes ) {
//		uint64_t sharedHets = 0, hets1 = 0, hets2 = 0, ibs0 = 0;
//		enum AlleleType {HET, HOM_AT, HOM_CG, UNKNOWN};
//
//		for (vector<unsigned>::const_iterator i = validIndexes.begin(); i != validIndexes.end();
//				++i) {
//			AlleleType type1 = UNKNOWN, type2 = UNKNOWN;
//			if (m_counts.at(index1).at(*i).first > opt::minCov) {
//				if(m_counts.at(index1).at(*i).second > opt::minCov){
//					type1 = HET;
//					hets1 += m_counts.at(index1).at(*i).first + m_counts.at(index1).at(*i).second;
//				}
//				else{
//					type1 = HOM_AT;
//				}
//			}
//			else if(m_counts.at(index1).at(*i).second > opt::minCov){
//				type1 = HOM_CG;
//			}
//			if (m_counts.at(index2).at(*i).first > opt::minCov) {
//				if(m_counts.at(index2).at(*i).second > opt::minCov){
//					type2 = HET;
//					hets2 += m_counts.at(index2).at(*i).first + m_counts.at(index2).at(*i).second;
//				}
//				else{
//					type2 = HOM_AT;
//				}
//			}
//			else if(m_counts.at(index2).at(*i).second > opt::minCov){
//				type2 = HOM_CG;
//			}
//			if (type1 == HET && type2 == HET) {
//				sharedHets += m_counts.at(index1).at(*i).first
//						+ m_counts.at(index1).at(*i).second
//						+ m_counts.at(index2).at(*i).first
//						+ m_counts.at(index2).at(*i).second;
//			} else if ((type1 == HOM_AT && type2 == HOM_CG)
//					|| (type1 == HOM_CG && type2 == HOM_AT)) {
//				if(m_counts.at(index1).at(*i).first > m_counts.at(index1).at(*i).second){
//					ibs0 += m_counts.at(index1).at(*i).first + m_counts.at(index2).at(*i).second;
//				}
//				else{
//					ibs0 += m_counts.at(index1).at(*i).second + m_counts.at(index2).at(*i).first;
//				}
//			}
//		}
//		cout << sharedHets << " " << ibs0 << " "  << (hets1 + hets2) << endl;
//		return (double(sharedHets)-2.0*double(ibs0))/double(hets1 + hets2);
//	}
};

#endif /* SRC_COMPARECOUNTS_HPP_ */

