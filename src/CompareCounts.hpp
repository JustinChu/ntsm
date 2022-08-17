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

#include "vendor/kfunc.c"

using namespace std;

class CompareCounts {
public:
	CompareCounts(const vector<string> &filenames) :
			m_filenames(filenames), m_sumlogPSingle(
					vector<double>(filenames.size())) {

		vector<uint64_t> maxCounts;
		vector<shared_ptr<vector<pair<unsigned, unsigned>>>> counts;

		//load in counts
		//open file
		//TODO parallelize loop -> may need to preallocate size of count vector
		for (vector<string>::const_iterator itr = m_filenames.begin();
				itr != m_filenames.end(); ++itr) {
			counts.push_back(
					shared_ptr<vector<pair<unsigned, unsigned>>>(
							new vector<pair<unsigned, unsigned>>));
			ifstream fh(itr->c_str());
			string line;
			if (fh.is_open()) {
				while (getline(fh, line)) {
					if (line.length() > 0) {
						std::string delimiter = "\t";
						size_t pos = line.find(delimiter);

						string locusID = line.substr(0, pos);
						line.erase(0, pos + delimiter.length());
						pos = line.find(delimiter);
						unsigned count1 = std::stoi(
								line.substr(0, pos).c_str());
						line.erase(0, pos + delimiter.length());
						pos = line.find(delimiter);
						unsigned count2 = std::stoi(
								line.substr(0, pos).c_str());
						counts.back()->push_back(
								std::make_pair(count1, count2));
						if(maxCounts.size() >= counts.back()->size()){
							if(maxCounts[counts.back()->size() - 1] < count1){
								maxCounts[counts.back()->size() - 1] = count1;
							}
							if(maxCounts[counts.back()->size() - 1] < count2){
								maxCounts[counts.back()->size() - 1] = count2;
							}
						}
						else{
							maxCounts.push_back(count1);
							assert(maxCounts.size() == counts.back()->size());
							if(maxCounts[counts.back()->size() - 1] < count2){
								maxCounts[counts.back()->size() - 1] = count2;
							}
						}
					}
				}
			}
		}

		//remove counts that exceed limits
		for (vector<shared_ptr<vector<pair<unsigned, unsigned>>>>::iterator i =
				counts.begin(); i != counts.end(); ++i) {
			m_totalCounts.push_back(0);
			m_counts.push_back(
					shared_ptr<vector<pair<unsigned, unsigned>>>(
							new vector<pair<unsigned, unsigned>>));
			for (size_t j = 0; j < maxCounts.size(); ++j) {
				unsigned count1 = (*i)->at(j).first;
				unsigned count2 = (*i)->at(j).second;

				if (!opt::maxCov || maxCounts[j] < opt::maxCov) {
					m_counts.back()->push_back(std::make_pair(count1, count2));
					m_totalCounts.back() += count1;
					m_totalCounts.back() += count2;
				}
			}
		}
	}

	void runLogLikelihood() {
		initLogPSum();
		string temp = "";
		for (unsigned i = 0; i < m_counts.size(); ++i) {
			for (unsigned j = i + 1; j < m_counts.size(); ++j) {
				double score = computeLogLikelihood(i, j)/m_counts[0]->size();
				temp += m_filenames[i];
				temp += "\t";
				temp += m_filenames[j];
				temp += "\t";
				temp += to_string(score);
				if (score < opt::scoreThresh) {
					temp += "\tY\t";
				} else {
					temp += "\tN\t";
				}
				temp += to_string(double(m_totalCounts[i])/double(m_counts[0]->size()));
				temp += "\t";
				temp += to_string(double(m_totalCounts[j])/double(m_counts[0]->size()));
				temp += "\n";
				cout << temp;
				temp.clear();
			}
		}
	}

	void runLogLikelihoodRemove() {
		vector<unsigned> validIndexes = gatherValidEntries();
		cerr << "Retained " << validIndexes.size()
				<< " SNP sites for evaluation" << endl;
		initLogPSum(validIndexes);
		string temp = "";
		for (unsigned i = 0; i < m_counts.size(); ++i) {
			for (unsigned j = i + 1; j < m_counts.size(); ++j) {
				double score = computeLogLikelihood(i, j, validIndexes)
						/ validIndexes.size();
				temp += m_filenames[i];
				temp += "\t";
				temp += m_filenames[j];
				temp += "\t";
				temp += to_string(score);
				if (score < opt::scoreThresh) {
					temp += "\tY\t";
				} else {
					temp += "\tN\t";
				}
				temp += to_string(
						double(m_totalCounts[i]) / double(m_counts[0]->size()));
				temp += "\t";
				temp += to_string(
						double(m_totalCounts[j]) / double(m_counts[0]->size()));
				temp += "\n";
				cout << temp;
				temp.clear();
			}
		}
	}

	void runLogLikelihoodRemovePairwise() {
		string temp = "";
#pragma omp parallel for
		for (unsigned i = 0; i < m_counts.size(); ++i) {
			for (unsigned j = i + 1; j < m_counts.size(); ++j) {
				unsigned indexesUsed = 0;
				double score = computeLogLikelihood(i, j, indexesUsed);
				if (opt::all || score < opt::scoreThresh) {
					temp += m_filenames[i];
					temp += "\t";
					temp += m_filenames[j];
					temp += "\t";
					temp += to_string(score);
					if (opt::all) {
						if (score < opt::scoreThresh) {
							temp += "\tY\t";
						} else {
							temp += "\tN\t";
						}
					}
					temp += to_string(
							double(m_totalCounts[i])
									/ double(m_counts[0]->size()));
					temp += "\t";
					temp += to_string(
							double(m_totalCounts[j])
									/ double(m_counts[0]->size()));
					temp += "\t";
					temp += to_string(indexesUsed);
					temp += "\n";
#pragma omp critical(cout)
					{
						cout << temp;
					}
					temp.clear();
				}
			}
		}
	}

	void runFET() {
		string temp = "";
		for (unsigned i = 0; i < m_counts.size(); ++i) {
			for (unsigned j = i + 1; j < m_counts.size(); ++j) {
				double pVal = runCombinedPval(i, j);
				temp += m_filenames[i];
				temp += "\t";
				temp += m_filenames[j];
				temp += "\t";
				temp += to_string(pVal);
				if (pVal < opt::scoreThresh) {
					temp += "\tY\t";
				} else {
					temp += "\tN\t";
				}
				temp += to_string(double(m_totalCounts[i])/double(m_counts[0]->size()));
				temp += "\t";
				temp += to_string(double(m_totalCounts[j])/double(m_counts[0]->size()));

				temp += "\n";
				cout << temp;
				temp.clear();
			}
		}
	}

	~CompareCounts() {
		// TODO Auto-generated destructor stub
	}

private:
	const vector<string> &m_filenames;
	vector<shared_ptr<vector<pair<unsigned, unsigned>>>> m_counts;
	vector<uint64_t> m_totalCounts;
	vector<double> m_sumlogPSingle;

	//new method for calculating similarity using log likelihood
	double computeSumLogPSingle(unsigned index) const{
		double sumLogP = 0;
		for (unsigned i = 0; i < m_counts[index]->size(); ++i) {
			double freqAT = 0;
			double freqCG = 0;
			if(m_counts[index]->at(i).first > opt::covThresh)
			{
				freqAT = double(m_counts[index]->at(i).first)
						/ double(
								m_counts[index]->at(i).first
										+ m_counts[index]->at(i).second);
			}
			if(m_counts[index]->at(i).second > opt::covThresh)
			{
				freqCG = double(m_counts[index]->at(i).second)
						/ double(
								m_counts[index]->at(i).first
										+ m_counts[index]->at(i).second);
			}
			sumLogP += m_counts[index]->at(i).first * freqAT
					+ m_counts[index]->at(i).second * freqCG;
		}
		return(sumLogP);
	}

	double computeSumLogPSingle(unsigned index, const vector<unsigned> &pos) const{
		double sumLogP = 0;
		for (vector<unsigned>::const_iterator i = pos.begin(); i != pos.end(); ++i) {
			double freqAT = 0;
			double freqCG = 0;
			if(m_counts[index]->at(*i).first > opt::covThresh)
			{
				freqAT = double(m_counts[index]->at(*i).first)
						/ double(
								m_counts[index]->at(*i).first
										+ m_counts[index]->at(*i).second);
			}
			if(m_counts[index]->at(*i).second > opt::covThresh)
			{
				freqCG = double(m_counts[index]->at(*i).second)
						/ double(
								m_counts[index]->at(*i).first
										+ m_counts[index]->at(*i).second);
			}
			sumLogP += m_counts[index]->at(*i).first * freqAT
					+ m_counts[index]->at(*i).second * freqCG;
		}
		return(sumLogP);
	}

	double computeSumLogPJoint(unsigned index1, unsigned index2) const{
		double sumLogP = 0;
		for (unsigned i = 0; i < m_counts[index1]->size(); ++i) {
			double freqAT = 0;
			double freqCG = 0;
			unsigned countAT = m_counts[index1]->at(i).first
					+ m_counts[index2]->at(i).first;
			unsigned countCG = m_counts[index1]->at(i).second
					+ m_counts[index2]->at(i).second;
			if (countAT > opt::covThresh) {
				freqAT = double(countAT) / double(countAT + countCG);
			}
			if (countCG > opt::covThresh) {
				freqCG = double(countCG) / double(countAT + countCG);
			}
			sumLogP += countAT * freqAT + countCG * freqCG;
		}
		return(sumLogP);
	}

	double computeSumLogPJoint(unsigned index1, unsigned index2, const vector<unsigned> &pos) const{
		double sumLogP = 0;
		for (vector<unsigned>::const_iterator i = pos.begin(); i != pos.end(); ++i) {
			double freqAT = 0;
			double freqCG = 0;
			unsigned countAT = m_counts[index1]->at(*i).first
					+ m_counts[index2]->at(*i).first;
			unsigned countCG = m_counts[index1]->at(*i).second
					+ m_counts[index2]->at(*i).second;
			if (countAT > opt::covThresh) {
				freqAT = double(countAT) / double(countAT + countCG);
			}
			if (countCG > opt::covThresh) {
				freqCG = double(countCG) / double(countAT + countCG);
			}
			sumLogP += countAT * freqAT + countCG * freqCG;
		}
		return(sumLogP);
	}

	//TODO: Easily paralizable
	vector<unsigned> gatherValidEntries() {
		vector<unsigned> valid;
		vector<bool> binValid(m_counts[0]->size(), true);
		unsigned count = 0;
		for (unsigned i = 0; i != m_counts.size(); ++i) {
			for (unsigned j = 0; j < m_counts[i]->size(); ++j) {
				if(m_counts[i]->at(j).first <= opt::covThresh && m_counts[i]->at(j).second <= opt::covThresh){
					count++;
					binValid[j] = false;
				}
			}
		}
		for(unsigned i = 0; i < binValid.size(); ++i){
			if(binValid[i]){
				valid.push_back(i);
			}
		}
		return(valid);
	}

	//TODO: Easily paralizable
	vector<unsigned> gatherValidEntries(unsigned index1, unsigned index2) {
		vector<unsigned> valid;
		vector<bool> binValid(m_counts[0]->size(), true);
		unsigned count = 0;
		for (unsigned j = 0; j < m_counts[0]->size(); ++j) {
			if(m_counts[index1]->at(j).first <= opt::covThresh && m_counts[index1]->at(j).second <= opt::covThresh){
				count++;
				binValid[j] = false;
			}
		}
		for (unsigned j = 0; j < m_counts[0]->size(); ++j) {
			if(m_counts[index2]->at(j).first <= opt::covThresh && m_counts[index2]->at(j).second <= opt::covThresh){
				count++;
				binValid[j] = false;
			}
		}
		for(unsigned i = 0; i < binValid.size(); ++i){
			if(binValid[i]){
				valid.push_back(i);
			}
		}
		return(valid);
	}

	//standard computation of Likelihood
	double computeLogLikelihood(unsigned index1, unsigned index2) const {
		return -2
				* (computeSumLogPJoint(index1, index2)
						- (m_sumlogPSingle[index1] + m_sumlogPSingle[index2]));
	}

	//compute only sites that aren't missing
	double computeLogLikelihood(unsigned index1, unsigned index2,
			unsigned &numRetained) {
		vector<unsigned> validIndexes = gatherValidEntries(index1, index2);
		numRetained = validIndexes.size();
		return -2
				* (computeSumLogPJoint(index1, index2, validIndexes)
						- (computeSumLogPSingle(index1, validIndexes)
								+ computeSumLogPSingle(index2, validIndexes)));
	}

	//computes sites that aren't missing across all datasets
	double computeLogLikelihood(unsigned index1, unsigned index2,
			const vector<unsigned> &pos) const {
		return -2
				* (computeSumLogPJoint(index1, index2, pos)
						- (m_sumlogPSingle[index1] + m_sumlogPSingle[index2]));
	}

	void initLogPSum(){
		for (unsigned i = 0; i < m_counts.size(); ++i) {
			m_sumlogPSingle[i] = computeSumLogPSingle(i);
		}
	}

	void initLogPSum(const vector<unsigned> &pos){
		for (unsigned i = 0; i < m_counts.size(); ++i) {
			m_sumlogPSingle[i] = computeSumLogPSingle(i, pos);
		}
	}


	//old method for calculating similarity using FET
	double runCombinedPval(unsigned index1, unsigned index2) {
		double sumLogPVal = 0.0;
		unsigned totalCount = 0;
		for (unsigned i = 0; i < m_counts[index1]->size(); ++i) {
			if ((m_counts[index1]->at(i).first + m_counts[index1]->at(i).second
					>= opt::covThresh)
					&& (m_counts[index2]->at(i).first
							+ m_counts[index2]->at(i).second >= opt::covThresh)) {
				double fisher_left_p, fisher_right_p, fisher_twosided_p;
				kt_fisher_exact(m_counts[index1]->at(i).first,
						m_counts[index2]->at(i).first,
						m_counts[index1]->at(i).second,
						m_counts[index2]->at(i).second, &fisher_left_p,
						&fisher_right_p, &fisher_twosided_p);
				sumLogPVal += log(fisher_twosided_p);
				totalCount++;
			}
		}
		return(exp(sumLogPVal/totalCount));
	}


};

#endif /* SRC_COMPARECOUNTS_HPP_ */

