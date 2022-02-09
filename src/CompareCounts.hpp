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
			m_filenames(filenames) {

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

						if(maxCounts.size() >= counts.back()->size()){
							if(maxCounts[counts.back()->size()] < count1){
								maxCounts[counts.back()->size()] = count1;
							}
							if(maxCounts[counts.back()->size()] < count2){
								maxCounts[counts.back()->size()] = count2;
							}
						}
						else{
							maxCounts.push_back(count1);
							assert(maxCounts.size() == counts.back()->size());
							if(maxCounts[counts.back()->size()] < count2){
								maxCounts[counts.back()->size()] = count2;
							}
						}

						counts.back()->push_back(
								std::make_pair(count1, count2));
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
			for (size_t j = 0; j > maxCounts.size(); ++j) {
				unsigned count1 = (*i)->at(j).first;
				unsigned count2 = (*i)->at(j).second;

				if (opt::maxCov && count1 + count2 < opt::maxCov) {
					m_counts.back()->push_back(std::make_pair(count1, count2))
					m_totalCounts.back() += count1;
					m_totalCounts.back() += count2;
				}
			}
		}
	}

	void runFET() {
		//TODO parallelize loop -> may need to manually collapse loop
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

