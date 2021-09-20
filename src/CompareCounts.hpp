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
		//load in counts
		//open file
		//TODO parallelize loop -> may need to preallocate size of count vector
		for (vector<string>::const_iterator itr = m_filenames.begin();
				itr != m_filenames.end(); ++itr) {
			m_counts.push_back(
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
						m_counts.back()->push_back(
								std::make_pair(count1, count2));
					}
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
					temp += "\tDifferent";
				} else {
					temp += "\tSimilar";
				}

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

	double runCombinedPval(unsigned index1, unsigned index2) {
		double sumLogPVal = 0.0;
		for (unsigned i = 0; i < m_counts[index1]->size(); ++i) {
			double fisher_left_p, fisher_right_p, fisher_twosided_p;
			kt_fisher_exact(m_counts[index1]->at(i).first,
					m_counts[index2]->at(i).first,
					m_counts[index1]->at(i).second,
					m_counts[index2]->at(i).second, &fisher_left_p,
					&fisher_right_p, &fisher_twosided_p);
			sumLogPVal += log(fisher_twosided_p);
		}
		return(exp(sumLogPVal/m_counts[index1]->size()));
	}

};

#endif /* SRC_COMPARECOUNTS_HPP_ */

