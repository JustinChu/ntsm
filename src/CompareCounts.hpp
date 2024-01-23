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
#include <algorithm>

#include "vendor/kfunc.c"
#include "vendor/tsl/robin_map.h"
#include "vendor/tsl/robin_set.h"
#include "vendor/nanoflann.hpp"
#include "KDTreeUtil.h"

using namespace std;

class CompareCounts {
public:
	CompareCounts(const vector<string> &filenames) :
			m_filenames(filenames), m_sumlogPSingle(
					vector<double>(filenames.size())), m_rawTotalCounts(
					vector<uint64_t>(filenames.size(), 0)), m_kmerSize(
					vector<unsigned>(filenames.size(), 0)), m_totalCounts(
					vector<uint64_t>(filenames.size(), 0)) {
		if(opt::verbose > 0){
			cerr << "Reading count files" << endl;
		}
		//read first file twice to init vectors
		{
			ifstream fh(m_filenames.at(0));
			string line;
			if (fh.is_open()) {
				while (getline(fh, line)) {
					if (line.length() > 0) {
						stringstream ss;
						ss.str(line);
						string item;
						getline(ss, item, '\t');
						if (line.at(0) != '#') {
							string locusID = item;
							m_locusIDToIndex[locusID] = m_locusIDs.size();
							m_locusIDs.emplace_back(locusID);
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
		m_cloud = vector<vector<double>>(filenames.size(), vector<double>(opt::dim));

#pragma omp parallel for
		for (unsigned i = 0; i < m_filenames.size(); ++i) {
			if(opt::verbose > 1){
#pragma omp critical(cerr)
				cerr << "Reading: "<< m_filenames.at(i) << endl;
			}
			ifstream fh(m_filenames.at(i));
#pragma omp critical(m_filenameToID)
			{
				m_filenameToID[m_filenames.at(i)] = i;
			}
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
							m_counts[i][m_locusIDToIndex.at(locusID)] = loadPair(ss,
									item);
							m_totalCounts[i] += m_counts[i][m_locusIDToIndex.at(
									locusID)].first
									+ m_counts[i][m_locusIDToIndex.at(locusID)].second;
							m_sum[i][m_locusIDToIndex.at(locusID)] = loadPair(ss,
									item);
						}
					}
				}
			}
		}
	}

	void projectPCs() {
		if(opt::verbose > 0){
			cerr << "Projecting samples onto PCA" << endl;
		}

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
		getline(fh, line);
		{
			stringstream ss(line);
			string val;
			//skip first entry
			ss >> val;
			//count number of components
			while (ss >> val) {
				++compNum;
			}
		}
		if(opt::verbose > 0){
			cerr << "Detected " << compNum << " components for " << normVals.size() << " sites" << endl;
		}
		assert(opt::dim <= compNum);

		vector<vector<long double>> rotVals(compNum,
				vector<long double>(normVals.size(), 0.0));
		unsigned index = 0;
		while (getline(fh, line)) {
			stringstream ss(line);
			//skip rsid
			string rsID;
			ss >> rsID;
			for (unsigned i = 0; i < compNum; ++i) {
				ss >> rotVals[i][index];
			}
			++index;
		}
		assert(index == normVals.size());
		//for each sample
#pragma omp parallel for
		for (unsigned i = 0; i < m_counts.size(); ++i) {
			//normalize
			vector<double> vals(m_counts.at(i).size(), 0.0);
			for (unsigned j = 0; j < vals.size(); ++j) {
				const pair<unsigned, unsigned> &tempCounts = m_counts.at(i).at(
						j);
				unsigned countAT = 0;
				unsigned countCG = 0;
				if(tempCounts.first > opt::minCov){
					countAT = tempCounts.first;
				}
				if(tempCounts.second > opt::minCov){
					countCG = tempCounts.second;
				}
				unsigned denom = countAT + countCG;
				//if value is missing set to zero (so it gets centered correctly)
				if(denom == 0){
					vals[j] = 0.0;
				}
				else {
					double nonNormGeno = double(countAT) / double(denom);
					vals[j] = ((nonNormGeno - 0.25) < 0.0 ? 0.0 :
								(nonNormGeno - 0.75) < 0.0 ? 0.5 : 1.0)
							- normVals[j];
//					vals[j] = nonNormGeno - normVals[j];
				}
			}
			if (opt::verbose > 2) {
#pragma omp critical(cerr)
				{
					cerr << "Normalizing " << m_filenames[i] << endl;
				}
			}
			//compute dot product for each PC
			for (unsigned j = 0; j < opt::dim; ++j) {
				m_cloud[i][j] = inner_product(vals.begin(), vals.end(),
						rotVals.at(j).begin(), 0.0);
			}
		}
		if (opt::verbose > 1) {
#pragma omp critical(cerr)
			{
				cerr << "Finished Normalization " << endl;
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

	/*
	 * Filters on 3 levels (search degree):
	 * 0 based on minimum count percent missing (first radius) < 0.01 & error rate < 0.01, radius = 2
	 * 1 based on minimum count percent missing (less stringent) < 0.3, radius = 15
	 * 2 if failing -> no radius (exhaustive search)
	 */
	void computeScorePCA() {
		if (opt::verbose > 1) {
			cerr << "Generating kd-tree" << endl;
		}
		cout << m_header;
		//build tree
		//max leaf size can be changed (10-50 seems to be fast for queries)
		kd_tree_t kdTree(opt::dim, m_cloud, { 10 });
		vector<GenotypeSummary> genotype(m_totalCounts.size());

		for (unsigned i = 0; i < m_totalCounts.size(); ++i) {
			genotype[i] = calcHomHetMiss(m_counts.at(i));
			genotype[i].errorRate = computeErrorRate(i);
			genotype[i].cov = double(m_totalCounts[i])
					/ double(m_distinct.size());
			double propMissing = double(genotype.at(i).miss)
					/ double(m_distinct.size());
			genotype[i].radius = std::numeric_limits<double>::max();
			if (genotype[i].errorRate < opt::pcErrorThresh
					&& propMissing < opt::pcMissSite1) {
				genotype[i].radius = pow(opt::pcSearchRadius1, 2);
			} else if (propMissing < opt::pcMissSite2) {
				genotype[i].radius = pow(opt::pcSearchRadius2, 2);
			}
		}
		if (opt::verbose > 1) {
			cerr << "Starting Score Computation with PCA" << endl;
		}
		if (opt::debug.empty()) {
			string temp = "\n";
			cout << temp;
#pragma omp parallel for private(temp)
			for (unsigned i = 0; i < m_counts.size(); ++i) {
				if (genotype.at(i).radius
						< std::numeric_limits<double>::max()) {
					std::vector<nanoflann::ResultItem<long unsigned, double>> ret_matches;
					size_t nMatches = kdTree.index->radiusSearch(
							&(m_cloud[i])[0], genotype.at(i).radius,
							ret_matches);
					for (size_t j = 0; j < nMatches; j++) {
						unsigned k = ret_matches[j].first;
						if (genotype.at(i).radius
								== genotype.at(k).radius) {
							if (k <= i) {
								continue;
							}
						} else if (genotype.at(i).radius
								< genotype.at(k).radius) {
							continue;
						}
						vector<unsigned> validIndexes = gatherValidEntries(i,
								k);
						double score = std::numeric_limits<double>::max();
						if (validIndexes.size() > 0) {
							score = skew(
									computeLogLikelihood(i, k, validIndexes),
									genotype[i].cov, genotype[k].cov);
							score /= double(validIndexes.size());
						}
						if (opt::all || (score < opt::scoreThresh)) {
							Relate info = calcRelatedness(i, k, validIndexes);
							resultsStr(temp, genotype, info,
									validIndexes.size(), score,
									to_string(
											calcDistance(m_cloud[i],
													m_cloud[k])), i, k);
							temp += "\n";
#pragma omp critical(cout)
							{
								cout << temp;
							}
						}
					}
				} else {
					//Search all
					for (size_t j = 0; j < m_counts.size(); j++) {
						if (std::numeric_limits<double>::max()
								== genotype.at(j).radius) {
							if (j <= i) {
								continue;
							}
						}
						vector<unsigned> validIndexes = gatherValidEntries(i,
								j);
						double score = std::numeric_limits<double>::max();
						if (validIndexes.size() > 0) {
							score = skew(
									computeLogLikelihood(i, j, validIndexes),
									genotype[i].cov, genotype[j].cov);
							score /= double(validIndexes.size());
						}
						if (opt::all || (score < opt::scoreThresh)) {
							Relate info = calcRelatedness(i, j, validIndexes);
							resultsStr(temp, genotype, info,
									validIndexes.size(), score,
									to_string(
											calcDistance(m_cloud[i],
													m_cloud[j])), i, j);
							temp += "\n";
#pragma omp critical(cout)
							{
								cout << temp;
							}
						}
					}
				}
			}
		} else {
			if(opt::verbose > 0){
				cerr << "Debug output enabled" << endl;
			}
			string temp = "\tpairs\tcandidates1\tcandidates2\tpossible\tradius1\tradius2\tcorrect\n";
			cout << temp;
			tsl::robin_set<std::pair<unsigned, unsigned>, pair_hash> truePairs;
			//load in debug file
			ifstream fh(opt::debug);
			string line;
			if (fh.is_open()) {
				while (getline(fh, line)) {
					stringstream ss(line);
					string fileID;
					vector<string> values;
					while(ss >> fileID){
						values.emplace_back(fileID);
					}
					for (unsigned i = 0; i < values.size(); ++i) {
						for (unsigned j = i + 1; j < values.size(); ++j) {
							if(!m_filenameToID.contains(values.at(i))){
								cerr << "missing file " << values.at(i) << endl;
							}
							if(!m_filenameToID.contains(values.at(j))){
								cerr << "missing file " << values.at(j) << endl;
							}
							unsigned x = m_filenameToID.at(values.at(i));
							unsigned y = m_filenameToID.at(values.at(j));
							if (x <= y) {
								truePairs.insert(std::make_pair(x,y));
							} else {
								truePairs.insert(std::make_pair(y,x));
							}
						}
					}
				}
			}
			if(opt::verbose > 0){
				cerr << "Finished creating ground truth pairs" << endl;
			}
			if(opt::all) {
				cerr << "Currently unable to output all pairs in debug mode." << endl;
				exit(1);
			}
			else {
				for (auto itr = truePairs.begin(); itr != truePairs.end();
						++itr) {
					unsigned x = itr->first;
					unsigned y = itr->second;
					vector<unsigned> validIndexes = gatherValidEntries(x, y);
					double score = std::numeric_limits<double>::max();
					if (validIndexes.size() > 0) {
						score = skew(computeLogLikelihood(x, y, validIndexes),
								genotype[x].cov, genotype[y].cov);
						score /= double(validIndexes.size());
					}
					double distance = calcDistance(m_cloud[x], m_cloud[y]);
					unsigned pairs = 0;
					for (unsigned i = 0; i < m_counts.size(); ++i) {
						std::vector<
								nanoflann::ResultItem<long unsigned int, double>> ret_matches;
						size_t nMatches = kdTree.index->radiusSearch(
								&(m_cloud[i])[0], distance, ret_matches);
						for (size_t j = 0; j < nMatches; j++) {
							auto k = ret_matches[j].first;
							if (k > i) {
								++pairs;
							}
						}
					}
					size_t candidates1 = 0;
					{
						std::vector<
								nanoflann::ResultItem<long unsigned int, double>> ret_matches;
						size_t nMatches = kdTree.index->radiusSearch(
								&(m_cloud[x])[0], genotype.at(x).radius,
								ret_matches);
						for (size_t j = 0; j < nMatches; j++) {
							unsigned k = ret_matches[j].first;
							if (genotype.at(x).radius
									== genotype.at(k).radius) {
								if (k <= x) {
									continue;
								}
							} else if (genotype.at(x).radius
									< genotype.at(k).radius) {
								continue;
							}
							++candidates1;
						}
					}
					size_t candidates2 = 0;
					{
						std::vector<
								nanoflann::ResultItem<long unsigned int, double>> ret_matches;
						size_t nMatches = kdTree.index->radiusSearch(
								&(m_cloud[y])[0], genotype.at(y).radius,
								ret_matches);
						for (size_t j = 0; j < nMatches; j++) {
							unsigned k = ret_matches[j].first;
							if (genotype.at(y).radius
									== genotype.at(k).radius) {
								if (k <= y) {
									continue;
								}
							} else if (genotype.at(y).radius
									< genotype.at(k).radius) {
								continue;
							}
							++candidates2;
						}
					}
					Relate info = calcRelatedness(x, y, validIndexes);
					resultsStr(temp, genotype, info, validIndexes.size(), score,
							to_string(calcDistance(m_cloud[x], m_cloud[y])), x,
							y);
					temp += "\t";
					temp += to_string(pairs);
					temp += "\t";
					temp += to_string(candidates1);
					temp += "\t";
					temp += to_string(candidates2);
					temp += "\t";
					temp += to_string(m_filenames.size() - 1);
					temp += "\t";
					temp += to_string(genotype.at(x).radius);
					temp += "\t";
					temp += to_string(genotype.at(y).radius);
					temp += "\t1\n";
#pragma omp critical(cout)
					{
						cout << temp;
					}
				}
			}
		}
	}

	/*
	 * Single file QC
	 * Outputs:
	 * sample: Filename for sample X
	 * cov: Coverage of sample X
	 * error_rate: Error rate of sample X. May underestimate error caused by long indels.
	 * miss: Total number of missing sites in sample X
	 * hom: Total number of homozygous sites in sample X
	 * het: Total number of heterozygous sites in sample X
	 * PCX: Principle components after projection
	 */
	void computeScoreSingle() {
		string temp =
				"sample\tcov\terrorRate\tmiss\thom\thet";
		vector<GenotypeSummary> genotype(m_totalCounts.size());
		for (unsigned i = 0; i < m_totalCounts.size(); ++i) {
			genotype[i] = calcHomHetMiss(m_counts.at(i));
			genotype[i].errorRate = computeErrorRate(i);
			genotype[i].cov = double(m_totalCounts[i])
					/ double(m_distinct.size());
		}
		//print header information for PCAs
		if(!opt::pca.empty()){
			projectPCs();
			for(size_t i = 1; i <= m_cloud[0].size(); ++i){
				temp += "\tPC";
				temp += to_string(i);
			}
		}
		cout << temp << endl;
#pragma omp parallel for private(temp)
		for (unsigned i = 0; i < m_totalCounts.size(); ++i) {
			temp.clear();
			temp += m_filenames[i];
			temp += "\t";
			temp += to_string(genotype.at(i).cov);
			temp += "\t";
			temp += to_string(genotype.at(i).errorRate);
			temp += "\t";
			temp += to_string(genotype.at(i).miss);
			temp += "\t";
			temp += to_string(genotype.at(i).homs);
			temp += "\t";
			temp += to_string(genotype.at(i).hets);
			if (!opt::pca.empty()) {
				for (size_t j = 0; j < m_cloud[i].size(); ++j) {
					temp += "\t";
					temp += to_string(m_cloud[i][j]);
				}
			}
#pragma omp critical(cout)
			{
				cout << temp;
			}
		}
	}

	/*
	 * Compute comparisons between all combinations
	 * Slower than using PCA to lower comparison numbers
	 */
	void computeScore() {
		cout << m_header;
		vector<GenotypeSummary> genotype(m_totalCounts.size());
		for (unsigned i = 0; i < m_totalCounts.size(); ++i) {
			genotype[i] = calcHomHetMiss(m_counts.at(i));
			genotype[i].errorRate = computeErrorRate(i);
			genotype[i].cov = double(m_totalCounts[i])
					/ double(m_distinct.size());
		}
		string temp = "\n";
		cout << temp;
#pragma omp parallel for private(temp)
		for (unsigned i = 0; i < m_counts.size(); ++i) {
			for (size_t j = i + 1; j < m_counts.size(); j++) {
				vector<unsigned> validIndexes = gatherValidEntries(i, j);
				double score = std::numeric_limits<double>::max();
				if (validIndexes.size() > 0) {
					score = skew(computeLogLikelihood(i, j, validIndexes),
							genotype[i].cov, genotype[j].cov);
					score /= double(validIndexes.size());
				}
				if (opt::all || (score < opt::scoreThresh)) {
					Relate info = calcRelatedness(i, j, validIndexes);
					resultsStr(temp, genotype, info, validIndexes.size(), score, "-1",
							i, j);
					temp += "\n";
#pragma omp critical(cout)
					{
						cout << temp;
					}
				}
			}
		}
	}

	void mergeCounts(){
		ofstream out(opt::merge);
		//merge tags
		uint64_t tk = 0;
		//check if k-mers are the same
		for (unsigned i = 0; i != m_kmerSize.size(); ++i) {
			for (unsigned j = i + 1; j != m_kmerSize.size(); ++j) {
				assert(m_kmerSize.at(i) == m_kmerSize.at(j));
			}
		}
		for (unsigned i = 0; i != m_rawTotalCounts.size(); ++i) {
			tk += m_rawTotalCounts.at(i);
		}
		string tmp = "#@TK\t";
		tmp += to_string(tk);
		tmp += "\n#@KS\t";
		tmp += to_string(m_kmerSize[0]);
		//TODO redundant string found in FingerPrint.hpp, refactor?
		tmp += "\n#locusID\tcountAT\tcountCG\tsumAT\tsumCG\tdistinctAT\tdistinctCG\n";
		out << tmp;
		for (unsigned i = 0; i != m_distinct.size(); ++i) {
			tmp.clear();
			unsigned countAT = 0;
			unsigned countCG = 0;
			unsigned sumAT = 0;
			unsigned sumCG = 0;
			for (unsigned j = 0; j != m_counts.size(); ++j) {
				countAT += m_counts.at(j).at(i).first;
				countCG += m_counts.at(j).at(i).second;
				sumAT += m_sum.at(j).at(i).first;
				sumCG += m_sum.at(j).at(i).second;
			}
			tmp += m_locusIDs.at(i);
			tmp += "\t";
			tmp += to_string(countAT);
			tmp += "\t";
			tmp += to_string(countCG);
			tmp += "\t";
			tmp += to_string(sumAT);
			tmp += "\t";
			tmp += to_string(sumCG);
			tmp += "\t";
			tmp += to_string(m_distinct.at(i).first);
			tmp += "\t";
			tmp += to_string(m_distinct.at(i).second);
			tmp += "\n";
			out << tmp;
		}
	}

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

	struct GenotypeSummary{
		unsigned hets = 0;
		unsigned homs = 0;
		unsigned miss = 0;
//		double madFreq = 0;
//		double median = 0;
//		double mean = 0;
//		double var50 = 0;
		double errorRate = 0;
		double cov = 0;
		double radius = 0;
	};

	const vector<string> &m_filenames;
	typedef vector<vector<pair<unsigned, unsigned>>> PairedCount;
	typedef std::vector<std::vector<double>> vector_of_vectors_t;
    typedef KDTreeVectorOfVectorsAdaptor<vector_of_vectors_t, double> kd_tree_t;

	vector<double> m_sumlogPSingle;
	PairedCount m_counts;
	PairedCount m_sum;
	vector<pair<unsigned, unsigned>> m_distinct;
	vector<uint64_t> m_rawTotalCounts;
	vector<unsigned> m_kmerSize;
	vector<uint64_t> m_totalCounts;
	vector<string> m_locusIDs;
	tsl::robin_map<string,unsigned> m_locusIDToIndex;
	vector_of_vectors_t m_cloud;
	tsl::robin_map<string,unsigned> m_filenameToID;
	const string m_header =
			"sample1\tsample2\tscore\tsame\tdist\trelate\tibs0\tibs2\thomConcord"
					"\thet1\thet2\tsharedHet\thom1\thom2\tsharedHom\tn"
					"\tcov1\tcov2\terrorRate1\terrorRate2\tmiss1\tmiss2"
					"\tallHom1\tallHom2\tallHet1\tallHet2";
//					"\tmedianGC1\tmedianGC2\tmad50GC1\tmad50GC2"
//					"\tmeanGC1\tmeanGC2\tvar50GC1\tvar50GC2";

	struct pair_hash
	{
	    template <class T1, class T2>
	    std::size_t operator() (const std::pair<T1, T2> &pair) const {
	        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	    }
	};

	GenotypeSummary calcHomHetMiss(const vector<pair<unsigned, unsigned>> &counts){
		GenotypeSummary count = { };
		vector<double> hetCount; //heterozygous frequency
		for (size_t i = 0; i < counts.size(); ++i) {
			if (counts.at(i).first > opt::minCov) {
				if (counts.at(i).second > opt::minCov) {
					++count.hets;
					hetCount.push_back(
							double(counts.at(i).first)
									/ double(
											counts.at(i).first
													+ counts.at(i).second));
				} else {
					++count.homs;
				}
			} else if (counts.at(i).second > opt::minCov) {
				++count.homs;
			} else {
				++count.miss;
			}
		}
//		count.median = 0.5;
//		count.madFreq = calculateMAD(hetCount, count.median);
//		count.mean = 0.5;
//		count.var50 = calculateVar50(hetCount, count.mean);
		return (count);
	}

	/*
	 * Calculates variance assuming equal AT vs CG counts
	 */
	double calculateVar50(const vector<double> &freq, double &mean){
		double sum = 0.0;
		for(size_t i = 0; i < freq.size(); ++i){
			sum += freq.at(i);
		}
		mean = sum/double(freq.size());
		double varSum = 0.0;

		for(size_t i = 0; i < freq.size(); ++i){
			varSum += pow(freq.at(i) - 0.5, 2);
		}
		return (varSum/freq.size());
	}

	/*
	 * Calculates the median absolute deviation
	 * Return 0.0 if there is no values within the frequency set
	 */
	double calculateMAD(const vector<double> &freq, double &median){
		if(freq.size() == 0) {
			return 0.0;
		}
				vector<double> sumCounts(freq.size(), 0);
				for(size_t i = 0; i < sumCounts.size(); ++i){
			sumCounts[i] = freq.at(i);
		}
				sort(sumCounts.begin(), sumCounts.end());
		median =
				sumCounts.size() % 2 == 0 ?
						(sumCounts[sumCounts.size() / 2]
								+ sumCounts[sumCounts.size() / 2 - 1]) / 2.0 :
						sumCounts[sumCounts.size() / 2];
				for(size_t i = 0; i < sumCounts.size(); ++i){
			sumCounts[i] = abs(sumCounts[i] - 0.5);
		}
				sort(sumCounts.begin(), sumCounts.end());
				return (sumCounts.size() % 2 == 0 ?
				(sumCounts[sumCounts.size() / 2]
						+ sumCounts[sumCounts.size() / 2 - 1]) / 2.0 :
				sumCounts[sumCounts.size() / 2]);
	}

	/*
	 * Calculates the median absolute deviation of the counts
	 * Takes the sum between each allele first before sorting
	 */
	double calculateMAD(const vector<pair<unsigned, unsigned>> &counts) const{
		vector<double> sumCounts(counts.size(), 0);
		for(size_t i = 0; i < counts.size(); ++i){
			sumCounts[i] = double(counts.at(i).first + counts.at(i).second);
		}
		sort(sumCounts.begin(), sumCounts.end());
		double median =
				sumCounts.size() % 2 == 0 ?
						(sumCounts[sumCounts.size() / 2]
								+ sumCounts[sumCounts.size() / 2 - 1]) / 2.0 :
						sumCounts[sumCounts.size() / 2];

		for(size_t i = 0; i < sumCounts.size(); ++i){
			sumCounts[i] = abs(sumCounts[i] - median);
		}
		sort(sumCounts.begin(), sumCounts.end());
		return (sumCounts.size() % 2 == 0 ?
				(sumCounts[sumCounts.size() / 2]
						+ sumCounts[sumCounts.size() / 2 - 1]) / 2.0 :
				sumCounts[sumCounts.size() / 2]);
	}

	/*
	 * prepare results string
	 */
	void resultsStr(string &temp, const vector<GenotypeSummary> &genotype,
			const Relate &info, uint64_t indexesUsed, double score, const string &dist, unsigned i,
			unsigned j) {
		temp.clear();
		temp += m_filenames[i];
		temp += "\t";
		temp += m_filenames[j];
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
		temp += dist;
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
		temp += to_string(genotype.at(i).cov);
		temp += "\t";
		temp += to_string(genotype.at(j).cov);
		temp += "\t";
		temp += to_string(genotype.at(i).errorRate);
		temp += "\t";
		temp += to_string(genotype.at(j).errorRate);
		temp += "\t";
		temp += to_string(genotype.at(i).miss);
		temp += "\t";
		temp += to_string(genotype.at(j).miss);
		temp += "\t";
		temp += to_string(genotype.at(i).homs);
		temp += "\t";
		temp += to_string(genotype.at(j).homs);
		temp += "\t";
		temp += to_string(genotype.at(i).hets);
		temp += "\t";
		temp += to_string(genotype.at(j).hets);
//		temp += "\t";
//		temp += to_string(genotype.at(i).median);
//		temp += "\t";
//		temp += to_string(genotype.at(j).median);
//		temp += "\t";
//		temp += to_string(genotype.at(i).madFreq);
//		temp += "\t";
//		temp += to_string(genotype.at(j).madFreq);
//		temp += "\t";
//		temp += to_string(genotype.at(i).mean);
//		temp += "\t";
//		temp += to_string(genotype.at(j).mean);
//		temp += "\t";
//		temp += to_string(genotype.at(i).var50);
//		temp += "\t";
//		temp += to_string(genotype.at(j).var50);
	}

	/*
	 * Computes squared Euclidian distance
	 */
	double calcDistance(const std::vector<double> &pos1, const std::vector<double> &pos2){
		double dist = 0.0;
		for (unsigned i = 0; i < opt::dim; ++i) {
			dist += pow(pos1.at(i) < pos2.at(i) ? pos2.at(i) - pos1.at(i) : pos1.at(i) - pos2.at(i), 2);
		}
		return(dist);
	}

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
			const vector<unsigned> &validIndexes) {
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
			return -1.0;
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

