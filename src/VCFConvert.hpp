/*
 * VCFConvert.hpp
 *
 *  Created on: Nov 8, 2022
 *      Author: cjustin
 *
 *  Original intent: Take in a multiVCF file, reference and outputs sequence to be counted or
 *  Rounds up coverage to value specified, even in the case of multiple counts
 *  If rsIDs exist convert to counts
 */

#ifndef SRC_VCFCONVERT_HPP_
#define SRC_VCFCONVERT_HPP_

#include <vector>
#include <string>
#include <zlib.h>
#include <math.h>

#include "FingerPrint.hpp"

#include "vendor/tsl/robin_map.h"
#include "vendor/tsl/robin_set.h"
#include "vendor/kseq.h"
#include "vendor/kseq_util.h"

#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class VCFConvert {
public:
	VCFConvert() {
		//load in reference
		gzFile fp;
		fp = gzopen(opt::ref.c_str(), "r");
		//read in seq
		kseq_t *seq = kseq_init(fp);
		int l = kseq_read(seq);
		unsigned size = 0;
		while (l >= 0) {
			m_chrIDs[seq->name.s] = m_ref.size();
			cpy_kseq(&m_ref[size++], seq);
			l = kseq_read(seq);
		}
		kseq_destroy(seq);
		gzclose(fp);
	}

	//loads in reference file
	void count(string filename) {
		//load in vcf file
		ifstream fh(filename.c_str());
		string line;
		while (getline(fh, line)) {
			if (line.at(0) == '#') {
				size_t pos = line.find("\t");
				//parse out header and IDs
				//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	...
				if (line.compare(1, pos - 1, "CHROM") == 0) {
					//skip first 9 entries
					for (unsigned i = 0; i < 9; ++i) {
						pos = line.find("\t");
						line.erase(0, pos + 1);
					}
					while (line != "") {
						pos = line.find("\t");
						m_sampleIDs.push_back(line.substr(0, pos));
						line.erase(0, pos + 1);
					}
					break;
				}
			}
		}
		vector<FingerPrint> counts(m_sampleIDs.size());
		//read remainder
		while (getline(fh, line)) {
			//chr1	911428	1:911428:C:T	C	T	.	PASS	AC=1655...	GT	0|0 ...
			size_t pos = line.find("\t");
			string chr = line.substr(0, pos);
			line.erase(0, pos + 1);
			pos = line.find("\t");
			size_t loc = stoi(line.substr(0, pos));
			line.erase(0, pos + 1);
			pos = line.find("\t");
			//skip
			line.erase(0, pos + 1);
			pos = line.find("\t");
//			string ref = line.substr(0, pos);
			line.erase(0, pos + 1);
			pos = line.find("\t");
			assert(pos == 1);
			char var = line.substr(0, pos).at(0);

			//skip 4 columns
			for (unsigned i = 0; i < 5; ++i) {
				pos = line.find("\t");
				line.erase(0, pos + 1);
			}
			//ref, var
			pair<string,string> seqs = getSeqFromSite(chr, loc, var);
			string genotype;
			unsigned index = 0;
			while (line != "") {
				pos = line.find("\t");
				genotype = line.substr(0, pos);
				if (genotype == "0|0") {
					counts[index].insertCount(seqs.first.c_str(), seqs.first.length(), opt::multi * 2);
				}
				else if (genotype == "0|1" || genotype ==  "1|0") {
					counts[index].insertCount(seqs.second.c_str(), seqs.second.length(), opt::multi);
					counts[index].insertCount(seqs.first.c_str(), seqs.first.length(), opt::multi);
				}
				else if (genotype == "1|0") {
					counts[index].insertCount(seqs.second.c_str(), seqs.second.length(), opt::multi * 2);
				}
				index++;
				line.erase(0, pos + 1);
			}
		}
		for(unsigned i = 0; i < m_sampleIDs.size(); ++i){
			ofstream out(m_sampleIDs.at(i) + ".counts.txt");
			counts[i].printCountsMax(out);
			out.close();
		}
	}

	~VCFConvert() {
		// TODO Auto-generated destructor stub
	}
private:
	vector<kseq_t> m_ref;
	vector<string> m_sampleIDs;
	tsl::robin_map<string,unsigned> m_chrIDs;

	pair<string,string> getSeqFromSite(string chr, size_t pos, char var){
		unsigned chrIndex = m_chrIDs.at(chr);
		char refStr[opt::window + 1] = {0};
		char varStr[opt::window + 1] = {0};
		size_t offset = pos - 1;
        size_t pos1 = ceil(offset - opt::window/ 2);
		strncpy(refStr, &(m_ref.at(chrIndex).seq.s[pos1]), opt::window );
		strncpy(varStr, &(m_ref.at(chrIndex).seq.s[pos1]), opt::window );
		varStr[opt::window/2] = var;
		return(std::make_pair(string(refStr), string(varStr)));
	}
};
#endif /* SRC_VCFCONVERT_HPP_ */
