/*
 * Compares count lists using Fisher's Exact Test per pair and produces a table
 * with the combinations of all pairs
 */

#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include "config.h"
#include "src/Options.h"
#include "src/Util.h"
#include "CompareCounts.hpp"

#include <omp.h>

using namespace std;

#define PROGRAM "ntsmEval"

void printVersion()
{
	const char VERSION_MESSAGE[] = PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	"Written by Justin Chu <cjustin@ds.dfci.harvard.edu>\n"
	"\n"
	"Copyright 2021 Dana-Farber Cancer Institute\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog(){
	const string dialog =
	"Usage: " PROGRAM " [FILES...]\n"
    "Processes sets of counts files and compares their similarity.\n"
    "If only a single file is provided general QC information returned.\n"
	"  -t, --threads              Number of threads to run.[1]\n"
	"  -s, --score_thresh = FLOAT Score threshold ["+ to_string(opt::scoreThresh)+"]\n"
	"  -a, --all                  Output results of all tests tried, not just those that\n"
	"                             pass the score threshold.\n"
	"  -w, --skew = FLOAT         Divides the score by coverage. Formula: (cov1*cov2)^skew\n"
	"                             Set to zero for no skew.[" + to_string(opt::covSkew) + "]\n"
	"  -c, --min_cov = INT        Keep only sites with this coverage and above.[" + to_string(opt::minCov) + "]\n"
	"  -g, --genome_size = INT    Diploid genome size for error rate estimation.\n"
	"                             ["+ to_string(opt::genomeSize)+"]\n"
	"  -e, --merge = STR          After analysis merge counts and output to file.\n"
	"  -o, --only_merge           Do not perform an analysis. Only functions when\n"
	"                             -e (--merge) option is specified.\n"
	"  -p, --pca = STR            Use PCA information to speed up analysis. Input is a\n"
	"                             set of rotational values from a PCA.\n"
	"  -d, --dim = INT            Number of dimensions to consider in PCA. [" + to_string(opt::dim) + "]\n"
	"  -n, --norm = STR           Set of values use to center the data before rotation\n"
	"                             during PCA. [Required if -p is enabled]\n"
	"  -r, --error_rate = FLOAT   Error rate threshold for PCA based search [" + to_string(opt::pcErrorThresh) + "]\n"
	"  -1, --miss_small = FLOAT   Missing site threshold small for PCA based search [" + to_string(opt::pcMissSite1) + "]\n"
	"  -2, --miss_large = FLOAT   Missing site threshold large PCA based search [" + to_string(opt::pcMissSite2) + "]\n"
	"  -S, --small = FLOAT        Search radius for small PCA based search [" + to_string(opt::pcSearchRadius1) + "]\n"
	"  -l, --large = FLOAT        Search radius for large PCA based search [" + to_string(opt::pcSearchRadius2) + "]\n"
	"  -h, --help                 Display this dialog.\n"
	"  -v, --verbose              Display verbose output.\n"
	"      --version              Print version information.\n";

	cerr << dialog << endl;
	exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
	//switch statement variable
	int c;

	//control variables
	bool die = false;
	int OPT_VERSION = 0;

	//long form arguments
	static struct option long_options[] = { {
		"score_thresh", required_argument, NULL, 's' }, {
		"all", no_argument, NULL, 'a' }, {
		"min_cov", required_argument, NULL, 'c' }, {
		"max_cov", required_argument, NULL, 'm' }, {
		"skew", required_argument, NULL, 'w' }, {
		"genome_size", required_argument, NULL, 'g' }, {
		"threads", required_argument, NULL, 't' }, {
		"merge", required_argument, NULL, 'e' }, {
		"only_merge", required_argument, NULL, 'o' }, {
		"help", no_argument, NULL, 'h' }, {
	    "pca", required_argument, NULL, 'p' }, {
		"norm", required_argument, NULL, 'n' }, {\
		"error_rate", required_argument, NULL, 'r' }, {
		"miss_small", required_argument, NULL, '1' }, {
		"miss_large", required_argument, NULL, '2' }, {
		"small", required_argument, NULL, 'k' }, {
		"large", required_argument, NULL, 'l' }, {
		"debug", required_argument, NULL, 'b' }, {
		"version", no_argument, &OPT_VERSION, 1 }, {
		"verbose", no_argument, NULL, 'v' }, {
		NULL, 0, NULL, 0 } };

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "t:vhs:c:m:aw:g:p:n:d:r:e:o1:2:S:l:b:", long_options,
			&option_index)) != -1)
	{
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 'h': {
			printHelpDialog();
			break;
		}
		case 'a': {
			opt::all = true;
			break;
		}
		case 's': {
			stringstream convert(optarg);
			if (!(convert >> opt::scoreThresh)) {
				cerr << "Error - Invalid parameter s: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'w': {
			stringstream convert(optarg);
			if (!(convert >> opt::covSkew)) {
				cerr << "Error - Invalid parameter w: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'c': {
			stringstream convert(optarg);
			if (!(convert >> opt::minCov)) {
				cerr << "Error - Invalid parameter c: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'm': {
			stringstream convert(optarg);
			if (!(convert >> opt::maxCov)) {
				cerr << "Error - Invalid parameter m: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'g': {
			stringstream convert(optarg);
			if (!(convert >> opt::genomeSize)) {
				cerr << "Error - Invalid parameter g: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				cerr << "Error - Invalid parameter t: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'e': {
			stringstream convert(optarg);
			if (!(convert >> opt::merge)) {
				cerr << "Error - Invalid parameter e: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'o': {
			opt::onlyMerge = true;
			break;
		}
		case 'p': {
			stringstream convert(optarg);
			if (!(convert >> opt::pca)) {
				cerr << "Error - Invalid parameter p: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'n': {
			stringstream convert(optarg);
			if (!(convert >> opt::norm)) {
				cerr << "Error - Invalid parameter n: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'r': {
			stringstream convert(optarg);
			if (!(convert >> opt::pcErrorThresh)) {
				cerr << "Error - Invalid parameter r: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case '1': {
			stringstream convert(optarg);
			if (!(convert >> opt::pcMissSite1)) {
				cerr << "Error - Invalid parameter 1: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case '2': {
			stringstream convert(optarg);
			if (!(convert >> opt::pcMissSite2)) {
				cerr << "Error - Invalid parameter 2: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'S': {
			stringstream convert(optarg);
			if (!(convert >> opt::pcSearchRadius1)) {
				cerr << "Error - Invalid parameter S: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'l': {
			stringstream convert(optarg);
			if (!(convert >> opt::pcSearchRadius2)) {
				cerr << "Error - Invalid parameter l: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'd': {
			stringstream convert(optarg);
			if (!(convert >> opt::dim)) {
				cerr << "Error - Invalid parameter d: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'b': {
			stringstream convert(optarg);
			if (!(convert >> opt::debug)) {
				cerr << "Error - Invalid parameter b: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'v': {
			opt::verbose++;
			break;
		}
		case '?': {
			die = true;
			break;
		}
		}
	}

#if defined(_OPENMP)
	if (opt::threads > 0)
	omp_set_num_threads(opt::threads);
#endif

	if (OPT_VERSION) {
		printVersion();
	}

	vector<string> inputFiles;
	while (optind < argc) {
		inputFiles.emplace_back(argv[optind]);
		assert(Util::fexists(inputFiles.back()));
		optind++;
	}

	//Check needed options
	if (inputFiles.size() == 0) {
		cerr << "Error: Need Input File" << endl;
		die = true;
	}
	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	double time = omp_get_wtime();

	CompareCounts comp(inputFiles);
	if(inputFiles.size() == 1){
		if(opt::verbose > 1){
			cerr << "Detected only 1 file, providing only QC information." << endl;
		}
		comp.computeScoreSingle();
	}
	else{
		if(opt::verbose > 1){
			cerr << "Finished loading files. Now comparing all samples." << endl;
		}
		if(opt::onlyMerge){
			if(opt::merge.empty()){
				cerr << "(-l) cannot be used without --merge (-e) option." << endl;
				exit(EXIT_FAILURE);
			}
			else{
				cerr << " (-l) option detected. Not performing analysis, only merging." << endl;
			}
		}
		else{
			if(opt::pca.empty()){
				cerr << "Performing all-to-all score computation.\nSpecify -p (--pca) to enable faster comparisons." << endl;
				comp.computeScore();
			}
			else{
				if (!Util::fexists(opt::norm)) {
					cerr << "Error: Need normalization file" << endl;
					die = true;
				}
				comp.projectPCs();
				comp.computeScorePCA();
			}
		}
		if(!opt::merge.empty()){
			comp.mergeCounts();
		}
	}

	cerr << "Time: " << omp_get_wtime() - time << " s Memory: "  << Util::getRSS() << " kbytes" << endl;
	return 0;
}




