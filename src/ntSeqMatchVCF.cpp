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
#include <omp.h>

#include "VCFConvert.hpp"

using namespace std;

#define PROGRAM "ntsmVCF"

void printVersion() {
	const char VERSION_MESSAGE[] =
	PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	"Written by Justin Chu <cjustin@ds.dfci.harvard.edu>\n"
	"\n"
	"Copyright 2020 Dana-Farber Cancer Institute\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog() {
	const char dialog[] =
			"Usage: " PROGRAM " -s [FASTA] -r [FASTA] [VCF]\n"
			"Converts a multi vcf file to a set of counts files.\n"
			"  -t, --threads = INT    Number of threads to run.[1]\n"
			"  -d, --dupes            Allow shared k-mers between sites to\n"
			"                         be counted.\n"
			"  -s, --snp = STR        Interleaved fasta of SNP sites to\n"
			"                         k-merize. [required]\n"
			"  -k, --kmer = INT       k-mer size used. [19]\n"
			"  -m, --multi = INT      Multiply counts by this.[1]\n"
			"  -w, --window = INT     Window size used. " + to_string(opt::window) + "\n"
			"  -r, --ref = STR        Reference fasta. [required]\n"
			"  -h, --help             Display this dialog.\n"
			"  -v, --verbose          Display verbose output.\n"
			"      --version          Print version information.\n";

	cerr << dialog << endl;
	exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[]) {
	//switch statement variable
	int c;

	//control variables
	bool die = false;
	int OPT_VERSION = 0;

	//long form arguments
	static struct option long_options[] = {
			{ "threads", required_argument, NULL, 't' },
			{ "dupes", required_argument, NULL, 'd' },
			{ "snp", required_argument, NULL, 's' },
			{ "kmer", required_argument, NULL, 'k' },
			{ "multi", required_argument, NULL, 'm' },
			{ "window", required_argument, NULL, 'w' },
			{ "ref", required_argument, NULL, 'r' },
			{ "help", no_argument, NULL, 'h' },
			{ "version", no_argument,&OPT_VERSION, 1 },
			{ "verbose", no_argument, NULL, 'v' },
			{NULL, 0, NULL, 0 } };

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "s:t:vhk:dr:w:m:", long_options,
			&option_index)) != -1) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 'h': {
			printHelpDialog();
			break;
		}

		case 'd': {
			opt::dupes = true;
			break;
		}
		case 's': {
			stringstream convert(optarg);
			if (!(convert >> opt::snp)) {
				cerr << "Error - Invalid parameter s: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'k': {
			stringstream convert(optarg);
			if (!(convert >> opt::k)) {
				cerr << "Error - Invalid parameter k: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'w': {
			stringstream convert(optarg);
			if (!(convert >> opt::window)) {
				cerr << "Error - Invalid parameter w: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'm': {
			stringstream convert(optarg);
			if (!(convert >> opt::multi)) {
				cerr << "Error - Invalid parameter m: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				cerr << "Error - Invalid parameter t: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'r': {
			stringstream convert(optarg);
			if (!(convert >> opt::ref)) {
				cerr << "Error - Invalid parameter r: " << optarg << endl;
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

	if (opt::k > 32) {
		die = true;
		cerr << "k cannot be greater than 32" << endl;
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

	if (!Util::fexists(opt::ref)){
		cerr << "Error: Unable to load reference file" << endl;
		die = true;
	}

	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	double time = omp_get_wtime();

	VCFConvert convert;
	assert(inputFiles.size() == 1);
	convert.count(inputFiles[0]);

	cerr << "Time: " << omp_get_wtime() - time << " s Memory: " << Util::getRSS()
			<< " kbytes" << endl;
	return 0;
}

