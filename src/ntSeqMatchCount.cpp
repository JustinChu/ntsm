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

#include "FingerPrint.hpp"

#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "ntsmCount"

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
			"Usage: " PROGRAM " -s [FASTA] [OPTION]... [FILES...]\n"
			"  -t, --threads = INT    Number of threads to run.[1]\n"
			"  -m, --maxCov = INT     k-mer coverage threshold for early\n"
			"                         termination. [inf]\n"
			"  -o, --output = STR     Output for summary file.\n"
			"  -d, --dupes            Allow shared k-mers between sites to\n"
			"                         be counted.\n"
			"  -s, --snp = STR        Interleaved fasta of SNP sites to\n"
			"                         k-merize. [required]\n"
			"  -k, --kmer = INT       k-mer size used. [19]\n"
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
			{ "maxCov", required_argument, NULL, 'm' },
			{ "output", required_argument, NULL, 'o' },
			{ "dupes", required_argument, NULL, 'd' },
			{ "snp", required_argument, NULL, 's' },
			{ "kmer", required_argument, NULL, 'k' },
			{ "help", no_argument, NULL, 'h' },
			{ "version", no_argument,&OPT_VERSION, 1 },
			{ "verbose", no_argument, NULL, 'v' },
			{NULL, 0, NULL, 0 } };

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "s:t:vhk:m:do:", long_options,
			&option_index)) != -1) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 'h': {
			printHelpDialog();
			break;
		}
		case 'o': {
			stringstream convert(optarg);
			if (!(convert >> opt::summary)) {
				cerr << "Error - Invalid parameter o: " << optarg << endl;
				return 0;
			}
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
		case 'm': {
			stringstream convert(optarg);
			if (!(convert >> opt::covThresh)) {
				cerr << "Error - Invalid parameter m: " << optarg << endl;
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
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				cerr << "Error - Invalid parameter t: " << optarg << endl;
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
		cerr << "Error: k cannot be greater than 32" << endl;
	}

	if (opt::snp.empty()) {
		die = true;
		cerr << "Error: Missing variants (-s) file" << endl;
	}

	vector<string> inputFiles;
	while (optind < argc) {
		inputFiles.emplace_back(argv[optind]);
		assert(Util::fexists(inputFiles.back()));
		optind++;
	}

	//Check needed options
	if (inputFiles.size() == 0) {
		cerr << "Error: Need input files" << endl;
		die = true;
	}

	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	double time = omp_get_wtime();

	FingerPrint fp;
	fp.computeCounts(inputFiles);
	fp.printOptionalHeader();
	fp.printCountsMax();
	cerr << fp.printInfoSummary() << endl;
	cerr << "Time: " << omp_get_wtime() - time << " s Memory: " << Util::getRSS()
			<< " kbytes" << endl;
	return 0;
}

