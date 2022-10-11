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
			"Usage: " PROGRAM " -r [FASTA] -a [FASTA] [OPTION]... [FILES...]\n"
			"  -t, --threads = INT    Number of threads to run.[1]\n"
			"  -m, --maxCov = INT     k-mer coverage threshold for early\n"
			"                         termination [inf].\n"
//			"  -c, --con_thread       Number of threads in consumer threading.\n"
//			"                         In this mode the number of threads used\n"
//			"                         will be equal to the number of files\n"
//			"                         plus the number of producer threads.[0]\n"
			"  -o, --output = STR     Output for summary file.\n"
			"  -d, --dupes            Allow shared k-mers between sites to be\n"
			"                         counted.\n"
			"  -r, --ref = STR        Wildtype reference fasta. [required]\n"
			"  -a, --var = STR        Variant reference fasta. [required]\n"
			"  -k, --kmer = INT       Kmer size use. [25]\n"
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
	static struct option long_options[] = { { "threads", required_argument, NULL, 't' },
			{ "maxCov", required_argument, NULL, 'm' },
			{ "output", required_argument, NULL, 'o' },
			{ "dupes", required_argument, NULL, 'd' },
			{ "ref", required_argument, NULL, 'r' },
			{ "var",required_argument, NULL, 'a' },
			{ "kmer", required_argument, NULL, 'k' },
			{ "help", no_argument, NULL, 'h' },
			{ "version", no_argument,&OPT_VERSION, 1 },
			{ "verbose", no_argument, NULL, 'v' },
			{NULL, 0, NULL, 0 } };

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "r:a:t:vhk:m:do:", long_options,
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
		case 'r': {
			stringstream convert(optarg);
			if (!(convert >> opt::ref)) {
				cerr << "Error - Invalid parameter r: " << optarg << endl;
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
		case 'a': {
			stringstream convert(optarg);
			if (!(convert >> opt::var)) {
				cerr << "Error - Invalid parameter a: " << optarg << endl;
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
		cerr << "k cannot be greater than 32" << endl;
		printVersion();
	}

	if (opt::k > 32) {
		die = true;
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

	FingerPrint fp(inputFiles);
	fp.computeCountsProducerConsumer();
	fp.printCountsMax();
	cerr << fp.printInfoSummary() << endl;
	cerr << "Time: " << omp_get_wtime() - time << " s Memory: " << Util::getRSS()
			<< " kbytes" << endl;
	return 0;
}

