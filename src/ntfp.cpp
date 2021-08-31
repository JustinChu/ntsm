
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

#define PROGRAM "ntfp"

void printVersion()
{
	const char VERSION_MESSAGE[] = PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	"Written by Justin Chu <cjustin@ds.dfci.harvard.edu>\n"
	"\n"
	"Copyright 2020 Dana-Farber Cancer Institute\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog(){
	const char dialog[] =
//	"Usage: " PROGRAM " build [OPTION]... [FASTA]...\n"
	"Usage: " PROGRAM " -r [FASTA] -a [FASTA] [OPTION]... [FILES...]\n"
	"  -t, --threads          Number of threads to run.[1]\n"
	"  -r, --ref              Wildtype reference fasta. [required]\n"
	"  -a, --var              Variant reference fasta. [required]\n"
	"  -k, --kmer             Kmer size use. [25]"
	"  -h, --help             Display this dialog.\n"
	"  -v, --verbose          Display verbose output.\n"
	"      --version          Print version information.\n";

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
		"threads", required_argument, NULL, 't' }, {
		"ref", required_argument, NULL, 'r' }, {
		"var", required_argument, NULL, 'a' }, {
		"kmer", required_argument, NULL, 'k' }, {
		"help", no_argument, NULL, 'h' }, {
		"version", no_argument, &OPT_VERSION, 1 }, {
		"verbose", no_argument, NULL, 'v' }, {
		NULL, 0, NULL, 0 } };

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "r:a:t:vhk:", long_options,
			&option_index)) != -1)
	{
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 'h': {
			printHelpDialog();
			break;
		}
		case 'r': {
			stringstream convert(optarg);
			if (!(convert >> opt::ref)) {
				cerr << "Error - Invalid parameter i: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'a': {
			stringstream convert(optarg);
			if (!(convert >> opt::var)) {
				cerr << "Error - Invalid parameter i: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'k': {
			stringstream convert(optarg);
			if (!(convert >> opt::k)) {
				cerr << "Error - Invalid parameter k: "
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

	FingerPrint fp(inputFiles);
	fp.computeCounts();
	fp.printCounts();
	cerr << "Memory: " << omp_get_wtime() - time << " Time:"  << Util::getRSS() << endl;
	return 0;
}




