/*
 * SeqGC.cpp
 *
 *  Created on: Dec 14, 2020
 *      Author: cjustin
 */

/*
 * PanGenomeRefBuild.cpp
 *
 *  Created on: Oct. 13, 2020
 *      Author: cjustin
 */

#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include "config.h"
#include "src/Options.h"
#include "src/Util.h"
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
	"Usage: " PROGRAM " [OPTION]...\n"
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
		"help", no_argument, NULL, 'h' }, {
		"version", no_argument, &OPT_VERSION, 1 }, {
		"verbose", no_argument, NULL, 'v' }, {
		NULL, 0, NULL, 0 } };

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "vh", long_options,
			&option_index)) != -1)
	{
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 'h': {
			printHelpDialog();
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

//#if defined(_OPENMP)
//	if (opt::threads > 0)
//	omp_set_num_threads(opt::threads);
//#endif

	if (OPT_VERSION) {
		printVersion();
	}

	return 0;
}




