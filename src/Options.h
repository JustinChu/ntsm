/*
 * Options.h
 *
 *  Created on: Oct 14, 2020
 *      Author: cjustin
 */

#ifndef OPTIONS_H
#define OPTIONS_H 1

#include <stdint.h>
#include <string>
#include <limits>

using namespace std;

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {
int verbose = 0;
unsigned threads = 1;
unsigned k = 25;
//double refRatioThrehold = 0.1;
//double minAlleleCount = 3;
string var = "";
string ref = "";
string summary = "";
float siteCovThreshold = 0.75;
double covThresh = std::numeric_limits<double>::max();

double scoreThresh = 0.5;
double covSkew = 0.2;
bool all = false;
unsigned maxCov = std::numeric_limits<unsigned>::max();
unsigned minCov = 1;
bool dupes = false;
//uint64_t minSites = 10000;
uint64_t genomeSize = 6200000000;

}
#endif
