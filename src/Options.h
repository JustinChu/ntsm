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
//Change at compile time to increase or lower number of dimensions
const unsigned dim = 20;

int verbose = 0;
unsigned threads = 1;
unsigned k = 19;
//double refRatioThrehold = 0.1;
//double minAlleleCount = 3;
string snp = "";
string summary = "";
float siteCovThreshold = 0.75;
double covThresh = std::numeric_limits<double>::max();

//PCA search criteria
double pcSearchRadius = 1.0;
double pcErrorThresh = 0.01;
double pcCovThresh = 10;
double pcLargeRadius = 20.0;
double pcMissThresh = 0.01;

string pca = "";
string norm = "";

//merge counts files filename
string merge = "";
bool onlyMerge = false;

double scoreThresh = 0.5;
double covSkew = 0.2;
bool all = false;
unsigned maxCov = std::numeric_limits<unsigned>::max();
unsigned minCov = 1;
bool dupes = false;
//uint64_t minSites = 10000;
uint64_t genomeSize = 6200000000;

string ref;
unsigned window = 31;
unsigned multi = 1;

string debug = "";
}
#endif
