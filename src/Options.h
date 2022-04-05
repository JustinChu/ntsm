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

double scoreThresh = 0.01;
double covThresh = 0;
double maxCov = 0;
bool dupes = false;

}
#endif
