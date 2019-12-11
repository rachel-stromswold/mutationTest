#ifndef MAIN_H
#define MAIN_H

//#include "include/bitstrings.hpp"
#include "include/percolation.hpp"
#include <gsl/gsl_fit.h>
#include <map>
#include <iostream>
#include <string.h>
#include <fstream>
//#include <climits>
#include <chrono>

#define NUM_TRIALS	1000000
#define NUM_STEPS	32768
#define PROBABILITY	0.01
#define COL_WIDTH	80
#define INTERVAL	500

#define N_STARS		200

#define DEF_LEN		32
#define DEF_SEED	3141593

#endif //MAIN_H
