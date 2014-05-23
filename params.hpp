#ifndef PARAMS_HPP_INCLUDED
#define PARAMS_HPP_INCLUDED

#include <cmath>
#include <string>

// Default parameter values, change with command line parameters, 
//  ./a.out -seed 1 -yLen 32 -xLen 32
//
// or pipe a file through stdin:
//  config.txt | ./a.out
//
// where config.txt looks like
//  START
//  seed 1
//  xLen 32
//  yLen 32
//  STOP
//
//  START
//  seed 2
//  xLen 32
//  yLen 32
//  STOP
// etc
struct Params
{
    int seed = 12345, xLen = 4, yLen = 8, totalSteps = 300, dStep = 300;
    int maxBubbles = 2;
    double saturationNW = 0.7, linkRadMin = 0.1, linkRadMax = 0.4;
    double muWET = 0.01, muNON = 0.01, Ca = 0.005, sTens = 3.0;
    double pressDiff = 0;
    int bubbleMergeType = 2, flowType= 0, boundaryType = 0;
    std::string preText = "";
    bool SAVE = true;
    bool saveBubbles = false;
};

// Global Parameters
const double PI = 4.0*atan(1.0);
const double ERR = 1.0E-16;  // minium error in conjugate gradient solver
const double TOL = 1.0E-12;  // zero tolerance

#endif // PARAMS_HPP_INCLUDED
