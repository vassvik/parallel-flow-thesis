#include "functions.hpp"
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>

// COMPILE USING c++ *.cpp -std=c++11 -O3


// funcion to parse the command line arguments
bool parseCommandLine(int argc, char * argv[], Params & params)
{
    if (argc == 0)
        return false;
    else
    {
        for (int i = 1; i < argc; i+=2)
        {
            std::cout << argv[i] << " " << argv[i+1] << std::endl;
            if (strcmp( argv[i], "-xLen") == 0)                  params.xLen = atoi(argv[i+1]);
            else if (strcmp( argv[i], "-yLen") == 0)             params.yLen = atoi(argv[i+1]);
            else if (strcmp( argv[i], "-totalSteps") == 0)       params.totalSteps = atoi(argv[i+1]);
            else if (strcmp( argv[i], "-dStep") == 0)            params.dStep = atoi(argv[i+1]);
            else if (strcmp( argv[i], "-maxBubbles") == 0)       params.maxBubbles = atoi(argv[i+1]);
            else if (strcmp( argv[i], "-seed") == 0)             params.seed = atoi(argv[i+1]);
            else if (strcmp( argv[i], "-bubbleMergeType") == 0)  params.bubbleMergeType = atoi(argv[i+1]);
            else if (strcmp( argv[i], "-flowType") == 0)         params.flowType = atoi(argv[i+1]);
            else if (strcmp( argv[i], "-boundaryType") == 0)     params.boundaryType = atoi(argv[i+1]);
            else if (strcmp( argv[i], "-saturationNW") == 0)     params.saturationNW = atof(argv[i+1]);
            else if (strcmp( argv[i], "-linkRadMin") == 0)       params.linkRadMin = atof(argv[i+1]);
            else if (strcmp( argv[i], "-linkRadMax") == 0)       params.linkRadMax = atof(argv[i+1]);
            else if (strcmp( argv[i], "-muWET") == 0)            params.muWET = atof(argv[i+1]);
            else if (strcmp( argv[i], "-muNON") == 0)            params.muNON = atof(argv[i+1]);
            else if (strcmp( argv[i], "-pressDiff") == 0)        params.pressDiff = atof(argv[i+1])*params.yLen;
            else if (strcmp( argv[i], "-Ca") == 0)               params.Ca = atof(argv[i+1]);
            else if (strcmp( argv[i], "-sTens") == 0)            params.sTens = atof(argv[i+1]);
            else if (strcmp( argv[i], "-preText") == 0)          params.preText = argv[i+1];
            else if (strcmp( argv[i], "-saveBubbles") == 0)      params.saveBubbles = atoi(argv[i+1]);
        }
    }

    return true;

}


// funciton to parse arguments from STDIN
//
// Example config.txt:
// cat config.txt | ./a.out
//
// START
// xLen 64
// yLen 64
// maxBubbles 3
// pressDiff 1.0
// saturationNW 0.50
// totalSteps 1000000
// dStep 1001
// seed 10
// bubbleMergeType 33
// flowType 0
// boundaryType 0
// saveBubbles 1
// STOP
//
// parseStdIn and parseCommandLine can be mixed 
// with parseCommandLine having priority.
bool parseStdIn(Params & params)
{
    bool found = false;

    std::string tmp;
    while (std::cin >> tmp)
    {
        if (tmp == "STOP")
            break;
        
        if (tmp == "xLen")
            std::cin >> params.xLen;
        else if (tmp == "yLen")
            std::cin >> params.yLen;
        else if (tmp == "linkRadMin")
            std::cin >> params.linkRadMin;
        else if (tmp == "linkRadMax")
            std::cin >> params.linkRadMax;
        else if (tmp == "maxBubbles")
            std::cin >> params.maxBubbles;
        else if (tmp == "pressDiff")
            std::cin >> params.pressDiff;
        else if (tmp == "saturationNW")
            std::cin >> params.saturationNW;
        else if (tmp == "totalSteps")
            std::cin >> params.totalSteps;
        else if (tmp == "dStep")
            std::cin >> params.dStep;
        else if (tmp == "seed")
            std::cin >> params.seed;
        else if (tmp == "bubbleMergeType")
            std::cin >> params.bubbleMergeType;
        else if (tmp == "flowType")
            std::cin >> params.flowType;
        else if (tmp == "boundaryType")
            std::cin >> params.boundaryType;
        else if (tmp == "saveBubbles")
            std::cin >> params.saveBubbles;
        else if (tmp == "preText")
            std::cin >> params.preText;
        else if (tmp == "Ca")
            std::cin >> params.Ca;
        else
        {
            found=false;
            break;
        }
        found = true;
    }

    return found;

}

int main(int argc, char * argv[]) {

    std::string tmp;
//    /*  //UNCOMMENT THIS
    bool found = false;
    while (std::cin >> tmp)
    {

        if (tmp == "START")
        {
            std::cout << "Found config" << std::endl;
            Params params;

            found = parseStdIn(params);
            std::cout << "seed = " << params.seed << std::endl;
            std::cout << params.xLen << std::endl;
            std::cout << params.yLen << std::endl;
            std::cout << params.dStep << std::endl;
            params.pressDiff *= params.yLen;

            parseCommandLine(argc,argv,params);


            Flow flow;
            flow.config(params);
            std::cout << "Running" << std::endl;
            flow.run();
            flow.mem_Free();
        }
        else
        {
            break;
        }
    }
    //  */ AND THIS LINE TO NOT READ FROM STDIN

/*  // COMMENT THIS LINE
    Params params;

    params.xLen = 200;
    params.yLen = 50;

    params.dStep = 2500;
    params.totalSteps = 5000000;

    params.bubbleMergeType = LINEAR_FORWARDS | MAX_BUBBLES;
//        params.bubbleMergeType = LINEAR_BACKWARDS | MAX_BUBBLES;
//        params.bubbleMergeType = CENTER_OF_MASS | MAX_BUBBLES;

    //params.flowType = FlowTypes::CONSTANT_FLOW_RATE;
    params.flowType = FlowTypes::CONSTANT_PRESSURE;

    // params.boundaryType = BoundaryTypes::OPEN_DOUBLE_NON_WETTING;
    // params.boundaryType = BoundaryTypes::OPEN_DOUBLE_WETTING;
    // params.boundaryType = BoundaryTypes::PERIODIC;
    params.boundaryType = BoundaryTypes::OPEN_SINGLE;

    params.maxBubbles = 2;
    params.seed = 10;

    params.Ca = 0.0;

    params.preText="";
    params.SAVE = true;
    params.saveBubbles = true;

    params.pressDiff = params.yLen*pressDiff;
    params.saturationNW = S;

    parseCommandLine(argc,argv,params);

    Flow flow;
    flow.config(params);
    flow.run();
    flow.mem_Free();

*/   // AND THIS LINE TO RUN JUST ONCE USING A PARAMETER OBJECT


	return 0;
}

