#ifndef FUNCTIONS_HPP_INCLUDED
#define FUNCTIONS_HPP_INCLUDED

#include <random>
#include "params.hpp"

enum FlowTypes
{
    CONSTANT_PRESSURE = 0,
    CONSTANT_FLOW_RATE = 1
};

enum BoundaryTypes
{
    OPEN_SINGLE = 0,              // NW left hand side
    OPEN_DOUBLE_NON_WETTING = 1,  // W band in middle
    OPEN_DOUBLE_WETTING = 2,      // NW band in middle
    PERIODIC = 3,                 // bi-periodic
    OPEN_SINGLE_OPPOSITE = 4,     // NW right hand side

};

enum bubbleMergeTypes
{
    MAX_BUBBLES = 1,
    MAX_2MENISCI = 2,  // test case, ignore
    MAX_3MENISCI = 4,  // test case, ignore
    LINEAR_BACKWARDS = 8,
    LINEAR_FORWARDS = 16,
    CENTER_OF_MASS = 32
};

class Flow
{
private:
    std::mt19937 generator;                                // mersenne twister RNG
    std::uniform_real_distribution<double> distribution;   // uniform distribution for radii

    Params params;

    int xLen, yLen, totalSteps, dStep, maxBubbles;
    int totalNodes, totalLinks, nStep;
    double linkRadMin, linkRadMax, pressDiff;
    double saturationNW, muWET, muNON, Ca, sTens;
    bool finish = false;                                    // simulation will quit if this variable becomes true.
    int seed;
    double M;
    double linkLen = 1.0;
    int timeStep = 0;

    double area;

    double volumeTL = 0.0;
    double flowTL2;

    double volumeNW, fluxTL, fluxNW, flowTL, flowNW, capA=-1, capB=-1, qRate; //qRate in mm^3/s.
    double globalPressure = -1, velMax, deltaT, totalT, flipT;
    int *linkNodeUp, *linkNodeDn,  *numberOfBubbles, *linkIsGhost;
    int **nodeNeighbour, **linkNeighbour;
    double *linkRad, *nodePressure, *linkWFirst, *nodeFluxNW, *nodeFluxTL, *mPos, *linkNWFirst;
    double *linkPc, *linkMob, *linkQ, *linkQ2, *linkAtil, *linkBtil, *linkSignPc, *linkSatu;
    double **bubbleStr, **bubbleEnd;
    double * p, * ap, * r;
    
    int totalits = 0;
    bool SAVE = false;
    bool saveBubbles = false;

    std::string preText = "";
    int bubbleMergeType = 9, flowType = 0, boundaryType = 0;

    const double globalPressureA = 32.0;    // dyne/mm^2 = 10 Pa
    const double globalPressureB = 64.0;    // dyne/mm^2 = 10 Pa

public:
    Flow();
    void run();
    void config(Params params);
    void mem_Aloc();
    void mem_Free();
    void node_Connect();
    void initSys_Sequential();
    void calc_qRate();
    void calc_Link_Mob_Pc();
    void solve_Conj_Grad(double p);
    void calc_Flow();
    void calc_Flow_const_flow();
    void calc_Flow_const_flow2();

    void bubbles_Kill(int link);           // Kill Bubble.
    void bubbles_Create(int link);         // Create Bubble.
    void bubbles_Merge(int link);          // Merge Bubbles
    void bubbles_Merge_Forwards(int link); // Merge Bubble (Old Rule).
    void bubbles_Merge_Backwards(int link);// Merge Bubble (New Rule).
    void bubbles_Merge_CM(int link);       // Merge Bubble (Keeping CM fixed).
    void bubbles_Merge_CM2(int link);      // New Merge Bubble (Keeping CM fixed).
    double bubbles_Move(int link);         // Move Bubble.
    void dynamic_Move();                   // One bubble step.

    void save_avgDist(std::string);
    void save_bblDist(std::string);        
    void save_Parameters(std::string);
    void loadChangeFile(std::string);
};

#endif // FUNCTIONS_HPP_INCLUDED
