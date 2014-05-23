#include <random>
#include <chrono>
#include <functional>
#include <iostream>
#include <fstream>
#include <deque>
#include "sslib.h"
#include "functions.hpp"
#include "params.hpp"

Flow::Flow()
{

}

// configure from parameter object and initialize
void Flow::config(Params params)
{
    this->params = params;
    xLen = params.xLen;
    yLen = params.yLen;
    totalSteps = params.totalSteps;
    dStep = params.dStep;
    maxBubbles = params.maxBubbles;
    seed = params.seed;
    saturationNW = params.saturationNW;
    linkRadMin = params.linkRadMin;
    linkRadMax = params.linkRadMax;
    muWET = params.muWET;
    muNON = params.muNON;
    Ca = params.Ca;
    sTens = params.sTens;
    pressDiff = params.pressDiff;

    bubbleMergeType = params.bubbleMergeType;
    flowType = params.flowType;
    boundaryType = params.boundaryType;
    preText = params.preText;

    totalNodes = xLen*yLen;
    totalLinks = 2*totalNodes;
    nStep = totalSteps/dStep;
    M = muNON/muWET;
    SAVE = params.SAVE;
    saveBubbles = params.saveBubbles;

    mem_Aloc(); node_Connect();

    generator = std::mt19937(seed);
    distribution = std::uniform_real_distribution<double>(0.0,1.0);

    initSys_Sequential();
    calc_qRate();
}


void calcStatistics(const std::deque<double> & cDeque, double & Avg, double & Std)
{
    Avg = 0.0;
    for (int i = 0; i < cDeque.size(); i++)
        Avg += cDeque[i];
    Avg /= cDeque.size();

    Std = 0.0;
    for (int i = 0; i < cDeque.size(); i++)
        Std += (cDeque[i] - Avg)*(cDeque[i] - Avg);
    Std /= cDeque.size();
    Std = sqrt(Std);
}


void Flow::run()
{
    char buf[120];

    std::string root_directory;
    std::string bubble_directory;
    std::string avg_directory;

    std::ofstream ofNeighbors;
    std::ofstream ofRadius;
    std::ofstream ofFracFlow;
    std::ofstream ofPressure;
    std::ofstream ofTestIncrease;
    std::ofstream ofTotalFlow;
    std::ofstream ofTotalFlow2;
    std::ofstream ofSaturation;
    std::ofstream ofTotalTime;

    std::deque<double> sDeque;
    std::deque<double> s2Deque;

    if (SAVE)
    {
        sprintf(buf,"%sd_%d_%dx%d_bmt%d_ft%d_bt%d_S%.2f_Ca%.4e_P%.4f_%1dB",preText.c_str(),seed,xLen,yLen,bubbleMergeType,flowType,boundaryType,saturationNW,Ca,pressDiff,maxBubbles);

        root_directory = buf;
        std::cout << root_directory << std::endl;
        bubble_directory = root_directory + "/bubbleDst";
        avg_directory = root_directory + "/avgDst";

        MKDIR(root_directory,0755);
        MKDIR(avg_directory,0755);
        if (saveBubbles)
            MKDIR(root_directory+"/bubbleDst",0755);



        ofFracFlow.open(root_directory+"/fracFlow.d");
        ofPressure.open(root_directory+"/pressure.d");
        ofTotalFlow.open(root_directory+"/totalFlow.d");
        ofTotalFlow2.open(root_directory+"/totalFlow2.d");
        ofTotalTime.open(root_directory+"/totalTime.d");
        ofTestIncrease.open(root_directory+"/testIncrease.d");
        ofSaturation.open(root_directory+"/saturation.d");

        if (saveBubbles)
        {
            ofNeighbors.open(root_directory+"/neighbors.d");
            for (int i = 0; i < totalNodes; i++)
                ofNeighbors << i << " " << nodeNeighbour[i][0] << " " << nodeNeighbour[i][1]
                            <<      " " << nodeNeighbour[i][2] << " " << nodeNeighbour[i][3]
                            <<      " " << linkNeighbour[i][0] << " " << linkNeighbour[i][1]
                            <<      " " << linkNeighbour[i][2] << " " << linkNeighbour[i][3] << std::endl;
            ofNeighbors.close();

            ofRadius.open(root_directory+"/radius.d");
            for (int i = 0; i < totalLinks; i++)
                ofRadius << i << " " << linkRad[i] << std::endl;
            ofRadius.close();
        }
        save_Parameters(root_directory);
        save_avgDist(avg_directory);

        if (saveBubbles)
            save_bblDist(bubble_directory);
    }

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    std::chrono::system_clock::time_point t2_prev = t1;

    int cTest = 0;

    double pressDiff2 = pressDiff;

    for (timeStep=1; timeStep<=totalSteps; timeStep++) {

        if (finish)
        {
            std::cout << "forcing finishing run" << std::endl;
            break;
        }
        if (flowType == 0)
            calc_Flow();
        else if (flowType == 1)
            calc_Flow_const_flow();
        else if (flowType == 2)
            calc_Flow_const_flow2();
        dynamic_Move();



        if (SAVE)
        {
            if (timeStep%dStep == 0) {

                std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
                double dt2 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
                double dt22 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t2_prev).count();
                t2_prev = t2;

                double maxsatu = -100.0;
                double totalsatu = 0;
                for (int i = 0 ; i < 2*xLen; i++)
                {
                    double satu = 0.0;
                    for (int j = 0; j < yLen; j++)
                    {
                        int k = j*xLen*2 + i;
                        for (int m = 0; m < numberOfBubbles[k]; m++)
                        {
                            totalsatu += bubbleEnd[k][m] - bubbleStr[k][m];
                            satu += bubbleEnd[k][m] - bubbleStr[k][m];
                        }

                    }
                    satu = satu/yLen;
                    if (i > 1 && i < 2*xLen-2)
                        maxsatu = std::max(satu,maxsatu);
                }


                pressDiff = pressDiff2;
                ofFracFlow   << timeStep << " " << flowNW/flowTL << std::endl;
                ofTotalFlow  << timeStep << " " << flowTL << " " << muNON*(flowTL/yLen)/sTens/area << std::endl;
                ofTotalTime  << timeStep << " " << totalT << std::endl;
                ofPressure   << timeStep << " " << globalPressure << " " << pressDiff << std::endl;
                ofSaturation << timeStep << " " << totalsatu/(2*xLen*yLen) << " " << volumeNW/volumeTL << std::endl;

                save_Parameters(root_directory);
                save_avgDist(avg_directory);

                double linkqx2 = 0.0;
                double A = 0;
                double A1 = 0;
                double A2 = 0;
                double B = 0;
                double Q = 0;
                
                if (flowType != 0)
                    pressDiff = globalPressure;
                for (int j = 0; j < totalLinks; j++)
                {
                    int x = j % (2*xLen);
                    int y = j / (2*xLen);
                    int i = j;

                    double mult = 1.0;
                    if (x%2 == 0) mult *= -1;
                    if (y%2 == 0) mult *= -1;
                    linkqx2 += mult*linkQ[j];
                    double b = -mult*linkMob[i]*linkPc[i];
                    double a = mult*linkQ[i] - b;
                    double q = mult*linkQ[i];
                    A += a;
                    B += b;
                    Q += q;
                }
                flowTL2 = linkqx2/totalLinks;

                ofTotalFlow2 << timeStep << " " << flowTL2 << " " << flowTL2*yLen << " " << A << " " << B << " " << Q << std::endl;

                if (saveBubbles)
                    save_bblDist(bubble_directory);

                sDeque.push_front(flowTL2);
                s2Deque.push_front(volumeNW/volumeTL);
                if (sDeque.size() > 200)
                    sDeque.pop_back();

                if (s2Deque.size() > 50)
                    s2Deque.pop_back();

                double sAvg, sStd;
                double s2Avg, s2Std;
                calcStatistics(sDeque,sAvg,sStd);
                calcStatistics(s2Deque,s2Avg,s2Std);

                char buf[200];
                sprintf(buf,"t%.7d T%.2e(%.2e) S%.3f(%.4f,%.4f) mS%.2f P%.2e Qt%9.2e(%9.2e,%9.2e) Ca%.4e i%d", timeStep
                                                                                                         , dt2, dt22
                                                                                                         , volumeNW/volumeTL
                                                                                                         , s2Avg, s2Std
                                                                                                         , maxsatu
                                                                                                         , (flowType != 0 ? (globalPressure)/yLen : (pressDiff)/yLen)
                                                                                                         , flowTL2
                                                                                                         , sAvg, sStd
                                                                                                         , muNON*(flowTL/yLen)/sTens/area
                                                                                                         , totalits/dStep);
                totalits = 0;
                std::cout << buf << std::endl;

                // calculate Q_\perp for different \Delta P
                if (timeStep % (dStep*10) == 0 && flowType == 0)
                {
                    
                    double first, last;
                    pressDiff = 0.25*yLen;
                    double A = 0;
                    double A1 = 0;
                    double A2 = 0;

                    double B = 0;
                    double Q = 0;
                    for (int i = 0; i < 100; i++)
                    {
                        calc_Flow();

                        double linkqx2 = 0.0;
                        A = 0;
                        A1 = 0;
                        A2 = 0;
                        B = 0;
                        Q = 0;
                        for (int j = 0; j < totalLinks; j++)
                        {
                            int x = j % (2*xLen);
                            int y = j / (2*xLen);
                            int i = j;

                            double mult = 1.0;
                            if (x%2 == 0) mult *= -1;
                            if (y%2 == 0) mult *= -1;
                            linkqx2 += mult*linkQ[j];
                            A += -mult*linkMob[i]*(nodePressure[linkNodeUp[i]] - nodePressure[linkNodeDn[i]] - linkIsGhost[i]*pressDiff);
                            A1 += mult*linkMob[i]*(linkIsGhost[i]*pressDiff);
                            A2 += -mult*linkMob[i]*(nodePressure[linkNodeUp[i]] - nodePressure[linkNodeDn[i]]);
                            B -= mult*linkMob[i]*linkPc[i];
                            Q += -mult*linkMob[i]*(nodePressure[linkNodeUp[i]] - nodePressure[linkNodeDn[i]] - linkIsGhost[i]*pressDiff + linkPc[i]);
                        }
                        flowTL2 = linkqx2/totalLinks;

                        // std::cout << "A " << A << " B " << B << " Qt " << A + B << " Q " << Q << " Q2 " << linkqx2 << " Q3 " << Q/totalLinks << " A1 " << A1 << " A2 " << A2   << std::endl;
                        if (i == 0)
                            first = A;
                        else if (i == 99)
                            last = A;
                        ofTestIncrease << cTest << " " << pressDiff << " " << flowTL/yLen/(2*xLen) << " " << flowTL2 << " " << A  << " " << B << " " << A + B << " " << A1 << " " << A2  << std::endl;

                        pressDiff += 0.25*yLen;
                        cTest++;


                        

                    }
                    // std::cout << "A = " << (last - first)/(99*0.25*yLen - 0.25*yLen) << "P + " << last - (last - first)/(99*0.25*yLen - 0.25*yLen)*(99*0.25*yLen) << ", B = " << B << std::endl;

                    pressDiff = pressDiff2;

                }


                // if (maxsatu < 0.9 || ( (fabs(totalsatu/totalLinks - avgS) < 0.0001 && stdS < 0.00001) && sDeque.size() == 10)  )
                if (  maxsatu < 0.97 || ((fabs((volumeNW/volumeTL) - s2Avg) < 0.001 || s2Std < 1e-8       // criterions to cancel prematurely
                    || fabs(sStd) > 10*fabs(sAvg) || (volumeNW/volumeTL > saturationNW*1.0*100) 
                    || fabs(sStd) < 1e-6 || fabs(sAvg) < 1e-6 ) && timeStep > 100000 && flowType == 0) )
                // if ( (maxsatu < 0.9) && timeStep > 50000)
                {
                    std::cout << "DONE premature" << std::endl;
                    break;
                }
                loadChangeFile(root_directory);
            }
        }
    }
    std::cout << "done last" << std::endl;
    if (SAVE)
    {
        ofFracFlow.close();
        ofPressure.close();
        ofTotalFlow.close();
        ofTotalFlow2.close();
        ofTestIncrease.close();
        ofTotalTime.close();
        ofSaturation.close();
    }
}

void Flow::save_avgDist(std::string fname)
{
    static char buf[120];
    static std::vector<double> linkfx(2*xLen);
    static std::vector<double> linkqx(2*xLen);
    static std::vector<double> linkqx2(2*xLen);
    static std::vector<double> linkvx(2*xLen);
    static std::vector<double> linksx(2*xLen);
    static std::vector<double> linksx2(2*xLen);
    static std::vector<double> linkareax(2*xLen);
    static std::vector<double> linkfy(yLen);
    static std::vector<double> linkqy(yLen);
    static std::vector<double> linkvy(yLen);
    static std::vector<double> linksy(yLen);
    static std::vector<double> linksy2(yLen);
    static std::vector<double> linkareay(yLen);

    std::fill(linkfx.begin(), linkfx.end(), 0.0);
    std::fill(linkqx.begin(), linkqx.end(), 0.0);
    std::fill(linkqx2.begin(), linkqx2.end(), 0.0);
    std::fill(linkvx.begin(), linkvx.end(), 0.0);
    std::fill(linksx.begin(), linksx.end(), 0.0);
    std::fill(linksx2.begin(), linksx2.end(), 0.0);
    std::fill(linkareax.begin(), linkareax.end(), 0.0);
    std::fill(linkfy.begin(), linkfy.end(), 0.0);
    std::fill(linkqy.begin(), linkqy.end(), 0.0);
    std::fill(linkvy.begin(), linkvy.end(), 0.0);
    std::fill(linksy.begin(), linksy.end(), 0.0);
    std::fill(linksy2.begin(), linksy2.end(), 0.0);
    std::fill(linkareay.begin(), linkareay.end(), 0.0);

    for (int i = 0; i < totalLinks; i++)
    {
        int x = i%(2*xLen);
        int y = i/(2*xLen);

        linksx[x]  += linkSatu[i];
        linkareax[x] += linkRad[i]*linkRad[i]*PI;
        linksx2[x] += linkSatu[i]*linkRad[i]*linkRad[i]*PI;
        linkqx[x] += linkQ[i];
        if (x%2 == 0)
            if (y%2 == 0)
                linkqx2[x] -= linkQ[i];
            else
                linkqx2[x] += linkQ[i];
        else
            if (y%2 == 0)
                linkqx2[x] += linkQ[i];
            else
                linkqx2[x] -= linkQ[i];
        linkvx[x] += linkQ[i]/(PI*linkRad[i]*linkRad[i]);
        linkfx[x] += linkQ[i]*linkSatu[i];

        linksy[y]  += linkSatu[i];
        linkareay[y] += linkRad[i]*linkRad[i]*PI;
        linksy2[y] += linkSatu[i]*linkRad[i]*linkRad[i]*PI;
        linkqy[y] += linkQ[i];
        linkvy[y] += linkQ[i]/(PI*linkRad[i]*linkRad[i]);
        linkfy[y] += linkQ[i]*linkSatu[i];
    }
    flowTL2 = linkqx2[0]/yLen;

    sprintf(buf,"%s/linkXdst_T-%08d.d", fname.c_str(), timeStep);
    std::ofstream xoutfile(buf);;
    xoutfile << "# sx sx_w qx vx qx2 fx" << std::endl;
    for (int i = 0; i < 2*xLen; i++)
        xoutfile << linksx[i]/yLen << " " << linksx2[i]/linkareax[i] << " "
                 << linkqx[i]/yLen << " " << linkvx[i]/yLen           << " " 
                 << linkqx2[i]/yLen << " " << linkfx[i]/linkqx[i]   << std::endl;
    xoutfile.close();

    sprintf(buf,"%s/linkYdst_T-%08d.d", fname.c_str(), timeStep);
    std::ofstream youtfile(buf);;
    youtfile << "# sy sy_w qy vy fy" << std::endl;
    for (int i = 0; i < yLen; i++)
        youtfile << linksy[i]/(2*xLen) << " " << linksy2[i]/linkareay[i] << " "
                 << linkqy[i]/(2*xLen) << " " << linkvy[i]/(2*xLen)      << " " 
                 << linkfy[i]/linkqy[i]       << std::endl;
    youtfile.close();

}



void Flow::mem_Aloc() {
    nodeNeighbour = init2D(totalNodes, 4, 0);
    linkNeighbour = init2D(totalNodes, 4, 0);
    linkNodeUp = init1D(totalLinks, 0);
    linkNodeDn = init1D(totalLinks, 0);
    linkIsGhost = init1D(totalLinks, 0);
    linkRad = init1D(totalLinks, 0.0);
    bubbleStr = init2D(totalLinks, maxBubbles+2, 0.0);
    bubbleEnd = init2D(totalLinks, maxBubbles+2, 0.0);
    numberOfBubbles = init1D(totalLinks, 0);
    nodePressure = init1D(totalNodes, 0.0);
    linkQ = init1D(totalLinks, 0.0);
    linkQ2 = init1D(totalLinks, 0.0);
    linkPc = init1D(totalLinks, 0.0);
    linkMob = init1D(totalLinks, 0.0);
    linkSatu = init1D(totalLinks, 0.0);
    linkAtil = init1D(totalLinks, 0.0);
    linkBtil = init1D(totalLinks, 0.0);
    linkNWFirst = init1D(totalLinks, 0.0);
    nodeFluxNW = init1D(totalNodes, 0.0);
    nodeFluxTL = init1D(totalNodes, 0.0);
    linkSignPc = init1D(totalLinks, 1.0);
    mPos = init1D((maxBubbles+4)*2, 0.0);

    p = init1D(totalNodes,0.0);
    r = init1D(totalNodes,0.0);
    ap = init1D(totalNodes,0.0);
}

void Flow::mem_Free() {
    free(linkNodeUp);
    free(linkNodeDn);
    free(linkIsGhost);
    free(linkRad);
    free(numberOfBubbles);
    free(nodePressure);
    free(linkQ);
    free(linkQ2);
    free(linkPc);
    free(linkMob);
    free(linkAtil);
    free(linkBtil);
    free(linkSatu);
    free(linkNWFirst);
    free(nodeFluxNW);
    free(nodeFluxTL);
    free(linkSignPc);
    free(mPos);
    free2D(nodeNeighbour, totalNodes);
    free2D(linkNeighbour, totalNodes);
    free2D(bubbleStr, totalLinks);
    free2D(bubbleEnd, totalLinks);

    free(p);
    free(r);
    free(ap);
}

void Flow::node_Connect() {
    for (int j = 0; j < yLen; j += 2)
    {
        for (int i = 0; i < xLen; i++)
        {
            int k = j*xLen + i;

            nodeNeighbour[k][0] = k - xLen - 1;
            nodeNeighbour[k][1] = k - xLen;
            nodeNeighbour[k][2] = k + xLen;
            nodeNeighbour[k][3] = k + xLen - 1;
            if (k < xLen)
            {
                nodeNeighbour[k][0] = nodeNeighbour[k][0] + totalNodes;
                nodeNeighbour[k][1] = nodeNeighbour[k][1] + totalNodes;
            }
            if (i == 0)
            {
                nodeNeighbour[k][0] = nodeNeighbour[k][0] + xLen;
                nodeNeighbour[k][3] = nodeNeighbour[k][3] + xLen;
            }
        }
    }

    for (int j = 1; j < yLen; j += 2)
    {
        for (int i = 0; i < xLen; i++)
        {
            int k = j*xLen + i;
            nodeNeighbour[k][0] = k - xLen;
            nodeNeighbour[k][1] = k - xLen + 1;
            nodeNeighbour[k][2] = k + xLen + 1;
            nodeNeighbour[k][3] = k + xLen;
            if (i == (xLen-1))
            {
                nodeNeighbour[k][1] = nodeNeighbour[k][1] - xLen;
                nodeNeighbour[k][2] = nodeNeighbour[k][2] - xLen;
            }
            if (k >= (yLen-1)*xLen)
            {
                nodeNeighbour[k][2] = nodeNeighbour[k][2] - totalNodes;
                nodeNeighbour[k][3] = nodeNeighbour[k][3] - totalNodes;
            }
        }
    }

    for (int i = 0; i < xLen; i++)
    {
        for (int j = 0; j < yLen; j += 2)
        {
            int k = j*xLen + i;
            linkNeighbour[k][0] = 2*k;
            linkNeighbour[k][1] = 2*k+1;
            linkNeighbour[k][2] = 2*(k+xLen)+1;
            linkNeighbour[k][3] = 2*(k+xLen);
        }
    }

    for (int i = 0; i < xLen; i++)
    {
        for (int j = 1; j < yLen; j += 2)
        {
            int k = j*xLen + i;
            linkNeighbour[k][0] = linkNeighbour[nodeNeighbour[k][0]][2];
            linkNeighbour[k][1] = linkNeighbour[nodeNeighbour[k][1]][3];
            linkNeighbour[k][2] = linkNeighbour[nodeNeighbour[k][2]][0];
            linkNeighbour[k][3] = linkNeighbour[nodeNeighbour[k][3]][1];
        }
    }

    for (int k = 0; k < totalNodes; k++) {
        linkNodeUp[linkNeighbour[k][0]] = k;
        linkNodeDn[linkNeighbour[k][0]] = nodeNeighbour[k][0];
        linkNodeUp[linkNeighbour[k][1]] = k;
        linkNodeDn[linkNeighbour[k][1]] = nodeNeighbour[k][1];
    }
    for (int k = 0; k < 2*xLen; k++)
        linkIsGhost[int(k)] = 1;
}

void Flow::initSys_Sequential() {
    auto rng = std::bind ( distribution, generator );

    for (int i = 0; i < totalNodes; i++)
        nodePressure[i] = 0;

    volumeTL = 0.0;
    for (int i = 0; i < totalLinks; i++)
    {
        linkRad[i] = linkLen * (linkRadMin + (linkRadMax - linkRadMin) * rng());
        volumeTL += PI*linkRad[i] * linkRad[i] * linkLen;
        linkNWFirst[i] = 0.0;
        numberOfBubbles[i] = 0;
        linkSignPc[i] = 1.0;
    }


    if (boundaryType == BoundaryTypes::OPEN_SINGLE || boundaryType == BoundaryTypes::PERIODIC || boundaryType == OPEN_SINGLE_OPPOSITE)
    {
        volumeNW = 0.0;
        for (int i = 0; i < (2*xLen); i++)
        // for (int i = 0; i < saturationNW*(2*xLen); i++)
        {
            for (int j = 0; j < yLen; j++)
            {
                int k = j*(2*xLen) + i;

                bubbleStr[k][0] = 0.0;
                if (boundaryType == OPEN_SINGLE_OPPOSITE)
                {
                    if (i > 2*(1 - saturationNW)*xLen - 8*cos(2*PI*j/yLen - 0*PI))
                    {
                        numberOfBubbles[k] = 1;
                        bubbleEnd[k][0] = 1.0;
                        linkSatu[k] = 1.0;
                    }    
                    else
                    {
                        numberOfBubbles[k] = 0;
                        bubbleEnd[k][0] = 0.0;
                        linkSatu[k] = 1.0;
                    }
                }
                else
                {
                    if (i < 2*saturationNW*xLen - 8*cos(2*PI*j/yLen - 0*PI))
                    {
                        numberOfBubbles[k] = 1;
                        bubbleEnd[k][0] = 1.0;
                        linkSatu[k] = 1.0;
                    }    
                    else
                    {
                        numberOfBubbles[k] = 0;
                        bubbleEnd[k][0] = 0.0;
                        linkSatu[k] = 1.0;
                    }
                }
                volumeNW += linkSatu[k]*(PI*linkRad[i]*linkRad[i]*linkLen);   
            }
        }
    }

    if (boundaryType == BoundaryTypes::OPEN_DOUBLE_WETTING)
    {
        double start = (1.0 - saturationNW)/2;
        double stop = 1 - (1.0 - saturationNW)/2;

        for (int j = 0; j < yLen; j++)
        {
            for (int i = start*(2*xLen)- 4*cos(2*PI*j/yLen - PI); i < stop*(2*xLen) + 4*cos(2*PI*j/yLen - PI); i++)
            {
                int k = j*(2*xLen) + i;

                bubbleStr[k][0] = 0.0;
                bubbleEnd[k][0] = 1.0;
                numberOfBubbles[k] = 1;
                linkSatu[k] = 1.0;
                volumeNW += linkSatu[k]*(PI*linkRad[i]*linkRad[i]*linkLen);
            }
        }
    }

    if (boundaryType == BoundaryTypes::OPEN_DOUBLE_NON_WETTING)
    {
        double start = (1.0 - saturationNW)/2;
        double stop = 1 - (1.0 - saturationNW)/2;

        for (int i = 0; i < start*(2*xLen); i++)
        {
            for (int j = 0; j < yLen; j++)
            {
                int k = j*(2*xLen) + i;

                bubbleStr[k][0] = 0.0;
                bubbleEnd[k][0] = 1.0;
                numberOfBubbles[k] = 1;
                linkSatu[k] = 1.0;
                volumeNW += linkSatu[k]*(PI*linkRad[i]*linkRad[i]*linkLen);
            }
        }

        for (int i = stop*(2*xLen); i < (2*xLen); i++)
        {
            for (int j = 0; j < yLen; j++)
            {
                int k = j*(2*xLen) + i;

                bubbleStr[k][0] = 0.0;
                bubbleEnd[k][0] = 1.0;
                numberOfBubbles[k] = 1;
                linkSatu[k] = 1.0;
                volumeNW += linkSatu[k]*(PI*linkRad[i]*linkRad[i]*linkLen);
            }
        }
    }

}


void Flow::calc_qRate()
{
    area = 0.0;
    for (int i = 0; i < totalLinks; i++)
    {
        area += PI * linkRad[i] * linkRad[i] * 1.0;
    }
    area = area/yLen;

    if (muNON > muWET)
        qRate = Ca * (sTens * area) / muNON;
    else
        qRate = Ca * (sTens * area) / muWET;

}

void Flow::calc_Link_Mob_Pc()
{
    for (int i = 0; i < totalLinks; i++)
    {
        linkPc[i] = 0.0;
        double radius = linkRad[i];
        double invRad = linkSignPc[i] / radius;
        double temp=0.0;
        for (int j = 0; j < numberOfBubbles[i]; j++)
        {
            double start = bubbleStr[i][j];
            double stop = bubbleEnd[i][j];
            temp += stop - start;
            linkPc[i] += 2.0 * sTens * invRad * (- cos(2.0 * PI * stop) + cos(2.0 * PI * start));
        }
        double muEff = muNON * temp + muWET * (1.0 - temp);
        linkMob[i] = PI * radius*radius*radius*radius / (8.0 * muEff);
    }
}

void Flow::solve_Conj_Grad(double pressure)
{
    int its = 0;

    double rpn = 0.0;
    //Initialization (Eq 32)
    for (int i = 0; i < totalNodes; i++)
    {
        p[i] = 0.0;
        for (int j = 0; j < 4; j++)
        {
            int linkNb = linkNeighbour[i][j];
            int nodeNb = nodeNeighbour[i][j];
            bool ghost = linkIsGhost[linkNb];
            if (j < 2)
                p[i] -= linkMob[linkNb] * (nodePressure[nodeNb] - nodePressure[i] - linkPc[linkNb] + pressure*ghost);
            else if (j >= 2)
                p[i] -= linkMob[linkNb] * (nodePressure[nodeNb] - nodePressure[i] + linkPc[linkNb] - pressure*ghost);
        }
        r[i] = p[i];
        rpn += r[i]*r[i];
    }
    double rpn0 = rpn;

    bool reinitialized = false;
    while(rpn/rpn0 >= ERR)
    {
        double rps = rpn;
        double am = 0.0;
        for (int i = 0; i < totalNodes; i++)
        {
            double api = 0.0;
            double pi = p[i];

            int linkNb1 = linkNeighbour[i][0];
            int linkNb2 = linkNeighbour[i][1];
            int linkNb3 = linkNeighbour[i][2];
            int linkNb4 = linkNeighbour[i][3];

            int nodeNb1 = nodeNeighbour[i][0];
            int nodeNb2 = nodeNeighbour[i][1];
            int nodeNb3 = nodeNeighbour[i][2];
            int nodeNb4 = nodeNeighbour[i][3];

            double p1 = p[nodeNb1];
            double p2 = p[nodeNb2];
            double p3 = p[nodeNb3];
            double p4 = p[nodeNb4];

            double lM1 = linkMob[linkNb1];
            double lM2 = linkMob[linkNb2];
            double lM3 = linkMob[linkNb3];
            double lM4 = linkMob[linkNb4];

            api = lM1 * (p1 - pi) + lM2 * (p2 - pi) + lM3 * (p3 - pi) + lM4 * (p4 - pi);

            ap[i] = api;
            am += pi * api;
        }
        am = rps / am;
        rpn = 0.0;
        for (int i = 0; i < totalNodes; i++)
        {
            nodePressure[i] += am * p[i];
            r[i] -= am * ap[i];
            rpn += r[i] * r[i];
        }

        double bm = rpn / rps;
        for (int i = 0; i < totalNodes; i++)
            p[i] = r[i] + bm*p[i];

        if (rpn/rpn0 < ERR && !reinitialized)
        {
            rps = rpn;
            rpn = 0.0;
            //Initialization (Eq 32)
            for (int i = 0; i < totalNodes; i++) {
                p[i] = 0.0;

                for (int j = 0; j < 4; j++) {
                    int linkNb = linkNeighbour[i][j];
                    int nodeNb = nodeNeighbour[i][j];
                    bool ghost = linkIsGhost[linkNb];
                    if (j < 2)
                        p[i] -= linkMob[linkNb] * (nodePressure[nodeNb] - nodePressure[i] - linkPc[linkNb] + pressure*ghost);
                    else if (j >= 2)
                        p[i] -= linkMob[linkNb] * (nodePressure[nodeNb] - nodePressure[i] + linkPc[linkNb] - pressure*ghost);
                }
                r[i] = p[i];
                rpn += r[i]*r[i];
            }
            reinitialized = true;
        }

        its++;

    }
    totalits += its;
    // std::cout << its << std::endl;
}

void Flow::calc_Flow()
{
    calc_Link_Mob_Pc();
    solve_Conj_Grad(pressDiff);

    flowTL = 0.0;
    for (int i = 0; i < totalLinks; i++)
    {

        linkQ[i] = -linkMob[i]*(nodePressure[linkNodeUp[i]] - nodePressure[linkNodeDn[i]] + linkPc[i] - pressDiff * linkIsGhost[i]);
        flowTL += linkQ[i];
    }
}

void Flow::calc_Flow_const_flow() {
    calc_Link_Mob_Pc();
    solve_Conj_Grad(globalPressureA);

    double qRate1 = 0.0;
    for (int i = 0; i < totalLinks; i++)
    {
        linkQ[i] = -linkMob[i]*(nodePressure[linkNodeUp[i]] - nodePressure[linkNodeDn[i]] + linkPc[i] - globalPressureA * linkIsGhost[i]);
        qRate1 = qRate1 + linkQ[i] * linkIsGhost[i];
    }

    solve_Conj_Grad(globalPressureB);

    double dP = globalPressureA - globalPressureB;
    double qRate2 = 0.0;
    for (int i = 0; i < totalLinks; i++)
    {
        linkQ2[i] = -linkMob[i]*(nodePressure[linkNodeUp[i]] - nodePressure[linkNodeDn[i]] + linkPc[i] - globalPressureB * linkIsGhost[i]);
        qRate2 = qRate2 + linkQ2[i] * linkIsGhost[i];
        linkAtil[i] = (linkQ[i] - linkQ2[i]) / dP;
    }

    capA = (qRate1 - qRate2) / dP;
    capB = qRate1 - capA * globalPressureA;
    globalPressure = qRate / capA - capB / capA;

    for (int i = 0; i < totalLinks; i++)
    {
        linkBtil[i] = linkQ[i] - linkAtil[i] * globalPressureA;
        linkQ[i] = linkAtil[i]*globalPressure + linkBtil[i];
    }
}

void Flow::calc_Flow_const_flow2() { // constant transverse flow rate
    area = 0.0;
    for (int i = 0; i < totalLinks; i++)
    {
        area += PI * linkRad[i] * linkRad[i] * 1.0;
    }
    area = area/(2*xLen);

    if (muNON > muWET)
        qRate = Ca * (sTens * area) / muNON;
    else
        qRate = Ca * (sTens * area) / muWET;

    calc_Link_Mob_Pc();
    solve_Conj_Grad(globalPressureA);

    double qRate1 = 0.0;
    for (int i = 0; i < totalLinks; i++)
    {
        int x = i % (2*xLen);
        int y = i / (2*xLen);

        double mult = 1.0;
        if (x%2 == 0) mult *= -1;
        if (y%2 == 0) mult *= -1;
        
        linkQ[i] = -linkMob[i]*(nodePressure[linkNodeUp[i]] - nodePressure[linkNodeDn[i]] + linkPc[i] - globalPressureA * linkIsGhost[i]);
        qRate1 = qRate1 + mult*linkQ[i] ;
    }
    qRate1 /= (2*xLen);

    solve_Conj_Grad(globalPressureB);

    double dP = globalPressureA - globalPressureB;
    double qRate2 = 0.0;
    for (int i = 0; i < totalLinks; i++)
    {
        int x = i % (2*xLen);
        int y = i / (2*xLen);

        double mult = 1.0;
        if (x%2 == 0) mult *= -1;
        if (y%2 == 0) mult *= -1;

        linkQ2[i] = -linkMob[i]*(nodePressure[linkNodeUp[i]] - nodePressure[linkNodeDn[i]] + linkPc[i] - globalPressureB * linkIsGhost[i]);
        qRate2 = qRate2 + mult*linkQ2[i];
        linkAtil[i] = (linkQ[i] - linkQ2[i]) / dP;
    }
    qRate2 /= (2*xLen);

    capA = (qRate1 - qRate2) / dP;
    capB = qRate1 - capA * globalPressureA;
    globalPressure = qRate / capA - capB / capA;
    // globalPressure = fabs(globalPressure);

    for (int i = 0; i < totalLinks; i++)
    {
        linkBtil[i] = linkQ[i] - linkAtil[i] * globalPressureA;
        linkQ[i] = linkAtil[i]*globalPressure + linkBtil[i];
    }
}

double Flow::bubbles_Move(int link)
{
    double area = PI * linkRad[link] * linkRad[link];
    double volNW = 0.0;
    linkSatu[link] = 0.0;
    if (linkQ[link] >= 0.0) {
        for (int i = 0; i < numberOfBubbles[link]; i++)
        {
            double x = bubbleEnd[link][i] - bubbleStr[link][i];
            linkSatu[link] += x/linkLen;
            volNW += x*area;
            bubbleStr[link][i] += linkQ[link] * deltaT / area;
            bubbleEnd[link][i] += linkQ[link] * deltaT / area;
            double lenNW = 0.0;
            if (bubbleStr[link][i]>linkLen && bubbleEnd[link][i]>linkLen)
                lenNW = bubbleEnd[link][i] - bubbleStr[link][i];
            else if (bubbleStr[link][i]<=linkLen && bubbleEnd[link][i]>linkLen)
                lenNW = bubbleEnd[link][i] - 1;
            nodeFluxNW[linkNodeUp[link]] += lenNW * area;
        }
        nodeFluxTL[linkNodeUp[link]] += linkQ[link] * deltaT;
    }
    else if (linkQ[link] < 0.0)
    {
        for (int i = 0; i < numberOfBubbles[link]; i++)
        {
            double x = bubbleEnd[link][i] - bubbleStr[link][i];
            linkSatu[link] += x/linkLen;
            volNW += x*area;
            bubbleStr[link][i] += linkQ[link] * deltaT / area;
            bubbleEnd[link][i] += linkQ[link] * deltaT / area;
            double lenNW = 0.0;
            if (bubbleStr[link][i]<0.0 && bubbleEnd[link][i]<0.0)
                lenNW = bubbleEnd[link][i] - bubbleStr[link][i];
            else if (bubbleStr[link][i]<0.0 && bubbleEnd[link][i]>=0.0)
                lenNW = -bubbleStr[link][i];
            nodeFluxNW[linkNodeDn[link]] += lenNW * area;
        }
        nodeFluxTL[linkNodeDn[link]] -= linkQ[link] * deltaT;
    }
    return volNW;
}

void Flow::bubbles_Kill(int link) {

    for (int i = 0; i < numberOfBubbles[link]; i++)
    {
        if (bubbleStr[link][i] < 0.0)
            bubbleStr[link][i] = 0.0;
        if (bubbleStr[link][i] > linkLen)
            bubbleStr[link][i] = linkLen; // ??
        if (bubbleEnd[link][i] < 0.0)
            bubbleEnd[link][i] = 0.0;
        if (bubbleEnd[link][i] > linkLen)
            bubbleEnd[link][i] = linkLen; /// happens to rest? in flux? zero length
    }

    int bubbles = numberOfBubbles[link];
    for (int i = numberOfBubbles[link]-1; i >= 0; i--)
    {
        if (fabs(bubbleEnd[link][i]-bubbleStr[link][i]) < TOL)
        {
            bubbles--;
            for (int j = i; j < bubbles; j++)
            {
                bubbleStr[link][j] = bubbleStr[link][j+1];
                bubbleEnd[link][j] = bubbleEnd[link][j+1];
            }
        }
    }
    numberOfBubbles[link] = bubbles;
}

void Flow::bubbles_Create(int link)
{
    double fracDn = 0.0, fracUp = 0.0;
    double area = PI * linkRad[link] * linkRad[link];

    if (nodeFluxTL[linkNodeDn[link]]>0)
        fracDn = nodeFluxNW[linkNodeDn[link]] / nodeFluxTL[linkNodeDn[link]];

    if (fracDn > 1.0)   /// ???????????????????????? WHAT HAPPENS TO THE REST???
        fracDn = 1.0;

    if (nodeFluxTL[linkNodeUp[link]]>0)
        fracUp = nodeFluxNW[linkNodeUp[link]] / nodeFluxTL[linkNodeUp[link]];

    if (fracUp > 1.0)
        fracUp = 1.0;

    if (fabs(linkQ[link]) > TOL)
    {
        if (linkQ[link] >= 0.0)
        {
            double str = linkNWFirst[link] * linkQ[link] * deltaT * (1.0-fracDn) / area;
            double end = str + linkQ[link] * deltaT * fracDn / area;
            for (int i = numberOfBubbles[link]-1; i >= 0; i--)  {
                bubbleStr[link][i+1] = bubbleStr[link][i];
                bubbleEnd[link][i+1] = bubbleEnd[link][i];
            }
            bubbleStr[link][0] = str;
            bubbleEnd[link][0] = end;
        }
        else if (linkQ[link] < 0.0)
        {
            double dif = linkQ[link] * deltaT * (1.0-fracUp) / area;
            double str = linkLen + linkNWFirst[link] * dif + linkQ[link] * deltaT * fracUp / area;
            double end = linkLen + linkNWFirst[link] * dif;
            bubbleStr[link][numberOfBubbles[link]] = str;
            bubbleEnd[link][numberOfBubbles[link]] = end;
        }

        // linkNWFirst[link] = 1.0 - linkNWFirst[link];
        linkNWFirst[link] = int(distribution(generator) + 0.5);

        numberOfBubbles[link] ++;
    }
}

void Flow::bubbles_Merge_CM(int link) {
    int i, j, k, iob, inb, imn, delPos;
    double x, dL, r1, r2, r12, m1, m2;
    iob = numberOfBubbles[link];
    imn = -1;
    delPos = -1;
    dL = linkLen + TOL;
    mPos[++imn] = 0.0;
    mPos[++imn] = 0.0;
    for (i=0; i<iob; i++) {
        mPos[++imn] = bubbleStr[link][i];
        mPos[++imn] = bubbleEnd[link][i];
        x = mPos[imn] - mPos[imn-1];
        if (dL>x && mPos[imn-1]>TOL && (linkLen-mPos[imn])>TOL) {
            dL = x;
            delPos = imn - 1;
        }
    }
    mPos[++imn] = linkLen;
    mPos[++imn] = linkLen;
    for (i=3;i<imn-3;i+=2) {
        x = mPos[i+1] - mPos[i];
        if (dL>x) {
            dL = x;
            delPos = i;
        }
    }
    if (delPos>-1) {
        r1 = mPos[delPos] - mPos[delPos-1];
        r2 = mPos[delPos+2] - mPos[delPos+1];
        r12 = r1+r2;
        m1 = dL*r2 / r12;
        m2 = dL*r1 / r12;
        mPos[delPos-1] += m1;
        mPos[delPos+2] -= m2;
    }
    for (i=delPos; i<imn-1; i++) mPos[i] = mPos[i+2];

    inb = 0;
    j=2; k=imn-4;
    if (delPos == 2) j = 0;
    else if (delPos == imn-3) k = imn-2;
    for (i=j;i<k;i+=2) {
        bubbleStr[link][inb] = mPos[i];
        bubbleEnd[link][inb] = mPos[i+1];
        inb++;
    }
    numberOfBubbles[link] = inb;
}

void Flow::bubbles_Merge_Backwards(int link)
{
    double mindist = 1 + TOL;
    int bubble1 = -1, bubble2 = -1;

    for (int i = 1; i < numberOfBubbles[link]; i++)
    {
        if (bubbleStr[link][i] - bubbleEnd[link][i-1] < mindist)
        {
            mindist = bubbleStr[link][i] - bubbleEnd[link][i-1];
            bubble1 = i-1;
            bubble2 = i;
        }
    }
    if (bubble1 > -1 && bubble2 > -1)
    {
        if (linkQ[link] >= 0.0)
        {
            bubbleEnd[link][bubble1] = bubbleEnd[link][bubble2];
            bubbleStr[link][bubble1] = bubbleStr[link][bubble1] + mindist;
            for (int i = bubble2; i < numberOfBubbles[link]-1; i++)
            {
                bubbleStr[link][i] = bubbleStr[link][i+1];
                bubbleEnd[link][i] = bubbleEnd[link][i+1];
            }
        }
        else if (linkQ[link] < 0.0)
        {
            bubbleStr[link][bubble2] = bubbleStr[link][bubble1];
            bubbleEnd[link][bubble2] = bubbleEnd[link][bubble2] - mindist;
            for (int i = bubble1; i < numberOfBubbles[link]-1; i++)
            {
                bubbleStr[link][i] = bubbleStr[link][i+1];
                bubbleEnd[link][i] = bubbleEnd[link][i+1];
            }
        }
        numberOfBubbles[link] --;
    }



    for (int k = 0; k < maxBubbles-1; k++)
    {
        mindist = 1.0 + TOL;
        for (int i = 1; i < numberOfBubbles[link]; i++)
        {
            if (bubbleStr[link][i] - bubbleEnd[link][i-1] < mindist)
            {
                mindist = bubbleStr[link][i] - bubbleEnd[link][i-1];
                bubble1 = i-1;
                bubble2 = i;
            }
        }


        if (bubble1 > -1 && bubble2 > -1 && mindist < 0.05)
        {
            if (linkQ[link] >= 0.0)
            {
                bubbleEnd[link][bubble1] = bubbleEnd[link][bubble2];
                bubbleStr[link][bubble1] = bubbleStr[link][bubble1] + mindist;
                for (int i = bubble2; i < numberOfBubbles[link]-1; i++)
                {
                    bubbleStr[link][i] = bubbleStr[link][i+1];
                    bubbleEnd[link][i] = bubbleEnd[link][i+1];
                }
            }
            else if (linkQ[link] < 0.0)
            {
                bubbleStr[link][bubble2] = bubbleStr[link][bubble1];
                bubbleEnd[link][bubble2] = bubbleEnd[link][bubble2] - mindist;
                for (int i = bubble1; i < numberOfBubbles[link]-1; i++)
                {
                    bubbleStr[link][i] = bubbleStr[link][i+1];
                    bubbleEnd[link][i] = bubbleEnd[link][i+1];
                }
            }
            numberOfBubbles[link] --;
        }
    }

}
void Flow::bubbles_Merge_Forwards(int link)
{
    double mindist = 1 + TOL;
    int bubble1 = -1;
    int bubble2 = -1;

    for (int i = 1; i < numberOfBubbles[link]; i++)
    {
        if (bubbleStr[link][i] - bubbleEnd[link][i-1] < mindist)
        {
            mindist = bubbleStr[link][i] - bubbleEnd[link][i-1];
            bubble1 = i-1;
            bubble2 = i;
        }
    }
    if (bubble1 > -1 && bubble2 > -1)
    {
        if (linkQ[link] >= 0.0)
        {
            bubbleEnd[link][bubble1] = bubbleEnd[link][bubble2] - mindist;
            for (int i = bubble2; i < numberOfBubbles[link]-1; i++)
            {
                bubbleStr[link][i] = bubbleStr[link][i+1];
                bubbleEnd[link][i] = bubbleEnd[link][i+1];
            }
        }
        else if (linkQ[link] < 0.0)
        {
            bubbleStr[link][bubble2] = bubbleStr[link][bubble1] + mindist;
            for (int i = bubble1; i < numberOfBubbles[link]-1; i++)
            {
                bubbleStr[link][i] = bubbleStr[link][i+1];
                bubbleEnd[link][i] = bubbleEnd[link][i+1];
            }
        }
        numberOfBubbles[link] --;
    }


    for (int k = 0; k < maxBubbles-1; k++)
    {
        mindist = 1.0 + TOL;
        for (int i = 1; i < numberOfBubbles[link]; i++)
        {
            if (bubbleStr[link][i] - bubbleEnd[link][i-1] < mindist)
            {
                mindist = bubbleStr[link][i] - bubbleEnd[link][i-1];
                bubble1 = i-1;
                bubble2 = i;
            }
        }
        if (bubble1 > -1 && bubble2 > -1 && mindist < 0.05)
        {
            if (linkQ[link] >= 0.0)
            {
                bubbleEnd[link][bubble1] = bubbleEnd[link][bubble2] - mindist;
                for (int i = bubble2; i < numberOfBubbles[link]-1; i++)
                {
                    bubbleStr[link][i] = bubbleStr[link][i+1];
                    bubbleEnd[link][i] = bubbleEnd[link][i+1];
                }
            }
            else if (linkQ[link] < 0.0)
            {
                bubbleStr[link][bubble2] = bubbleStr[link][bubble1] + mindist;
                for (int i = bubble1; i < numberOfBubbles[link]-1; i++)
                {
                    bubbleStr[link][i] = bubbleStr[link][i+1];
                    bubbleEnd[link][i] = bubbleEnd[link][i+1];
                }
            }
            numberOfBubbles[link] --;
        }
    }

}

void Flow::bubbles_Merge_CM2(int link)
{
    //Merge two nearest non-wetting bubbles.
    int i, bubble1, bubble2;
    double dist, mindist, r1, r2, r12, dL1, dL2;
    mindist = linkLen + TOL;
    bubble1 = -1;
    bubble2 = -1;

    for (i=1; i<numberOfBubbles[link]; i++) {
        dist = bubbleStr[link][i] - bubbleEnd[link][i-1];
        if (dist < mindist) {
            mindist = dist;
            bubble1 = i-1;
            bubble2 = i;
        }
    }
    if (bubble1 > -1 && bubble2 > -1) {
        r1 = bubbleEnd[link][bubble1] - bubbleStr[link][bubble1];
        r2 = bubbleEnd[link][bubble2] - bubbleStr[link][bubble2];
        r12 = r1 + r2;
        dL1 = mindist * r2 / r12;
        dL2 = mindist * r1 / r12;
        bubbleStr[link][bubble1] = bubbleStr[link][bubble1] + dL1;
        bubbleEnd[link][bubble1] = bubbleEnd[link][bubble2] - dL2;

        for(i=bubble2; i<numberOfBubbles[link]-1; i++) {
            bubbleStr[link][i] = bubbleStr[link][i+1];
            bubbleEnd[link][i] = bubbleEnd[link][i+1];
        }
        numberOfBubbles[link] --;
    }
}

void Flow::bubbles_Merge(int link)
{
    int i = link;
    while (numberOfBubbles[i] > maxBubbles)
    {
        if (bubbleMergeType & LINEAR_BACKWARDS)
            bubbles_Merge_Backwards(i);
        else if (bubbleMergeType & LINEAR_FORWARDS)
            bubbles_Merge_Forwards(i);
        else if (bubbleMergeType & CENTER_OF_MASS)
            bubbles_Merge_CM2(i);
        else
            std::cerr << "UNKNOWN BUBBLE MERGE RULE" << std::endl;
    }

    // enforce symmetric bubble merge rules
    if (numberOfBubbles[i] == maxBubbles)
    {
        if ((bubbleMergeType & MAX_3MENISCI) && (bubbleStr[i][0] > TOL && bubbleEnd[i] [maxBubbles-1] < 1 - TOL))
        {   // if both of the ends are empty (no NW menisci)
            if (bubbleMergeType & LINEAR_BACKWARDS)
                bubbles_Merge_Backwards(i);
            else if (bubbleMergeType & LINEAR_FORWARDS)
                bubbles_Merge_Forwards(i);
            else if (bubbleMergeType & CENTER_OF_MASS)
                bubbles_Merge_CM2(i);
            else
                std::cerr << "UNKNOWN BUBBLE MERGE RULE" << std::endl;
        }
        else if ((bubbleMergeType & MAX_2MENISCI) && (bubbleStr[i][0] > TOL || bubbleEnd[i] [maxBubbles-1] < 1 - TOL) )
        {   // if one of the ends are empty (one or none NW menisci)
            if (bubbleMergeType & LINEAR_BACKWARDS)
                bubbles_Merge_Backwards(i);
            else if (bubbleMergeType & LINEAR_FORWARDS)
                bubbles_Merge_Forwards(i);
            else if (bubbleMergeType & CENTER_OF_MASS)
                bubbles_Merge_CM2(i);
            else
                std::cerr << "UNKNOWN BUBBLE MERGE RULE" << std::endl;
        }
        else if (bubbleMergeType & MAX_BUBBLES)
            ;
        else
            ;
    }
}

void Flow::dynamic_Move() {
    double velMax = 0.0;
    int velMaxPos = -1;
    for (int i = 0; i < totalLinks; i++)
    {
        double area = PI * linkRad[i] * linkRad[i];
        double vel = linkQ[i] / area;
        if (velMax < fabs(vel) && numberOfBubbles[i] > 0)
        {
            velMax = fabs(vel);
            velMaxPos = i;
        }
    }
    if (velMaxPos == -1)
    {
        printf("No maximum velocity found.\n");
        throw -123;
    }
    deltaT = 0.1 * linkLen / velMax;
    totalT = totalT + deltaT;

    flowNW = 0.0;
    flowTL = 0.0;
    volumeNW = 0.0;

    for(int i=0; i<totalLinks; i++)
    {
        volumeNW = volumeNW + bubbles_Move(i);
        flowNW += linkSatu[i]*linkQ[i];
        flowTL += linkQ[i];
        bubbles_Kill(i);
    }
    for(int i=0; i<totalLinks; i++)
    {
        bubbles_Create(i);
        bubbles_Kill(i);
        bubbles_Merge(i);
    }

    fluxTL = 0.0;
    fluxNW = 0.0;
    for(int i=0; i<totalNodes; i++) {
        fluxTL = fluxTL + nodeFluxTL[i];
        fluxNW = fluxNW + nodeFluxNW[i];
        nodeFluxTL[i] = 0.0;
        nodeFluxNW[i] = 0.0;
    }


    if (boundaryType == BoundaryTypes::OPEN_SINGLE)
    {
        for (int j = 0; j < yLen; j++)
        {
            int k = j*(2*xLen) + 0;

            for (int b = 0; b < numberOfBubbles[k]; b++)
            {
                bubbleStr[k][b] = 0;
                bubbleEnd[k][b] = 0;
            }
            numberOfBubbles[k] = 1;
            bubbleStr[k][0] = 0.0;
            bubbleEnd[k][0] = 1.0;


            k = j*(2*xLen) + (2*xLen - 1);

            for (int b = 0; b < numberOfBubbles[k]; b++)
            {
                bubbleStr[k][b] = 0;
                bubbleEnd[k][b] = 0;
            }
            numberOfBubbles[k] = 0;
        }
    }

    if (boundaryType == BoundaryTypes::OPEN_SINGLE_OPPOSITE)
    {
        for (int j = 0; j < yLen; j++)
        {
            int k = j*(2*xLen) + (2*xLen - 1);

            for (int b = 0; b < numberOfBubbles[k]; b++)
            {
                bubbleStr[k][b] = 0;
                bubbleEnd[k][b] = 0;
            }
            numberOfBubbles[k] = 1;
            bubbleStr[k][0] = 0.0;
            bubbleEnd[k][0] = 1.0;


            k = j*(2*xLen) + 0;

            for (int b = 0; b < numberOfBubbles[k]; b++)
            {
                bubbleStr[k][b] = 0;
                bubbleEnd[k][b] = 0;
            }
            numberOfBubbles[k] = 0;
        }
    }


    if (boundaryType == BoundaryTypes::OPEN_DOUBLE_WETTING)
    {
        for (int j = 0; j < yLen; j++)
        {
            int k = j*(2*xLen) + 0;

            for (int b = 0; b < numberOfBubbles[k]; b++)
            {
                bubbleStr[k][b] = 0;
                bubbleEnd[k][b] = 0;
            }
            numberOfBubbles[k] = 0;

            k = j*(2*xLen) + (2*xLen - 1);

            for (int b = 0; b < numberOfBubbles[k]; b++)
            {
                bubbleStr[k][b] = 0;
                bubbleEnd[k][b] = 0;
            }
            numberOfBubbles[k] = 0;
        }
    }

    if (boundaryType == BoundaryTypes::OPEN_DOUBLE_NON_WETTING)
    {
        for (int j = 0; j < yLen; j++)
        {
            int k = j*(2*xLen) + 0;

            for (int b = 0; b < numberOfBubbles[k]; b++)
            {
                bubbleStr[k][b] = 0;
                bubbleEnd[k][b] = 0;
            }
            numberOfBubbles[k] = 1;
            bubbleStr[k][0] = 0.0;
            bubbleEnd[k][0] = 1.0;


            k = j*(2*xLen) + (2*xLen - 1);

            for (int b = 0; b < numberOfBubbles[k]; b++)
            {
                bubbleStr[k][b] = 0;
                bubbleEnd[k][b] = 0;
            }
            numberOfBubbles[k] = 1;
            bubbleStr[k][0] = 0.0;
            bubbleEnd[k][0] = 1.0;
        }
    }

    if (boundaryType == BoundaryTypes::PERIODIC)
    {

    }

}


// Menisci positions are saved as binary files, stored as follows
// 
// for each link
// flow rate (float)
// number of bubbles (char)
// linkNWFirst (char)
// (start and end position of bubble)*(number of bubbles) (double)

void Flow::save_bblDist(std::string fname) {
    char buf[120];
    sprintf(buf,"%s/bbdst_T-%08d.d", fname.c_str(), timeStep);
    std::ofstream outfile(buf,std::ios::binary);;

    int size_from_ascii = 0;
    for (int link=0; link<totalLinks; link++)
        size_from_ascii += 6 + numberOfBubbles[link]*16;

    char * binary_array_from_ascii = new char[size_from_ascii];
    int current_index = 0;

    for (int i = 0; i < totalLinks; i++)
    {
        float Q = linkQ[i];
        char * linkQ_c = (char*)&Q;
        char num_c = numberOfBubbles[i];

        binary_array_from_ascii[current_index + 0] = linkQ_c[0];
        binary_array_from_ascii[current_index + 1] = linkQ_c[1];
        binary_array_from_ascii[current_index + 2] = linkQ_c[2];
        binary_array_from_ascii[current_index + 3] = linkQ_c[3];
        binary_array_from_ascii[current_index + 4] = num_c;
        binary_array_from_ascii[current_index + 5] = (char)linkNWFirst[i];

        current_index += 6;

        if (num_c > 0)
        {
            char * bStart_c = (char*)bubbleStr[i];
            char * bEnd_c = (char*)bubbleEnd[i];

            for (int j = 0; j < num_c; j++)
            {
                for (int k = 0; k < 8; k++)
                {
                    binary_array_from_ascii[current_index + k + 16*j + 0] = bStart_c[k + 8*j];
                    binary_array_from_ascii[current_index + k + 16*j + 8] = bEnd_c[k + 8*j];
                }
            }
            current_index += num_c*16;
        }
    }

    outfile.write(binary_array_from_ascii,sizeof(char)*size_from_ascii);
    outfile.close();

    delete [] binary_array_from_ascii;


    sprintf(buf,"%s/pdst_T-%08d.d", fname.c_str(), timeStep);
    std::ofstream poutfile(buf);;
    for (int i = 0; i < totalNodes; i++)
    {
        poutfile << nodePressure[i] << "\n";
    }
    poutfile.close();

}

void Flow::save_Parameters(std::string fname) {
    std::ofstream outfile(fname+"/parameters.d");

    timeStep = std::min(timeStep,totalSteps);
    outfile << "xLen " << xLen << std::endl;
    outfile << "yLen " << yLen << std::endl;
    outfile << "linkLen " << linkLen << std::endl;
    outfile << "linkRadMin " << linkRadMin << std::endl;
    outfile << "linkRadMax " << linkRadMax << std::endl;
    outfile << "maxBubbles " << maxBubbles << std::endl;
    outfile << "muWET " << muWET << std::endl;
    outfile << "muNON " << muNON << std::endl;
    outfile << "pressDiff " << pressDiff << std::endl;
    outfile << "sTens " <<  sTens << std::endl;
    outfile << "saturationNW " << saturationNW << std::endl;
    outfile << "Ca " << Ca << std::endl;
    outfile << "Ca2 " << flowTL*muNON/(sTens*area)/yLen << std::endl;
    outfile << "totalSteps " << totalSteps << std::endl;
    outfile << "dStep " << dStep << std::endl;
    outfile << "nStep " << timeStep/dStep << std::endl;
    outfile << "timeStep " << std::min(timeStep,totalSteps) << std::endl;
    outfile << "totalT " << totalT << std::endl;
    outfile << "area " << area << std::endl;
    outfile << "seed " << seed << std::endl;
    outfile << "bubbleMergeType " << bubbleMergeType << std::endl;
    outfile << "flowType " << flowType << std::endl;
    outfile << "boundaryType " << boundaryType << std::endl;
    outfile << "saveBubbles " << saveBubbles << std::endl;
    outfile << "preText '" << preText << "'" << std::endl;

    outfile.close();

}

 
// looks for a file name "changeFile.txt" in the folder
// creating a file containg the word "RESET" in it will 
// reset the parameters to their original state.
// typing CONTINUE cancels the current simulation 
// and continues to the next one (if any).
//
// timeStep 100000 pressDiff/bubbleMergeType/flowType x
// will change pressDiff/bubbleMergeType/flowType to x at time step 100000
//
void Flow::loadChangeFile(std::string fname)
{
    try
    {
        std::ifstream infile(fname+"/changeFile.txt");

        std::string tmp;

        while (infile >> tmp)
        {
            if (tmp == "RESET")
            {
                std::cout << "RESET" << std::endl;
                pressDiff = params.pressDiff;

                bubbleMergeType = params.bubbleMergeType;
                flowType = params.flowType;
                std::ofstream outfile(fname+"/changeFile.txt");
                outfile.close();
            }
            else if (tmp == "CONTINUE")
            {
                std::cout << "CONTINUE" << std::endl;

                finish = true;

                std::ofstream outfile(fname+"/changeFile.txt");
                outfile.close();
            }
            else if (tmp == "timeStep")
            {
                std::cout << "timeStep" << std::endl;
                int ts;
                infile >> ts;
                if (ts == timeStep)
                {
                    infile >> tmp;
                    std::cout << tmp << " changed from ";
                    if (tmp == "pressDiff")
                    {
                        std::cout << pressDiff << " to ";
                        infile >> pressDiff;
                        pressDiff *= yLen;
                        std::cout << pressDiff;
                    }
                    else if (tmp == "bubbleMergeType")
                    {
                        std::cout << bubbleMergeType << " to ";
                        infile >> bubbleMergeType;
                        std::cout << flowType;
                    }
                    else if (tmp == "flowType")
                    {
                        std::cout << flowType << " to ";
                        infile >> flowType;
                        std::cout << flowType;
                    }
                    std::cout << std::endl;

                    std::ofstream outfile(fname+"/changeFile.txt");
                    outfile.close();
                }
                    
                break;
            }
        }
    }
    catch (...)
    {
        std::cout << "ERROR" << std::endl;
    }

}