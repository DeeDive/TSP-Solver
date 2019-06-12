#ifndef AMGRAPH_H
#define AMGRAPH_H

#include<bits/stdc++.h>
#include "utility.h"
using namespace std;
struct AMGraph{
    double **g;
    int vnum;
    size_t * idx;
    Point * vex;

    AMGraph(int vnum);
   // void ReadIn();
    void buildGraph();
    void showAM();
    AMGraph * MST_Prim();
    void preOrderTraverse(int &, int root,int* record,bool *vis);
    bool isHamiltonianCycle(int a[]);  //return cost : -1 not HamiltonianCycle
    double HamilCycleCost(int a[]);
    void reverseSegmentIfBetter(int *,int s,int e);
    double lowerbound(const state &);
    double getOnemin(size_t idx, double except);
    double getTwomin(size_t idx);
};
#endif // AMGRAPH_H
