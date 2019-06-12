#ifndef TSPSOLVER_H
#define TSPSOLVER_H

#include<bits/stdc++.h>
#include"AMGraph.h"
#include "Node.h"
#include "utility.h"
using namespace std;
class TSPSolver{
private:
    AMGraph * mg;
    //ALGraph * lg;
    string filename;
    int graphChoice;  //0 -> Adjacency Matrix : 1 -> Adjacency List
public:
    TSPSolver(string filename, int graphChoice = 0);
    pair<int *,double> BruteForce();  //蛮力法，返回tour顶点下标和消耗时间
    pair<int *,double> NearestNeighbor(); //最近邻点
    pair<int *,double> RepeatedNearestNeighbor(); //不同起点最近邻点
    pair<int *,double> AlteredNearestNeighbor(); //反转片段的最近邻点算法
    pair<int *,double> GreedyAlgorithm();  //最短链接算法
    pair<int *,double> DynamicProgramming(); // 动态规划
    pair<int *,double> Backtracking(); //
    pair<int *,double> BranchandBound(); //
    pair<int *,double> Approx_TSP_Tour();  //2-近似算法
    pair<vector<Node> *,double> processResult(int *);  //返回用于可视化的vector和总代价 of type double
    void outputResult(pair<int *,double>);
};
#endif // TSPSOLVER_H
