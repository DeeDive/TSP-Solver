#ifndef UTILITY_H
#define UTILITY_H
#include<bits/stdc++.h>
using namespace std;
//#define VNUM 8
#define INF (1<<30)-1
#define eps 1e-6
struct Edge{
    int s,e;
    double dist;
    Edge(int s,int e,double dist):s(s),e(e),dist(dist){}
    friend bool operator<(Edge &a,Edge &b);
};
struct state{
    int * H;
    int idx = 0;
    bool * vis;
    double currentCost=0;
    state(int vnum)
    {
        currentCost=0;
        H=new int[vnum];
        vis = new bool[vnum];
        memset(vis,0,sizeof(vis[0])*vnum);
    }
    state(const state & b,int vnum)
    {
        H=new int[vnum];
        idx=b.idx;
        currentCost = b.currentCost;
        vis = new bool[vnum];
        memcpy(H,b.H,sizeof(b.H[0])*vnum);
        memcpy(vis,b.vis,sizeof(b.vis[0])*vnum);
    }
    friend bool operator<(const state & a, const state & b);
};
struct Point{
    double x,y;
    friend double abs(Point & a ,Point & b);
    friend istream & operator>>(istream & in, Point & a);
};

#endif // UTILITY_H
