#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "utility.h"

double abs(Point & a ,Point & b)
{
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
istream & operator>>(istream & in, Point & a)
{
    in>>a.x>>a.y;
}

bool operator<(const state & a, const state & b)
{
    return a.currentCost>b.currentCost;
}
bool operator<(const Edge &a,const Edge &b)
{
    return a.dist>b.dist;
}

//disjoint set
int root(int parent[],int i)
{
    //return i==parent[i]?i:parent[i]=root(parent[i]);
    while(parent[i]!=i)
    {
        parent[i]=parent[parent[i]];
        i=parent[i];
    }
    return i;
}
void unionop(int parent[],int sz[],int i,int j)
{
    int a = root(parent,i);
    int b = root(parent,j);
    if(a!=b)
        if(sz[a]<sz[b])
        {
            parent[a]=parent[b];
            sz[b]+=sz[a];
        }
        else
        {
            parent[b]=parent[a];
            sz[a]+=sz[b];
        }
}

#endif // FUNCTIONS_H
