#include "AMGraph.h"
AMGraph::AMGraph(int vnum):vnum(vnum)
{
    g = new double*[vnum];
    for(int i=0;i<vnum;++i)
        g[i] = new double[vnum];
    vex = new Point[vnum];
    idx = new size_t[vnum];
}

/*
void AMGraph::ReadIn()
{
    int idx;
    for(int i=0;i<vnum;++i)
        cin>>idx>>vex[i];
}
*/
void AMGraph::buildGraph()
{
    for(int i=0;i<vnum;++i)
        for(int j=i;j<vnum;++j)
        {
            if(i==j)
            {
                g[i][j]=INF;
            }
            else
            {
                g[i][j]=g[j][i]=abs(vex[i],vex[j]);
            }
        }
}
void AMGraph::showAM()
{
    for(int i=0;i<vnum;++i,cout<<endl)
        for(int j=0;j<vnum;++j)
            cout<<setw(10)<<g[i][j]<<' ';
}
AMGraph * AMGraph::MST_Prim()
{
    AMGraph * T = new AMGraph(vnum);
    for(int i=0;i<T->vnum;++i)
        for(int j=i;j<T->vnum;++j)
            T->g[i][j]=T->g[j][i]=INF;
    bool inMST[vnum];
    memset(inMST,0,sizeof(inMST));
    inMST[0]=1;
    int numDone=1;
    double minCost[vnum];
    int edgeTo[vnum];
    for(int i=0;i<vnum;++i)
    {
        minCost[i]=g[0][i];
        edgeTo[i]=0;
    }
    while(numDone<vnum)
    {
        int u = -1;
        for(int i = 0;i<vnum;++i )
        {
            if(!inMST[i]&&(u==-1||minCost[i]<minCost[u]))
            {
                u=i;
            }
        }
        //cout<<"edgeTo[u]"<<u<<' '<<edgeTo[u]<<endl;
        T->g[u][edgeTo[u]]=T->g[edgeTo[u]][u]= g[u][edgeTo[u]];

        //copy_element data except AM graph
        for(size_t i=0;i<vnum;++i)
        {
            T->idx[i]=idx[i];
            T->vex[i]=vex[i];
        }
        inMST[u]=1;
        for(int i = 0;i<vnum;++i )
            if(minCost[i]>g[u][i])
            {
                minCost[i] = g[u][i];
                edgeTo[i]=u;
            }
        ++numDone;
    }
    return T;
}
void AMGraph::preOrderTraverse(int & idx, int root,int* record,bool *vis)
{
    vis[root]=1;
    record[idx++]=root;
    //cout<<'('<<idx[root]<<','<<vex[root].x<<','<<vex[root].y<<')'<<endl;
    for(int i = 0;i<vnum;++i)
        if(!vis[i]&&g[root][i]<INF)
        {
            preOrderTraverse(idx,i,record,vis);
        }
}
bool AMGraph::isHamiltonianCycle(int a[])
{
    for(size_t i=0;i<vnum-1;++i)
        if(g[a[i]][a[i+1]]==INF)
            return false;
    return g[a[vnum-1]][a[0]]!=INF;
}
double AMGraph::HamilCycleCost(int a[])
{
     double ret=0;
     for(size_t i=0;i<vnum-1;++i)
        ret+=g[a[i]][a[i+1]];
     ret+=g[a[vnum-1]][a[0]];
     return ret;
}
double AMGraph::getOnemin(size_t idx,double except)
{
    double minn=INF;
    for(size_t i=0;i<vnum;++i)
        if(except==g[idx][i])
            except=-1;//保证下次再见到except参与比较
        else if(g[idx][i]<minn)
            minn=g[idx][i];
    return minn;
}
double AMGraph::getTwomin(size_t idx)
{
    double a=INF,b=INF;
    for(size_t i=0;i<vnum;++i)
        if(g[idx][i]<a)
        {
            b=a;
            a=g[idx][i];
        }
        else if(g[idx][i]<b)
            b=g[idx][i];
    return a+b;
}
double AMGraph::lowerbound(const state & a)
{
    double ret = 0;
    for(int i=1;i<a.idx;++i)
        ret+=2*g[a.H[i-1]][a.H[i]];
    for(int i=0;i<vnum;++i)
        if(!a.vis[i])
            ret+=getTwomin(i);
    if(a.idx>=2)
    {
        ret+=getOnemin(a.H[0],g[a.H[0]][a.H[1]]);
        ret+=getOnemin(a.H[a.idx-1],g[a.H[a.idx-2]][a.H[a.idx-1]]);
    }
    return ret/2;
}
void AMGraph::reverseSegmentIfBetter(int H[],int s,int e)  //[s,e)
{
    //...a_b...c_d...  ->  ...a_c...b_d...
    int a = s-1>=0?s-1:s-1+vnum;
    int b = s;
    int c = e-1;
    int d = e%vnum;
    if(g[H[a]][H[b]]+g[H[c]][H[d]]>g[H[a]][H[c]]+g[H[b]][H[d]])
    {//reverse the segment
        for(int i=b;i<b+(c-b)/2;++i)
            swap(H[i],H[(e+s)-i-1]);
    }
}
