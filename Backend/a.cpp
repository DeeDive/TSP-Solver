#include<bits/stdc++.h>
using namespace std;
//#define VNUM 8
#define INF 1<<30-1
struct Point{
    double x,y;
    friend double abs(Point & a ,Point & b);
    friend istream & operator>>(istream & in, Point & a);
};
double abs(Point & a ,Point & b)
{
   // cout<<'*'<<a.x<<' '<<a.y<<' '<<b.x<<' '<<b.y<<'}'<<sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y))<<endl;
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));

}
istream & operator>>(istream & in, Point & a)
{
    in>>a.x>>a.y;
}
struct AMGraph{
    double **g;
    int vnum;
    Point * vex;

    AMGraph(int vnum);
    void ReadIn();
    void CalcuDist();
    void showAM();
};
AMGraph::AMGraph(int vnum):vnum(vnum)
{
    g = new double*[vnum];
    for(int i=0;i<vnum;++i)
        g[i] = new double[vnum];
    vex = new Point[vnum];
}
void AMGraph::ReadIn()
{
    int idx;
    for(int i=0;i<vnum;++i)
        cin>>idx>>vex[i];
}
void AMGraph::CalcuDist()
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
AMGraph * MST_Prim(AMGraph &G)
{
    AMGraph * T = new AMGraph(G.vnum);
    for(int i=0;i<T->vnum;++i)
        for(int j=i;j<T->vnum;++j)
            T->g[i][j]=T->g[j][i]=INF;
    bool inMST[G.vnum]={0};
    inMST[0]=1;
    int numDone=1;
    double minCost[G.vnum];
    int edgeTo[G.vnum];
    for(int i=0;i<G.vnum;++i)
    {
        minCost[i]=G.g[0][i];
        edgeTo[i]=0;
    }
    while(numDone<G.vnum)
    {
        int u = -1;
        for(int i = 0;i<G.vnum;++i )
        {
            if(!inMST[i]&&(u==-1||minCost[i]<minCost[u]))
            {
                u=i;
            }
        }
        //cout<<"edgeTo[u]"<<u<<' '<<edgeTo[u]<<endl;
        T->g[u][edgeTo[u]]=T->g[edgeTo[u]][u]= G.g[u][edgeTo[u]];
        inMST[u]=1;
        for(int i = 0;i<G.vnum;++i )
            if(minCost[i]>G.g[u][i])
            {
                minCost[i] = G.g[u][i];
                edgeTo[i]=u;
            }
        ++numDone;
    }
    return T;
}
int idx=-1;
void preOrderTraverse(AMGraph *T,int root,int * record,bool *vis)
{
    vis[root]=1;
    record[++idx]=root;
    for(int i = 0;i<T->vnum;++i)
        if(!vis[i]&&T->g[root][i]<INF)
        {
            preOrderTraverse(T,i,record,vis);
        }
}
int * Approx_TSP_Tour(AMGraph &G)
{
    int * ret = new int [G.vnum+1];
    AMGraph * T = MST_Prim(G);
    //T->showAM();
    bool vis[G.vnum]={0};
    int root = 0;
    preOrderTraverse(T,root,ret,vis);
    ret[G.vnum]= root;
    return ret;
}
int main()
{
    freopen("xqf131.tsp","r",stdin);
    string str;
    int vnum;
    while(cin>>str)
    {
        if(str=="NODE_COORD_SECTION")
            break;
        else if(str=="DIMENSION")
        {
            cin>>str;
            cin>>vnum;
            std::cout<<vnum;
        }
    }
    AMGraph G(vnum);
    G.ReadIn();
    //cout<<"ha"<<G.vex[2].x<<' '<<G.vex[2].y<<endl;
    G.CalcuDist();
    //G.showAM();
    int * H = Approx_TSP_Tour(G);
    for(int i=0;i<G.vnum+1;++i)
        cout<<H[i]<<' ';
    return 0;
}
