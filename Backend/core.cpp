#include<bits/stdc++.h>
#include"Node.h"
using namespace std;
//#define VNUM 8
#define INF 1<<30-1
#define eps 1e-6
struct Point{
    double x,y;
    friend double abs(Point & a ,Point & b);
    friend istream & operator>>(istream & in, Point & a);
};
double abs(Point & a ,Point & b)
{
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
istream & operator>>(istream & in, Point & a)
{
    in>>a.x>>a.y;
}
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
        memset(vis,0,sizeof(vis));
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
bool operator<(const state & a, const state & b)
{
    return a.currentCost>b.currentCost;
}
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
    double lowerbound(const state & a);
    double getOnemin(size_t idx, double except);
    double getTwomin(size_t idx);
};
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
    bool inMST[vnum]={0};
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
TSPSolver::TSPSolver(string filename, int graphChoice)  //0 导入文件 1 随机生成
{
    ifstream in(filename);
    if(!in)
    {
        cerr<<"Error in opening file"<<endl;
        exit(EXIT_FAILURE);
    }
    if(graphChoice==0)
    {
        cout<<"in";
        string str;
        int vnum;
        while(in>>str)
            if(str=="NODE_COORD_SECTION")
            {
                break;
            }
            else if(str=="DIMENSION")
            {
                in>>str; //':'
                in>>vnum;
                cout<<vnum;
            }
        mg = new AMGraph(vnum);

        //read in data
        for(int i=0;i<vnum;++i)
            in>>mg->idx[i]>>mg->vex[i];
        for(int i=0;i<vnum;++i)
            cout<<'('<<mg->idx[i]<<','<<mg->vex[i].x<<','<<mg->vex[i].y<<')'<<endl;

        //build AMGraph
        mg->buildGraph();
        //mg->showAM();
    }
    else if(graphChoice==1)
    {
        ;//新功能留有余地
    }
}
pair<int *,double> TSPSolver::BruteForce()
{
    int * ans = new int[mg->vnum];
    clock_t st = clock();
    int a[mg->vnum];
    double minCost = INF,cost;
    for(size_t i=0;i<mg->vnum;++i)
        a[i]=i;
    do
    {
        if(mg->isHamiltonianCycle(a)==true)
        {
            cost = mg->HamilCycleCost(a);
            if(cost<minCost)
            {
                minCost = cost;
                memcpy(ans,a,sizeof(a));
            }
        }
    }
    while(next_permutation(a+1,a+mg->vnum));
    clock_t ed = clock();
    return make_pair(ans,double(ed-st)/CLOCKS_PER_SEC);
}
pair<int *,double> TSPSolver::NearestNeighbor()
{
    clock_t st = clock();
    bool used[mg->vnum]={0};
    used[0]=1;
    int  * H = new int[mg->vnum];  //hamiltonian cycle
    int numDone=0;
    H[numDone++]=0;
    int from =0;
    while(numDone<mg->vnum)
    {
        int nn = -1; //nearest neighbor
        int nn_mincost = INF;
        for(int i=0;i<mg->vnum;++i)
           if(!used[i]&&mg->g[from][i]<nn_mincost)
           {
               nn_mincost = mg->g[from][i];
               nn = i;
           }
       used[nn]=1;
       H[numDone++]=nn;
       from = nn;
    }
    clock_t ed = clock();
    return make_pair(H,double(ed-st)/CLOCKS_PER_SEC);
}
pair<int *,double> TSPSolver::RepeatedNearestNeighbor()
{
    clock_t st = clock();
    double minC = INF;
    int * FinalH = new int [mg->vnum];
    for(int j = 0;j<mg->vnum;++j)
    {
        bool used[mg->vnum]={0};
        used[j]=1;
        int H[mg->vnum];  //hamiltonian cycle
        int numDone=0;
        H[numDone++]=j;
        int from =j;
        while(numDone<mg->vnum)
        {
            int nn = -1; //nearest neighbor
            int nn_mincost = INF;
            for(int i=0;i<mg->vnum;++i)
            if(!used[i]&&mg->g[from][i]<nn_mincost)
            {
                nn_mincost = mg->g[from][i];
                nn = i;
            }
            used[nn]=1;
            H[numDone++]=nn;
            from = nn;
        }
        double cost = mg->HamilCycleCost(H);
        if(cost<minC)
        {
            minC = cost;
            memcpy(FinalH,H,sizeof(H));
        }
    }
    clock_t ed = clock();
    return make_pair(FinalH,double(ed-st)/CLOCKS_PER_SEC);
}
pair<int *,double> TSPSolver::AlteredNearestNeighbor()
{
    static clock_t st = clock();
    static int * H = NearestNeighbor().first; //use nn to find a tour
    double original_cost = mg->HamilCycleCost(H);
    cout<<"alter"<<original_cost<<endl;

    for(int len = mg->vnum-2;len>1;--len)
        for(int start = 0;start< mg->vnum-len+1;++start)  //应该可以除以2--不行
            mg->reverseSegmentIfBetter(H,start,start+len);
    double newCost = mg->HamilCycleCost(H);
    //cout<<newCost<<' '<<original_cost<<endl;
    if(newCost<original_cost)
        return AlteredNearestNeighbor();
    else
    {
        //cout<<newCost<<' '<<original_cost<<endl;
        //assert(fabs(newCost-original_cost)<eps);
        clock_t ed = clock();
        return make_pair(H,double(ed-st)/CLOCKS_PER_SEC);
    }
}
struct Edge{
    int s,e;
    double dist;
    Edge(int s,int e,double dist):s(s),e(e),dist(dist){}
    friend bool operator<(Edge &a,Edge &b);
};
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
pair<int *,double>  TSPSolver::GreedyAlgorithm()
{
    clock_t st = clock();
    int degree[mg->vnum]={0}; //记录结点的度数
    int * H = new int[mg->vnum];
    int numDone =0;
    priority_queue<Edge> q;
    for(int i=0;i<mg->vnum;++i)
        for(int j=i+1;j<mg->vnum;++j)
            q.push(Edge(i,j,mg->g[i][j]));

    //union-find init
    int parent[mg->vnum], sz[mg->vnum];
    for(size_t i =0;i<mg->vnum;++i)
    {
        parent[i]=i;
        sz[i]=1;
    }
    vector<Edge> vec;
    while(numDone<mg->vnum-1&&!q.empty())
    {
        Edge e = q.top();
        q.pop();
        if(root(parent,e.s)!=root(parent,e.e)&&degree[e.s]<2&&degree[e.e]<2)
        {
            vec.push_back(e);
            ++degree[e.s];
            ++degree[e.e];
            unionop(parent,sz,e.e,e.s);
            assert(root(parent,e.s)==root(parent,e.e));
            ++numDone;
        }
    }
    int start;//查找起点
    int cnt[mg->vnum]={0};
    for(auto ele:vec)
    {
        //cout<<'('<<ele.e<<','<<ele.s<<')'<<' ';
        cnt[ele.s]++;
        cnt[ele.e]++;
    }
    //cout<<endl;
    for(int i=0;i<mg->vnum;++i)
        if(cnt[i]==1)
        {
            start=i;
            break;
        }
    int index=0;
    bool used[mg->vnum]={0};
    used[start]=1;
    assert(numDone==vec.size());
//    for(int i=0;i<numDone;++i)
//        cout<<'('<<vec[i].s<<','<<vec[i].e<<')';
    while(index<numDone) //n-1先加入n-1个点
    {
        for(size_t i=0;i<vec.size();++i)
        {
            if(start==vec[i].s&&!used[vec[i].e])
            {
                used[vec[i].s]=1;
                H[index++]=start;
                start=vec[i].e;
                break;
            }
            if(start ==vec[i].e&&!used[vec[i].s])
            {
                used[vec[i].e]=1;
                H[index++]=start;
                start=vec[i].s;
                break;
            }
        }
    }

    H[index++]=start; //这里不是加圈，而是把最后一个点加进来
    clock_t ed = clock();
    return make_pair(H,double(ed-st)/CLOCKS_PER_SEC);
}
pair<int *,double> TSPSolver::DynamicProgramming()
{
    clock_t st = clock();
    vector<int> ** path = new vector<int> *[1<<mg->vnum];
    double ** dp = new double*[1<<mg->vnum];
    for(int i=0;i<(1<<mg->vnum);++i)
    {
        path[i]=new vector<int>[mg->vnum];
        dp[i]= new double[mg->vnum];
    }
    int * H = new int [mg->vnum],idx=0;
    int cnt = (1 << mg->vnum); //总共的状态数
    double ans =INF;
    vector<int> ans_path;
    for(int i=0;i<cnt;++i)
        for(int j=0;j<mg->vnum;++j)
            dp[i][j]=INF;
    dp[1][0]=0;
    path[1][0].push_back(0);
    //for(int i=0;i<cnt;++i)
    //    dp[i][0]=mg->g[i][0];
    for(int i=1;i<cnt;++i)
        for(int j=1;j<mg->vnum;++j)
        {
            if(i&(1<<j))
                continue;
            if(!(i&1))
                continue;
            for(int k=0;k<mg->vnum;++k)
                if(i&(1<<k))
                {
                    if(dp[i][k]+mg->g[k][j]<dp[(1<<j)|i][j])
                    {
                        dp[(1<<j)|i][j] = dp[i][k]+mg->g[k][j];
                        path[(1 << j) | i][j] = path[i][k];
                        path[(1 << j) | i][j].push_back(j);
                    }
                    dp[(1 << j) | i][j] = min(dp[(1 << j) | i][j],
                                              dp[i][k] + mg->g[k][j]); // 转移方程
                }
        }
    for (int i = 0; i < mg->vnum; i++){
        if(dp[cnt - 1][i] + mg->g[i][0]< ans){
            ans=dp[cnt - 1][i] + mg->g[i][0];
            ans_path = path[cnt-1][i];
        }
    }
    for(int i=0;i<ans_path.size();++i)
        H[i]=ans_path[i];
    clock_t ed = clock();
    return make_pair(H,double(ed-st)/CLOCKS_PER_SEC);
}
pair<int *,double> TSPSolver::Backtracking()
{
    clock_t st = clock();
    int H[mg->vnum]={0};
    int * ans = new int[mg->vnum];
    int k=1;
    H[0]=0;
    bool used[mg->vnum]={0};
    used[0]=1;
    double minCost = INF;
    double currentCost = 0;
    while(k>=1)
    {
        //cout<<"k = "<<k;
        //cout<<"H[k] = "<<H[k];
        //for(int i =0 ;i<mg->vnum;++i)
          //  cout<<used[i]<<' ';
        //cout<<endl;
        if(H[k]!=0)  //是回溯回来的
        {
            used[H[k]]=0;//回溯法中的恢复
            currentCost-=mg->g[H[k-1]][H[k]];
        }
        H[k]++;
        while(H[k]<mg->vnum&&used[H[k]]||currentCost+mg->g[H[k-1]][H[k]]>=minCost) //剪枝
        {
            H[k]++;
        }
        if(H[k]<mg->vnum&&k==mg->vnum-1)
        {
           // cout<<"Done"<<endl;
            used[H[k]]=1;
            currentCost += mg->g[H[k-1]][H[k]];
            /*for(int i=0;i<mg->vnum;++i)
                cout<<H[i]<<' ';*/
            double cost =mg->HamilCycleCost(H);
            //cout<<"cost\&currentcost"<<cost<<' '<<currentCost+mg->g[H[k]][H[0]]<<endl;
            assert(fabs(cost-(currentCost+mg->g[H[k]][H[0]]))<eps);
            if(minCost>cost)
            {
                minCost=cost;
                memcpy(ans,H,sizeof(H));
            }
        }
        if(H[k]<mg->vnum&&k<mg->vnum-1)
        {
            currentCost+=mg->g[H[k-1]][H[k]];
            used[H[k++]]=1;
        }
        else
        {
            //cout<<"back"<<k<<' '<<H[k]<<endl;;
            used[H[k]]=0;
            currentCost-=mg->g[H[k-1]][H[k]];
            H[k--]=0;
        }
    }
    clock_t ed = clock();
    return make_pair(ans,double(ed-st)/CLOCKS_PER_SEC);
}
pair<int *,double> TSPSolver::BranchandBound()
{
    clock_t st = clock();
    int * H = new int[mg->vnum];  //用于返回最终结果
    state s(mg->vnum);
    s.H[s.idx++]=0;
    s.vis[0]=1;
    s.currentCost = mg->lowerbound(s);
    double lb=0;
    for(int i=0;i<mg->vnum;++i)
        lb+=mg->getTwomin(i);
    cout<<lb/2<<' '<<s.currentCost <<endl;
    priority_queue<state> q;
    double up=mg->HamilCycleCost(NearestNeighbor().first);  //upper bound
    cout<<"up"<<up<<endl;
    q.push(s);
    while(!q.empty())
    {
        //cout<<q.size()<<endl;
        state a = q.top();
        q.pop();
        cout<<'\n'<<a.currentCost<<' '<<up<<'\n'<<endl;
        if(a.currentCost>=up)
        {
            cout<<"##########Cut#######"<<endl;
            continue;
        }
        if(a.idx==mg->vnum)  //找到最短
        {
            memcpy(H,a.H,sizeof(a.H[0])*mg->vnum);
            up=min(up,mg->HamilCycleCost(H));
            //这里还不能返回结果，说不定队列中有更短的
            continue;
        }
        for(int i=0;i<mg->vnum;++i)
            if(!a.vis[i])  //lb<up
            {
                state b(a,mg->vnum);
                b.H[b.idx++]=i;
                b.vis[i]=1;
                b.currentCost = mg->lowerbound(b);
                if(b.currentCost<up)
                    q.push(b);
                else {
                    cout<<"##########Cut#######"<<endl;
                    continue;
                }
            }
    }
    clock_t ed = clock();
    return make_pair(H,double(ed-st)/CLOCKS_PER_SEC);
    //cout<<"NO ANSWER"; 如果比上界更小的可输出贪心的解
}
pair<int *,double> TSPSolver::Approx_TSP_Tour()
{
    clock_t st = clock();
    int * H = new int[mg->vnum];
    AMGraph * T = mg->MST_Prim();
    //T->showAM();
    bool * vis = new bool[mg->vnum]{0};
    size_t root = 0;
    int idx = 0;
    T->preOrderTraverse(idx,root,H,vis);
    clock_t ed = clock();
    return make_pair(H,double(ed-st)/CLOCKS_PER_SEC);
}
pair<vector<Node> *,double> TSPSolver::processResult(int * H)  //返回用于可视化的vector和总代价 of type double
{
    vector<Node> * ret = new vector<Node>;
    for(int i=0;i<mg->vnum;++i)
        ret->push_back(Node(mg->idx[H[i]],mg->vex[H[i]].x,mg->vex[H[i]].y));
    return make_pair(ret,mg->HamilCycleCost(H));
}
void TSPSolver::outputResult(pair<int *,double> tour)
{
    for(size_t i=0;i<mg->vnum;++i)
        cout<<tour.first[i]<<' ';
    cout<<"\nTotal cost: "<<mg->HamilCycleCost(tour.first)<<endl;
    cout<<"Total running time: "<<tour.second<<'\n'<<endl;
}
int main()
{
    TSPSolver tsp("../data/random.tsp",0);
   // TSPSolver tsp("../data/xqf131.tsp",0);  //alter alg.针对数据xqf131有一个小bug
    //TSPSolver tsp("default.tsp",0);
    pair<int *,double> tour;
    tour = tsp.BruteForce();
    tsp.outputResult(tour);
    tour = tsp.Approx_TSP_Tour();
    tsp.outputResult(tour);
    tour = tsp.NearestNeighbor();
    tsp.outputResult(tour);
    tour = tsp.RepeatedNearestNeighbor();
    tsp.outputResult(tour);
    tour = tsp.AlteredNearestNeighbor();
    tsp.outputResult(tour);
    tour = tsp.GreedyAlgorithm();
    tsp.outputResult(tour);

    //tour = tsp.DynamicProgramming();
    //tsp.outputResult(tour);
    //tour = tsp.Backtracking();
    //tsp.outputResult(tour);
    //tour = tsp.BranchandBound();
    //tsp.outputResult(tour);

    return 0;
}
