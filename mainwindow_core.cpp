#include "mainwindow.h"

#include "UI/GraphWidget.h"
#include <fstream>
#include <random>
#include <iostream>
#include <QDir>
#include <string>
#include <algorithm>
#include <vector>
#include <QSpinBox>
#include <windows.h>
#include <QCheckBox>
#include "qdebug.h"
#include <vector>
#include <ctime>
using namespace std;
void MainWindow::runApprox()
{
    pair<int *,double> tour = tsp->Approx_TSP_Tour();
    pair<vector<Node> *,double>  p= tsp->processResult(tour.first);
    vector<Node> vec=*(p.first);
    vec.push_back(vec[0]);
    std::vector<Node> showedge;
    for(size_t i = 0;i< vec.size();++i)
    {
        showedge.clear();
        for(size_t j=0;j<=i;++j)
            showedge.push_back(vec[j]);
        graphWidget_->setEdges(showedge);
        graphWidget_->viewport()->update();
        graphWidget_->repaint();
        Sleep(3000/vec.size());
    }
    result->setText("算法所求出的最短距离： "+QString::number(p.second));
    timec->setText("算法运行时间： "+QString::number(tour.second)+"s");
}
void MainWindow::runNN()
{
    pair<int *,double> tour = tsp->NearestNeighbor();
    pair<vector<Node> *,double>  p= tsp->processResult(tour.first);
    vector<Node> vec=*(p.first);
    vec.push_back(vec[0]);
    std::vector<Node> showedge;
    for(size_t i = 0;i< vec.size();++i)
    {
        showedge.clear();
        for(size_t j=0;j<=i;++j)
            showedge.push_back(vec[j]);
        graphWidget_->setEdges(showedge);
        graphWidget_->viewport()->update();
        graphWidget_->repaint();
        Sleep(3000/vec.size());
    }
    result->setText("算法所求出的最短距离： "+QString::number(p.second));
    timec->setText("算法运行时间： "+QString::number(tour.second)+"s");
}
void MainWindow::runRNN()
{
    pair<int *,double> tour = tsp->RepeatedNearestNeighbor();
    pair<vector<Node> *,double>  p= tsp->processResult(tour.first);
    vector<Node> vec=*(p.first);
    vec.push_back(vec[0]);
    std::vector<Node> showedge;
    for(size_t i = 0;i< vec.size();++i)
    {
        showedge.clear();
        for(size_t j=0;j<=i;++j)
            showedge.push_back(vec[j]);
        graphWidget_->setEdges(showedge);
        graphWidget_->viewport()->update();
        graphWidget_->repaint();
        Sleep(3000/vec.size());
    }
    result->setText("算法所求出的最短距离： "+QString::number(p.second));
    timec->setText("算法运行时间： "+QString::number(tour.second)+"s");
}
void MainWindow::runANN()
{
    pair<int *,double> tour = tsp->AlteredNearestNeighbor();
    pair<vector<Node> *,double>  p= tsp->processResult(tour.first);
    vector<Node> vec=*(p.first);
    vec.push_back(vec[0]);
    std::vector<Node> showedge;
    for(size_t i = 0;i< vec.size();++i)
    {
        showedge.clear();
        for(size_t j=0;j<=i;++j)
            showedge.push_back(vec[j]);
        graphWidget_->setEdges(showedge);
        graphWidget_->viewport()->update();
        graphWidget_->repaint();
        Sleep(3000/vec.size());
    }
    result->setText("算法所求出的最短距离： "+QString::number(p.second));
    timec->setText("算法运行时间： "+QString::number(tour.second)+"s");
}
void MainWindow::runGreedy()
{
    pair<int *,double> tour = tsp->GreedyAlgorithm();
    pair<vector<Node> *,double>  p= tsp->processResult(tour.first);
    vector<Node> vec=*(p.first);
    vec.push_back(vec[0]);
    std::vector<Node> showedge;
    for(size_t i = 0;i< vec.size();++i)
    {
        showedge.clear();
        for(size_t j=0;j<=i;++j)
            showedge.push_back(vec[j]);
        graphWidget_->setEdges(showedge);
        graphWidget_->viewport()->update();
        graphWidget_->repaint();
        Sleep(3000/vec.size());
    }
    result->setText("算法所求出的最短距离： "+QString::number(p.second));
    timec->setText("算法运行时间： "+QString::number(tour.second)+"s");
}
void MainWindow::runBF()
{
    pair<int *,double> tour = tsp->BruteForce();
    pair<vector<Node> *,double>  p= tsp->processResult(tour.first);
    vector<Node> vec=*(p.first);
    vec.push_back(vec[0]);
    std::vector<Node> showedge;
    for(size_t i = 0;i< vec.size();++i)
    {
        showedge.clear();
        for(size_t j=0;j<=i;++j)
            showedge.push_back(vec[j]);
        graphWidget_->setEdges(showedge);
        graphWidget_->viewport()->update();
        graphWidget_->repaint();
        Sleep(3000/vec.size());
    }
    result->setText("算法所求出的最短距离： "+QString::number(p.second));
    timec->setText("算法运行时间： "+QString::number(tour.second)+"s");
}
void MainWindow::runDP()
{
    pair<int *,double> tour = tsp->DynamicProgramming();
    pair<vector<Node> *,double>  p= tsp->processResult(tour.first);
    vector<Node> vec=*(p.first);
    vec.push_back(vec[0]);
    std::vector<Node> showedge;
    for(size_t i = 0;i< vec.size();++i)
    {
        showedge.clear();
        for(size_t j=0;j<=i;++j)
            showedge.push_back(vec[j]);
        graphWidget_->setEdges(showedge);
        graphWidget_->viewport()->update();
        graphWidget_->repaint();
        Sleep(3000/vec.size());
    }
    result->setText("算法所求出的最短距离： "+QString::number(p.second));
    timec->setText("算法运行时间： "+QString::number(tour.second)+"s");
}
void MainWindow::runBacktracking()
{
    pair<int *,double> tour = tsp->Backtracking();
    pair<vector<Node> *,double>  p= tsp->processResult(tour.first);
    vector<Node> vec=*(p.first);
    vec.push_back(vec[0]);
    std::vector<Node> showedge;
    for(size_t i = 0;i< vec.size();++i)
    {
        showedge.clear();
        for(size_t j=0;j<=i;++j)
            showedge.push_back(vec[j]);
        graphWidget_->setEdges(showedge);
        graphWidget_->viewport()->update();
        graphWidget_->repaint();
        Sleep(3000/vec.size());
    }
    result->setText("算法所求出的最短距离： "+QString::number(p.second));
    timec->setText("算法运行时间： "+QString::number(tour.second)+"s");
}
void MainWindow::runBranchandBound()
{
    pair<int *,double> tour = tsp->BranchandBound();
    pair<vector<Node> *,double>  p= tsp->processResult(tour.first);
    vector<Node> vec=*(p.first);
    vec.push_back(vec[0]);
    std::vector<Node> showedge;
    for(size_t i = 0;i< vec.size();++i)
    {
        showedge.clear();
        for(size_t j=0;j<=i;++j)
            showedge.push_back(vec[j]);
        graphWidget_->setEdges(showedge);
        graphWidget_->viewport()->update();
        graphWidget_->repaint();
        Sleep(3000/vec.size());
    }
    result->setText("算法所求出的最短距离： "+QString::number(p.second));
    timec->setText("算法运行时间： "+QString::number(tour.second)+"s");
}
void MainWindow::generateNodes(int choice)
{
    std::ifstream in;
    if(choice == 0 )
    {
        const double minX{ - (graphSizeX_->value() / 2.) };
        const double maxX{ graphSizeX_->value() / 2. };
        const double minY{ - (graphSizeY_->value() / 2.) };
        const double maxY{ graphSizeY_->value() / 2. };
        const unsigned count{ static_cast<unsigned>(graphItemCount_->value()) };

        static std::random_device rd;
        static std::mt19937 g(rd());
        std::uniform_real_distribution<double> distX(minX, maxX);
        std::uniform_real_distribution<double> distY(minY, maxY);

        ofstream out("data/random.tsp");
        out<<"NAME : RANDOM\nDIMENSION : "<<count<<"\nEDGE_WEIGHT_TYPE : EUC_2D\nNODE_COORD_SECTION\n";
        nodes_.clear();
        for(auto i = 0u; i < count; ++i)
        {
            double posX{ distX(g) };
            double posY{ distY(g) };
            out<<i+1<<' '<<posX<<' '<<posY<<'\n';
            nodes_.push_back(Node{i,posX, posY, "City " + std::to_string(i)});
        }
        out<<"EOF\n0\n"<<endl;
        out.close();
        scaleFactor =1;
        in.open("data/random.tsp");
//        graphWidget_->setNodes(nodes_);
//        graphWidget_->adjustSceneRect();
        tsp = new TSPSolver("data/random.tsp",0);
        ans->setText("实际最短哈密顿回路长度： N/A");
        result->setText("算法所求出的最短距离： 待求解");
        timec->setText("算法运行时间： N/A");
    }
    else if(choice == 1)
    {
        scaleFactor =100;
        tsp = new TSPSolver("data/default.tsp",0);
        in.open("data/default.tsp");
        if(!in)
            std::cerr<<"Error in opening file xqf131.tsp"<<std::endl;
    }
    else if(choice==2)
    {
        scaleFactor =100;
        tsp = new TSPSolver("data/ulysses16.tsp",0);
        //qDebug()<<QDir::currentPath();
        in.open("data/ulysses16.tsp");
        if(!in)
            std::cerr<<"Error in opening file pbl395.tsp"<<std::endl;
    }
    else if(choice==3)
    {
        scaleFactor =1;
        tsp = new TSPSolver("data/dj38.tsp",0);

        //qDebug()<<QDir::currentPath();
        in.open("data/dj38.tsp");
        if(!in)
            std::cerr<<"Error in opening file pbl395.tsp"<<std::endl;
    }
    else if(choice==4)
    {
        scaleFactor =100;
        tsp = new TSPSolver("data/xqf131.tsp",0);

        //qDebug()<<QDir::currentPath();
        in.open("data/xqf131.tsp");
        if(!in)
            std::cerr<<"Error in opening file pbl395.tsp"<<std::endl;
    }
    else if(choice==5)
    {
        scaleFactor =100;
        tsp = new TSPSolver("data/pbl395.tsp",0);

        //qDebug()<<QDir::currentPath();
        in.open("data/pbl395.tsp");
        if(!in)
            std::cerr<<"Error in opening file pbl395.tsp"<<std::endl;
    }
    std::string str;
    size_t count=0;
    while(in>>str)
    {
        std::cout<<str<<std::endl;
        if(str=="NODE_COORD_SECTION")
            break;
        else if(str=="DIMENSION")
        {
            in>>str;
            in>>count;
            std::cout<<count;
        }
    }
    nodes_.clear();
    for(auto i = 0u; i < count; ++i)
    {
        size_t idx;
        double posX;
        double posY;
        in>>idx>>posX>>posY;
        //进行缩放，最终结果相应的也需要除以缩放因子
        posX*=scaleFactor;
        posY*=scaleFactor;
        nodes_.push_back(Node{idx,posX, posY, "City " + std::to_string(i)});
    }
    in>>str;
    std::cout<<str<<std::endl;
    double optimum;
    in>>optimum;
    std::cout<<optimum<<std::endl;
    if(optimum);
    ans->setText("实际最短哈密顿回路长度： "+QString::number(optimum,10,5));
    graphWidget_->setNodes(nodes_);
    graphWidget_->adjustSceneRect();

}
