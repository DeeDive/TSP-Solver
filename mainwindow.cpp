#include "MainWindow.h"
#include "UI/GraphWidget.h"
#include <string>
#include <qlayout.h>
#include <QSpinBox>
#include <QPushButton>
#include <QGroupBox>
#include <QSlider>
#include <QLabel>
#include <QSplitter>
#include <QAction>
#include <QMenu>
#include <QToolBar>
#include <QApplication>
#include <qevent.h>
#include <QMenuBar>
#include <QMessageBox>
#include <QTabWidget>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QSettings>
#include <QCheckBox>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow{parent}, settings_{"settings.ini", QSettings::IniFormat}
{
   createActions_();
  // createToolBars_();
   //createMenus_();
   createWidgets_();

   generateNodes(0);
   loadSettings();
}



// EXIT / SETTINGS

void MainWindow::quit()
{
    //saveSettings();
    qApp->quit();
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    //saveSettings();
    event->accept();
}

// Settings constants.
constexpr auto window_state_key = "MainWindow/state";
constexpr auto window_geometry_key = "MainWindow/geometry";
constexpr auto toolbar_file_key = "MainWindow/toolbars/file/visible";
constexpr auto toolbar_edit_key = "MainWindow/toolbars/edit/visible";
constexpr auto toolbar_run_key = "MainWindow/toolbars/run/visible";

void MainWindow::saveSettings()
{
    settings_.setValue(window_state_key, saveState());
    settings_.setValue(window_geometry_key, saveGeometry());

    settings_.setValue(toolbar_file_key, toolbars_.file->isVisible());
    settings_.setValue(toolbar_edit_key, toolbars_.edit->isVisible());
    settings_.setValue(toolbar_run_key, toolbars_.run->isVisible());
}

void MainWindow::loadSettings()
{
//    if(!settings_.value(window_state_key).isNull())
//        restoreState(settings_.value(window_state_key).toByteArray());
//    if(!settings_.value(window_geometry_key).isNull())
//        restoreGeometry(settings_.value(window_geometry_key).toByteArray());
//    else
        resize(1500, 800); // Default size.

    // Toolbars
//    toolbars_.file->setVisible(settings_.value(toolbar_file_key, true).toBool());
//    actions_.showTBs.file->setChecked(settings_.value(toolbar_file_key, true).toBool());
//    toolbars_.edit->setVisible(settings_.value(toolbar_edit_key, true).toBool());
//    actions_.showTBs.edit->setChecked(settings_.value(toolbar_edit_key, true).toBool());
//    toolbars_.run->setVisible(settings_.value(toolbar_run_key, true).toBool());
//    actions_.showTBs.run->setChecked(settings_.value(toolbar_run_key, true).toBool());

}

// SETUP.

void MainWindow::createActions_()
{
    // File
        actions_.newSimulation = new QAction(QIcon(":/img/file"), tr("&New Simulation"), this);
            actions_.newSimulation->setShortcut(QKeySequence::New);
        actions_.save = new QAction(QIcon(":/img/save"), tr("&Save"), this);
            actions_.save->setShortcut(QKeySequence::Save);
        actions_.saveAs = new QAction(QIcon(":/img/saveAs"), tr("Save &as ..."), this);
            actions_.saveAs->setShortcut(QKeySequence::SaveAs);
        actions_.quit = new QAction(QIcon(":/img/quit"), tr("E&xit"), this);
            actions_.quit->setShortcut(QKeySequence::Quit);
            connect(actions_.quit, &QAction::triggered, this, &MainWindow::quit);

    // Edit
        actions_.generateNodes = new QAction(QIcon(":/img/generate"), tr("&Generate new graph"), this);
            actions_.generateNodes->setShortcut(QKeySequence("Ctrl+G"));
            connect(actions_.generateNodes, &QAction::triggered, this, [&]() {

               generateNodes(combo->currentIndex());
            });

    // Run
            actions_.start = new QAction(QIcon(":/img/run"), tr("&Run simulation"), this);
                actions_.start->setShortcut(QKeySequence("Ctrl+R"));
               /* connect(actions_.start, &QAction::triggered, this, [&]() {
                    runGA();
                 });  */     // actions_.start = new QAction(QIcon(":/img/run"), tr("&Run simulation"), this);
//                connect(actions_.runapprox, &QAction::triggered, this, [&]() {
//                    runApprox();
//                 });        actions_.start = new QAction(QIcon(":/img/run"), tr("&Run simulation"), this);
//                actions_.start->setShortcut(QKeySequence("Ctrl+R"));
//                connect(actions_.start, &QAction::triggered, this, [&]() {
//                    runGA();
//                 });        actions_.start = new QAction(QIcon(":/img/run"), tr("&Run simulation"), this);
//                actions_.start->setShortcut(QKeySequence("Ctrl+R"));
//                connect(actions_.start, &QAction::triggered, this, [&]() {
//                    runGA();
//                 });
        actions_.startSBS = new QAction(QIcon(":/img/runSBS"), tr("&Run one step"), this);
            actions_.startSBS->setShortcut(QKeySequence("Ctrl+P"));
        actions_.stop = new QAction(QIcon(":/img/stop"), tr("&Stop simulation"), this);
            actions_.stop->setShortcut(QKeySequence("Ctrl+H"));
            // connect.
        actions_.clear = new QAction(QIcon(":/img/delete"), tr("&Clear"), this);
            connect(actions_.clear, &QAction::triggered, this, [&](){
                graphWidget_->clearEdges();
            });
        actions_.mutate = new QAction(QIcon(":/img/mutate"), tr("&Mutate"), this);
            // connect.
        actions_.crossover = new QAction(QIcon(":/img/crossover"), tr("Cr&ossover"), this);
            // connect.

    // View
    actions_.showTBs.file = new QAction(QIcon(":/DEBUG/ico/purple"), tr("Show &file toolBar"), this);
                actions_.showTBs.file->setCheckable(true);
        connect(actions_.showTBs.file, &QAction::triggered, this, [=](bool checked) {
                                                                            toolbars_.file->setVisible(checked);
                                                                            });
    actions_.showTBs.edit = new QAction(QIcon(":/DEBUG/ico/blue"), tr("Show &edit toolBar"), this);
        actions_.showTBs.edit->setCheckable(true);
        connect(actions_.showTBs.edit, &QAction::triggered, this, [=](bool checked) {
                                                                            toolbars_.edit->setVisible(checked);
                                                                            });
    actions_.showTBs.run = new QAction(QIcon(":/DEBUG/ico/green"), tr("Show &run toolBar"), this);
        actions_.showTBs.run->setCheckable(true);
        connect(actions_.showTBs.run, &QAction::triggered, this, [=](bool checked) {
                                                                            toolbars_.run->setVisible(checked);
                                                                            });
    // Help
    actions_.help = new QAction(QIcon(":/img/help"), tr("&Help"), this);
        actions_.help->setShortcut(QKeySequence::HelpContents);

    // Not implemented.
    actions_.newSimulation->setEnabled(false);
    actions_.saveAs->setEnabled(false);
    actions_.save->setEnabled(false);
    actions_.startSBS->setEnabled(false);
    actions_.stop->setEnabled(false);
    actions_.mutate->setEnabled(false);
    actions_.crossover->setEnabled(false);
}

void MainWindow::createToolBars_()
{
    // File toolbar.
    toolbars_.file = addToolBar("File");
    toolbars_.file->setObjectName("ToolBarFile");
    toolbars_.file->setFloatable(true);

    toolbars_.file->addAction(actions_.newSimulation);
    toolbars_.file->addSeparator();
    toolbars_.file->addAction(actions_.save);
    toolbars_.file->addAction(actions_.saveAs);
    toolbars_.file->addSeparator();
    toolbars_.file->addAction(actions_.quit);

    // Edit toolbar.
    toolbars_.edit = addToolBar("Edit");
    toolbars_.edit->setObjectName("ToolBarEdit");
    toolbars_.edit->addAction(actions_.generateNodes);
    toolbars_.edit->setFloatable(true);

    // Run toolbar.
    toolbars_.run = addToolBar("Run");
    toolbars_.run->setObjectName("ToolBarRun");
    toolbars_.run->setFloatable(true);

    toolbars_.run->addAction(actions_.start);
    toolbars_.run->addAction(actions_.startSBS);
    toolbars_.run->addAction(actions_.stop);
    toolbars_.run->addSeparator();
    toolbars_.run->addAction(actions_.mutate);
    toolbars_.run->addAction(actions_.crossover);
    toolbars_.run->addAction(actions_.clear);
}

void MainWindow::createMenus_()
{
    // File menu.
    menus_.file = new QMenu(tr("&File"), this);
    menuBar()->addMenu(menus_.file);

    menus_.file->addAction(actions_.newSimulation);
    menus_.file->addSeparator();
    menus_.file->addAction(actions_.save);
    menus_.file->addAction(actions_.saveAs);
    menus_.file->addSeparator();
    menus_.file->addAction(actions_.quit);

    // Edit menu.
    menus_.edit = new QMenu(tr("&Edit"), this);
    menus_.edit->addAction(actions_.generateNodes);
    menuBar()->addMenu(menus_.edit);

    // Run menu.
    menus_.run = new QMenu(tr("&Run"), this);
    menuBar()->addMenu(menus_.run);

    menus_.run->addAction(actions_.start);
    menus_.run->addAction(actions_.startSBS);
    menus_.run->addAction(actions_.stop);
    menus_.run->addSeparator();
    menus_.run->addAction(actions_.mutate);
    menus_.run->addAction(actions_.crossover);
    menus_.run->addAction(actions_.clear);

    // Display menu.
    menus_.view = new QMenu(tr("&View"), this);
    menuBar()->addMenu(menus_.view);

        QMenu* menusDisplay_toolbars{ new QMenu(tr("ToolBars"), this) };
            menusDisplay_toolbars->addAction(actions_.showTBs.file);
            menusDisplay_toolbars->addAction(actions_.showTBs.edit);
            menusDisplay_toolbars->addAction(actions_.showTBs.run);

    menus_.view->addMenu(menusDisplay_toolbars);

    // Help menu.
    menus_.help = new QMenu(tr("&Help"), this);
    menuBar()->addMenu(menus_.help);

    menus_.help->addAction(actions_.help);
    menus_.help->addSeparator();
    menus_.help->addAction(QIcon(":/img/qt"), tr("About &Qt"), qApp, SLOT(aboutQt()));
}

#include <QStandardItemModel>

void util_disableItemAt(QStandardItemModel* model, int index)
{
   QStandardItem* item = model->item(index);
   item->setFlags(item->flags() & ~Qt::ItemIsEnabled);
}

void MainWindow::createWidgets_()
{
    // 左边部分.
    graphWidget_ = new GraphWidget(this);

    // 右边功能.

        // Genetic Algorithm.
            // Parameters.
                // Population.
                QLabel* populationLabel = new QLabel(tr("Population"), this);
                populationSpinBox_ = new QSpinBox(this);
                    populationSpinBox_->setRange(5, 100000);
                    populationSpinBox_->setValue(50);
                QHBoxLayout* populationHLayout = new QHBoxLayout;
                    populationHLayout->addWidget(populationLabel);
                    populationHLayout->addWidget(populationSpinBox_);

                // Generation count.
                QLabel* generationCountLabel = new QLabel(tr("Generations count"), this);
                generationCount_ = new QSpinBox(this);
                    generationCount_->setRange(5, 100000);
                    generationCount_->setValue(100);
                QHBoxLayout* generationHLayout = new QHBoxLayout;
                    generationHLayout->addWidget(generationCountLabel);
                    generationHLayout->addWidget(generationCount_);

                // Mutation.
                    // Mutation checkbox.
                    QLabel* enableMutationLabel = new QLabel(tr("Enable mutations"), this);
                    enableMutationCheckBox_ = new QCheckBox(this);
                        enableMutationCheckBox_->setChecked(true);
                    QHBoxLayout* enableMutationLayout = new QHBoxLayout;
                        enableMutationLayout->addWidget(enableMutationLabel);
                        enableMutationLayout->addWidget(enableMutationCheckBox_);

                    // Mutation percent.
                    QLabel* mutationPercentLabel = new QLabel(tr("Mutation rate"), this);
                    mutationPercentSpinBox_ = new QSpinBox(this);
                        mutationPercentSpinBox_->setRange(0, 100);
                        mutationPercentSpinBox_->setSuffix("%");
                        mutationPercentSpinBox_->setValue(15);
                        mutationPercentSpinBox_->setSingleStep(5);
                    QHBoxLayout* mutationPercentHLayout = new QHBoxLayout;
                        mutationPercentHLayout->addWidget(mutationPercentLabel);
                        mutationPercentHLayout->addWidget(mutationPercentSpinBox_);

                    // Connect.
                    connect(enableMutationCheckBox_, &QCheckBox::toggled, this, [&](bool checked) {
                       mutationPercentSpinBox_->setEnabled(checked);
                    });

                QVBoxLayout* mutationVLayout = new QVBoxLayout;
                    mutationVLayout->addLayout(enableMutationLayout);
                    mutationVLayout->addLayout(mutationPercentHLayout);
                QGroupBox* mutationGroupBox = new QGroupBox(tr("Mutation"), this);
                mutationGroupBox->setLayout(mutationVLayout);

                // Crossover
                    // Crossover checkbox.
                    QLabel* enableCrossoverLabel = new QLabel(tr("Enable crossovers"), this);
                    enableCrossoverCheckBox_ = new QCheckBox(this);
                        enableCrossoverCheckBox_->setChecked(true);
                    QHBoxLayout* enableCrossoverLayout = new QHBoxLayout;
                        enableCrossoverLayout->addWidget(enableCrossoverLabel);
                        enableCrossoverLayout->addWidget(enableCrossoverCheckBox_);

                    // Crossover percent.
                    QLabel* crossoverPercentLabel = new QLabel(tr("Crossover rate"), this);
                    crossoverPercentSpinBox_ = new QSpinBox(this);
                        crossoverPercentSpinBox_->setRange(0, 100);
                        crossoverPercentSpinBox_->setSuffix("%");
                        crossoverPercentSpinBox_->setValue(85);
                        crossoverPercentSpinBox_->setSingleStep(5);
                    QHBoxLayout* crossoverPercentHLayout = new QHBoxLayout;
                        crossoverPercentHLayout->addWidget(crossoverPercentLabel);
                        crossoverPercentHLayout->addWidget(crossoverPercentSpinBox_);

                    // Mutation branch count.
                    QLabel* selectionModeLabel = new QLabel(tr("Selection mode"), this);
                    selectionMode_ = new QComboBox(this);
                        selectionMode_->addItem(tr("Game of death"));
                        selectionMode_->addItem(tr("Rank"));
                        selectionMode_->addItem(tr("Tournament 5"));
                        selectionMode_->addItem(tr("Tournament 10"));
                        selectionMode_->addItem(tr("Tournament 25"));
                        // Disable not implemented features.
                            selectionMode_->setCurrentIndex(2);
                            util_disableItemAt(qobject_cast<QStandardItemModel*>(selectionMode_->model()), 0);
                            util_disableItemAt(qobject_cast<QStandardItemModel*>(selectionMode_->model()), 1);
                            util_disableItemAt(qobject_cast<QStandardItemModel*>(selectionMode_->model()), 3);
                            util_disableItemAt(qobject_cast<QStandardItemModel*>(selectionMode_->model()), 4);
                    QHBoxLayout* selectionModeLayout = new QHBoxLayout;
                        selectionModeLayout->addWidget(selectionModeLabel);
                        selectionModeLayout->addWidget(selectionMode_);

                    // Connect.
                    connect(enableCrossoverCheckBox_, &QCheckBox::toggled, this, [&](bool checked) {
                       crossoverPercentSpinBox_->setEnabled(checked);
                    });

                QVBoxLayout* crossoverVLayout = new QVBoxLayout;
                    crossoverVLayout->addLayout(enableCrossoverLayout);
                    crossoverVLayout->addLayout(crossoverPercentHLayout);
                    crossoverVLayout->addLayout(selectionModeLayout);
                QGroupBox* crossoverGroupBox = new QGroupBox(tr("Crossover"), this);
                crossoverGroupBox->setLayout(crossoverVLayout);

                // Elitism.
                    // Elitism checkbox.
                    QLabel* enableElitismLabel = new QLabel(tr("Enable elitism"), this);
                    enableElitismCheckBox_ = new QCheckBox(this);
                        enableElitismCheckBox_->setChecked(true);
                    QHBoxLayout* enableElitismLayout = new QHBoxLayout;
                        enableElitismLayout->addWidget(enableElitismLabel);
                        enableElitismLayout->addWidget(enableElitismCheckBox_);

                    // Elitism crossover.
                    QLabel* elitismCrossoverLabel = new QLabel(tr("Crossover elitism"), this);
                    elitismCrossover_ = new QSpinBox(this);
                        elitismCrossover_->setValue(5);
                        elitismCrossover_->setMinimum(0);
                        elitismCrossover_->setMaximum(populationSpinBox_->value());
                    QHBoxLayout* elitismCrossoverLayout = new QHBoxLayout;
                        elitismCrossoverLayout->addWidget(elitismCrossoverLabel);
                        elitismCrossoverLayout->addWidget(elitismCrossover_);

                    // Elitism keeped.
                    QLabel* elitismKeepedLabel = new QLabel(tr("Mutation elitism"), this);
                    elitismKeeped_ = new QSpinBox(this);
                        elitismKeeped_->setValue(10);
                        elitismKeeped_->setMinimum(0);
                        elitismKeeped_->setMaximum(populationSpinBox_->value());
                    QHBoxLayout* elitismKeepedLayout = new QHBoxLayout;
                        elitismKeepedLayout->addWidget(elitismKeepedLabel);
                        elitismKeepedLayout->addWidget(elitismKeeped_);

                    // Connect.
                    connect(enableElitismCheckBox_, &QCheckBox::toggled, this, [&](bool checked) {
                       elitismCrossover_->setEnabled(checked);
                       elitismKeeped_->setEnabled(checked);
                    });
                    connect(populationSpinBox_, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged), this, [&](int value) {
                        elitismCrossover_->setMaximum(value);
                        elitismKeeped_->setMaximum(value);
                    });


                QVBoxLayout* elitismVLayout = new QVBoxLayout;
                    elitismVLayout->addLayout(enableElitismLayout);
                    elitismVLayout->addLayout(elitismCrossoverLayout);
                    elitismVLayout->addLayout(elitismKeepedLayout);
                QGroupBox* elitismGroupBox = new QGroupBox(tr("Elitism"), this);
                elitismGroupBox->setLayout(elitismVLayout);

            QVBoxLayout* parameterVLayout = new QVBoxLayout;
                parameterVLayout->addLayout(populationHLayout);
                parameterVLayout->addLayout(generationHLayout);
                parameterVLayout->addWidget(mutationGroupBox);
                parameterVLayout->addWidget(crossoverGroupBox);
                parameterVLayout->addWidget(elitismGroupBox);

            startGAButton_ = new QPushButton(tr("&Start"), this);
                connect(startGAButton_, &QPushButton::pressed, this, [&]() {
                    result->setText("算法所求出的最短距离： N/A");
                    timec->setText("算法运行时间： N/A");

                    //ans->setText("实际最短路线长度：");
                    //runGA();
                 });

            QVBoxLayout* gaMainLayout = new QVBoxLayout;
                gaMainLayout->addLayout(parameterVLayout);
                gaMainLayout->addWidget(startGAButton_);

            QWidget* gaWidget = new QWidget(this);
                gaWidget->setLayout(gaMainLayout);
        //精确求解accurWidget
                accurbtn1 = new QPushButton;
                accurbtn1->setText("蛮力算法");
                accurbtn2 = new QPushButton;
                accurbtn2->setText("动态规划算法");
                accurbtn3 = new QPushButton;
                accurbtn3->setText("回溯算法");
                accurbtn4 = new QPushButton;
                accurbtn4->setText("分支限界算法");

                QVBoxLayout* accurlay = new QVBoxLayout;
                accurlay->addWidget(accurbtn1);
                accurlay->addWidget(accurbtn2);
                accurlay->addWidget(accurbtn3);
                accurlay->addWidget(accurbtn4);
//                connect(combo,SIGNAL(accurbtn1::clicked()),this,SLOT(runApprox()));
                connect(accurbtn1, &QPushButton::pressed, this, [&]() {
                    runBF();
                 });
                connect(accurbtn2, &QPushButton::pressed, this, [&]() {
                    runDP();
                 });
                connect(accurbtn3, &QPushButton::pressed, this, [&]() {
                    runBacktracking();
                 });
                connect(accurbtn4, &QPushButton::pressed, this, [&]() {
                    runBranchandBound();
                 });
                QWidget* accurWidget = new QWidget(this);
                    accurWidget->setLayout(accurlay);
                    accurWidget->setEnabled(true);
        //近似求解
          // QLabel* approx = new QLabel("近似求解",this);
            approxbtn1 = new QPushButton;
           approxbtn1->setText("2-近似算法");
           approxbtn1->resize(30,50);
                           approxbtn2 = new QPushButton;
                           approxbtn2->setText("最近邻点算法");
                           approxbtn3 = new QPushButton;
                           approxbtn3->setText("重复最近邻点算法");
                           approxbtn4 = new QPushButton;
                           approxbtn4->setText("变换最近邻点算法");
                           approxbtn5 = new QPushButton;
                           approxbtn5->setText("最短链接算法");
           QVBoxLayout* approxlay = new QVBoxLayout;
           approxlay->addWidget(approxbtn1);
           approxlay->addWidget(approxbtn2);
           approxlay->addWidget(approxbtn3);
           approxlay->addWidget(approxbtn4);
           approxlay->addWidget(approxbtn5);
           connect(approxbtn1, &QPushButton::pressed, this, [&]() {
               runApprox();
            });
           connect(approxbtn2, &QPushButton::pressed, this, [&]() {
               runNN();
            });
           connect(approxbtn3, &QPushButton::pressed, this, [&]() {
               runRNN();
            });
           connect(approxbtn4, &QPushButton::pressed, this, [&]() {
               runANN();
            });
           connect(approxbtn5, &QPushButton::pressed, this, [&]() {
               runGreedy();
            });
           QWidget* approxWidget = new QWidget(this);
               approxWidget->setLayout(approxlay);
               approxWidget->setEnabled(true);
        // Ant Colony Algorithm.
            // Population.
            QLabel* antCountLabel = new QLabel(tr("Ant count"), this);
            antCount_ = new QSpinBox(this);
                antCount_->setRange(5, 100000);
                antCount_->setValue(200);
            QHBoxLayout* antCountLayout = new QHBoxLayout;
                antCountLayout->addWidget(antCountLabel);
                antCountLayout->addWidget(antCount_);

            // Evaporation.
            QLabel* evapRateLabel = new QLabel(tr("Evaporation rate"), this);
            evaporationRate_ = new QSpinBox(this);
                mutationPercentSpinBox_->setRange(0, 100);
                mutationPercentSpinBox_->setSuffix("%");
                mutationPercentSpinBox_->setValue(15);
                mutationPercentSpinBox_->setSingleStep(5);
            QHBoxLayout* evapRateLayout = new QHBoxLayout;
                evapRateLayout->addWidget(evapRateLabel);
                evapRateLayout->addWidget(evaporationRate_);

            // Minimum pheromone
            QLabel* minPheroLabel = new QLabel(tr("Minimum pheromone"), this);
            minimumPheromone_ = new QDoubleSpinBox(this);
                minimumPheromone_->setRange(1, 100000);
                minimumPheromone_->setValue(1.0);
            QHBoxLayout* minPheroLayout = new QHBoxLayout;
                minPheroLayout->addWidget(minPheroLabel);
                minPheroLayout->addWidget(minimumPheromone_);

            // Maximum pheromone
            QLabel* maxPheroLabel = new QLabel(tr("Maximum pheromone"), this);
            maximumPheromone_ = new QDoubleSpinBox(this);
                maximumPheromone_->setRange(1, 100000);
                maximumPheromone_->setValue(200.);
            QHBoxLayout* maxPheroLayout = new QHBoxLayout;
                maxPheroLayout->addWidget(maxPheroLabel);
                maxPheroLayout->addWidget(maximumPheromone_);

            startACAButton_ = new QPushButton(tr("&Start"), this);

            QVBoxLayout* acaMainLayout = new QVBoxLayout;
                acaMainLayout->addLayout(antCountLayout);
                acaMainLayout->addLayout(evapRateLayout);
                acaMainLayout->addLayout(minPheroLayout);
                acaMainLayout->addLayout(maxPheroLayout);
                acaMainLayout->addItem(new QSpacerItem(0,10, QSizePolicy::Expanding, QSizePolicy::Expanding));
                acaMainLayout->addWidget(startACAButton_);


            QWidget* acWidget = new QWidget(this);
                acWidget->setLayout(acaMainLayout);
                acWidget->setEnabled(true);

        algorithmsTab_ = new QTabWidget(this);
        algorithmsTab_->addTab(accurWidget, tr("精确求解"));
            //algorithmsTab_->addTab(gaWidget, tr("Genetic algorithm"));
        algorithmsTab_->addTab(approxWidget, tr("近似求解"));
            //algorithmsTab_->addTab(acWidget, tr("Ant colony algorithm"));

        // Graph数据选择.
           combo =new QComboBox;
           combo->addItem("生成随机数据");
           combo->addItem("default.tsp");
           combo->addItem("ulysses16.tsp");
           combo->addItem("dj38.tsp");
           combo->addItem("xqf131.tsp");
           combo->addItem("pbl395.tsp");
           connect(combo,SIGNAL(currentIndexChanged(int)),this,SLOT(comboChangeValidity(int)));
            QLabel* graphCountLabel = new QLabel(tr("城市数量："), this);
            graphItemCount_ = new QSpinBox(this);
                graphItemCount_->setRange(3, 1000);
                graphItemCount_->setValue(10);
            QHBoxLayout* graphCountLayout = new QHBoxLayout;
                graphCountLayout->addWidget(graphCountLabel);
                graphCountLayout->addWidget(graphItemCount_);

            QLabel* graphSizeXLabel = new QLabel(tr("地图高度："), this);
            graphSizeX_ = new QSpinBox(this);
                graphSizeX_->setRange(100, 10000);
                graphSizeX_->setValue(800);
            QHBoxLayout* graphSizeXLayout = new QHBoxLayout;
                graphSizeXLayout->addWidget(graphSizeXLabel);
                graphSizeXLayout->addWidget(graphSizeX_);

            QLabel* graphSizeYLabel = new QLabel(tr("地图宽度："), this);
            graphSizeY_ = new QSpinBox(this);
                graphSizeY_->setRange(100, 10000);
                graphSizeY_->setValue(800);
            QHBoxLayout* graphSizeYLayout = new QHBoxLayout;
                graphSizeYLayout->addWidget(graphSizeYLabel);
                graphSizeYLayout->addWidget(graphSizeY_);

            generateNodes_ = new QPushButton(tr("导入/生成 数据"), this);
            connect(generateNodes_, &QPushButton::clicked, actions_.generateNodes, &QAction::trigger);

            QVBoxLayout* graphLayout = new QVBoxLayout;
                graphLayout->addWidget(combo);
                graphLayout->addLayout(graphCountLayout);
                graphLayout->addLayout(graphSizeXLayout);
                graphLayout->addLayout(graphSizeYLayout);
                graphLayout->addWidget(generateNodes_);
            QGroupBox* graphGroupBox = new QGroupBox(tr("数据选择"), this);
                graphGroupBox->setLayout(graphLayout);

        // Results
         result = new QLabel("算法所求出的最短距离： ",this);
         result->resize(sizeHint());
         ans = new QLabel("实际最短路线长度： ",this);
           timec=new QLabel("算法运行时间： ",this);
        startAllButton_ = new QPushButton(tr("Start &all"), this);

    QVBoxLayout* UILayout = new QVBoxLayout;
        UILayout->addWidget(graphGroupBox);
        UILayout->addStretch(1);
        UILayout->addWidget(result);
        UILayout->addWidget(ans);
        UILayout->addWidget(timec);
        UILayout->addStretch(1);
        UILayout->addWidget(algorithmsTab_);
        UILayout->addItem(new QSpacerItem(0,10, QSizePolicy::Expanding, QSizePolicy::Expanding));


    QWidget* UIWidget = new QWidget(this);
        UIWidget->setLayout(UILayout);

    // General
    QSplitter* centralSpliter = new QSplitter(this);
        centralSpliter->addWidget(graphWidget_);
        centralSpliter->addWidget(UIWidget);


    setCentralWidget(centralSpliter);
}
void MainWindow::comboChangeValidity(int index)
{
    if(index!=0)//load data
    {
        graphItemCount_->setEnabled(false);
        graphSizeX_->setEnabled(false);
        graphSizeY_->setEnabled(false);
    }
    else {
        graphItemCount_->setEnabled(true);
        graphSizeX_->setEnabled(true);
        graphSizeY_->setEnabled(true);
    }
}
