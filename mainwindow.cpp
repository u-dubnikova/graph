#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "convert.h"

#include <QFileDialog>
#include <QMessageBox>
#include <iostream>
#include <sstream>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    plot = ui->plot;
    plot->xAxis->setLabel("Epsilon [-]");
    plot->yAxis->setLabel("Sigma [GPa]");
    plot->xAxis->setRange(-1, 1);
    plot->yAxis->setRange(-1, 1);
    plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);

    connect(ui->actionOpen, SIGNAL(triggered()), this, SLOT(openTriggered()));
    connect(ui->actionConvert, SIGNAL(triggered()), this, SLOT(convertTriggered()));
    connect(ui->actionSave, SIGNAL(triggered()), this, SLOT(saveTriggered()));
    connect(ui->action_temp, SIGNAL(triggered()), this, SLOT(setTemperature()));
    connect(ui->action_report, SIGNAL(triggered()), this, SLOT(saveReport()));
    connect(ui->action_delta, SIGNAL(triggered()), this, SLOT(setDelta()));
    dEpsilon = 0.0002;

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::openTriggered()
{
    auto fileName = QFileDialog::getOpenFileName(this,
                                                 "Open file with results", "",
                                                 "ALV (*.alv);;All Files (*)");
    if (fileName.isEmpty())
           return;

    paintGraph(fileName);
}

bool MainWindow::loadPrepared(std::vector<PreparedResult> & results,const QString &filename)
{
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::information(this, "Error",
            "Unable to open file ");
        return false;
    }
    QTextStream in(&file);
    while (!in.atEnd())
    {
        PreparedResult result;
        if (in.status() != QTextStream::Ok)
        {
            QMessageBox::information(this, "Error",
                "Unable to read file, maybe wrong format");
            return false;
        }
        in >> result;
        if (!in.atEnd())
            results.push_back(result);
    }
    return true;
}

void MainWindow::paintGraph(const QString& filename)
{
    plot->clearPlottables();

    QVector<double> x, y;
    std::vector<PreparedResult> results;
    if (!loadPrepared(results,filename))
    {
        plot->replot();
        return;
    }
    for (PreparedResult & result: results)
    {
        std::cout<<result.epsilon<<" "<<result.sigma<<std::endl;
        x.append(result.epsilon);
        y.append(result.sigma);
    }

    auto graph = new QCPCurve(plot->xAxis, plot->yAxis);
    graph->setData(x, y);

    plot->rescaleAxes();
    plot->replot();
}


void MainWindow::convertTriggered()
{
    QStringList inputFileNames = QFileDialog::getOpenFileNames(this,
                                                        "Open file", "",
                                                        "ALC (*.alc);;All Files (*)");
    if (inputFileNames.isEmpty())
           return;
    if (inputFileNames.size() == 1)
    {
        QFileInfo inputInfo(inputFileNames[0]);
        QString new_name= inputInfo.absolutePath()+"/"+inputInfo.baseName()+".alv";

        auto outputFileName = QFileDialog::getSaveFileName(this,
                                                           "Save results", new_name,
                                                           "ALV (*.alv);;All Files (*)");
        if (outputFileName.isEmpty())
               return;
        QFileInfo info(outputFileName);
        if (info.suffix().isEmpty())
            outputFileName += ".alv";
        try
        {
            convertFile(inputFileNames[0].toStdString(), outputFileName.toStdString());
        }
        catch (std::runtime_error& e)
        {
            QMessageBox::information(this, "Error",
                QString("Unable to convert files: ") + e.what());
            return;
        }
        paintGraph(outputFileName);
    }
    else
    {
           for (auto & inputFileName: inputFileNames)
           {
               QFileInfo inputInfo(inputFileName);
               QString outputFileName = inputInfo.absolutePath()+"/"+inputInfo.baseName()+".alv";
               try
               {
                   convertFile(inputFileName.toStdString(), outputFileName.toStdString());
               }
               catch (std::runtime_error& e)
               {
                   QMessageBox::information(this, "Error",
                       QString("Unable to convert files: ") + e.what());
               }
           }

    }
}

void MainWindow::saveTriggered()
{
    auto fileName = QFileDialog::getSaveFileName(this,
                                                 "Save image", "",
                                                 "PNG (*.png);;All Files (*)");
    if (fileName.isEmpty())
           return;

    if (!plot->savePng(fileName)) {
        QMessageBox::information(this, "Error",
            "Unable to save graph image");
    }
}

void MainWindow::setTemperature() {
    bool bOk;
    double temp = QInputDialog::getDouble( this,
                                     "Введите температуру",
                                     "Температура:",
                                     0,
                                     -30,
                                     5000,
                                     0,
                                     &bOk
                                    );
    if (!bOk) {
        // Была нажата кнопка Cancel
        return;
    }
    std::cout<<"Set Temperature:"<<temp<<std::endl;
    QFileDialog dialog(this);
    QStringList fileNames = QFileDialog::getOpenFileNames(this,
                                                                           "Open file", "",
                                                                           "ALV (*.alv);;All Files (*)");
    for (QString & fname: fileNames)
    {
	ReportEntry re;
        if (!loadPrepared(re.results,fname))
            continue;
	re.FileName = fname;
	re.valid=findElas(re.results,dEpsilon,re.E,re.S,re.approx_good);
	report.emplace(temp,re);
    }
    temps.emplace(temp);
}
void MainWindow::setDelta()
{
    bool bOk;
    double delta = QInputDialog::getDouble( this,
                                     "Введите дельту",
                                     "Дельта:",
                                     dEpsilon*100,
                                     0,
                                     0.03,
                                     3,
                                     &bOk
                                    );
    if (!bOk) {
        // Была нажата кнопка Cancel
        return;
    }
    dEpsilon = delta/100;
}

void MainWindow::saveReport() {
    auto fileName = QFileDialog::getSaveFileName(this,
                                                 "Save Report", "",
                                                 "RPT (*.rpt);;All Files (*)");
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly)) {
        QMessageBox::information(this, "Error",
            "Unable to open file ");
        return;
    }
    QTextStream out(&file);
    out<<"Delta:"<<dEpsilon<<"\n";
    for (double temp: temps)
    {
        out<<"Temperature: "<<temp<<"\n";
        auto range = report.equal_range(temp);
	int count=0;
	double sE=0.,sS=0.,s2=0.,E2;
        for (auto it=range.first;it!=range.second;it++)
        {
	    ReportEntry & re=it->second;
	    out<<re.FileName<<"::";
	    if (!re.valid)
	    {
		out<<"not found\n";
		continue;
	    }
	    sE+=re.E;
	    sS+=re.S;
	    count++;
	    out<<re.E<<"\t"<<re.S<<"\n";
	    if (re.approx_good)
		    out<<"Approximation is good\n";
	    else
		    out<<"Approximation is not sufficient, try to increase delta\n";
        }
	sE/=count;
	sS/=count;
	if (!count)
	    continue;
	out<<"Average:: "<<sE<<"\t"<<sS<<"\n";
	count = 0;
        for (auto it=range.first;it!=range.second;it++)
        {
	    ReportEntry & re=it->second;
	    if (!re.valid)
		continue;
	    s2+=findSigma2(re.results, dEpsilon, sE);
	    count++;
	}
	s2/=count;
	out<<"Refined s:"<<s2<<'\n';
	E2 = 1./(1./sE+dEpsilon/s2);
	out<<"Refined E:"<<E2<<'\n';
    }
}