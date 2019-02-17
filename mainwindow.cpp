#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "convert.h"

#include <QFileDialog>
#include <QMessageBox>
#include <iostream>
#include <sstream>
#include <limits>

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
    plot->setAutoAddPlottableToLegend(true);
    plot->legend->setVisible(true);

    connect(ui->actionOpen, SIGNAL(triggered()), this, SLOT(openTriggered()));
    connect(ui->actionConvert, SIGNAL(triggered()), this, SLOT(convertTriggered()));
    connect(ui->actionSave, SIGNAL(triggered()), this, SLOT(saveTriggered()));
    connect(ui->action_temp, SIGNAL(triggered()), this, SLOT(setTemperature()));
    connect(ui->action_report, SIGNAL(triggered()), this, SLOT(saveReport()));
    connect(ui->action_delta, SIGNAL(triggered()), this, SLOT(setDelta()));
    connect(ui->action_clear, SIGNAL(triggered()), this, SLOT(doClear()));
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
                                                 "ALV (*.alv);;Cut ALV(*.acv);;Chi (*.chi);;All Files (*)");
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

void MainWindow::paintGraph(std::vector<PreparedResult> & data, const QString& name, QColor color)
{
    plot->xAxis->setLabel("Epsilon [-]");
    plot->yAxis->setLabel("Sigma [GPa]");
    plot->xAxis->setRange(-1, 1);
    plot->yAxis->setRange(-1, 1);
    plot->xAxis->setScaleType(QCPAxis::stLinear);
    plot->yAxis->setScaleType(QCPAxis::stLinear);
    QVector<double> x, y;
    for (PreparedResult & result: data)
    {
        x.append(result.epsilon);
        y.append(result.sigma);
    }
    auto graph = new QCPCurve(plot->xAxis, plot->yAxis);

    graph->setData(x, y);
    graph->setName(name);
    graph->setPen(QPen(color));
}
void MainWindow::paintChi(const QString & filename)
{
    QVector<double> yorig,ycut,x; 
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::information(this, "Error",
            "Unable to open file ");
        return;
    }
    QTextStream in(&file);
    while (!in.atEnd())
    {
	double ncyc,chiorig,chicut;
        if (in.status() != QTextStream::Ok)
        {
            QMessageBox::information(this, "Error",
                "Unable to read file, maybe wrong format");
            return;
        }
        in >> ncyc >>chiorig>>chicut;
        if (in.atEnd())
	    break;
	x.push_back(ncyc+1);
	yorig.push_back(chiorig);
	ycut.push_back(chicut);
    }
    auto graph = new QCPCurve(plot->xAxis, plot->yAxis);
    graph->setData(x, yorig);
    graph->setName("Original");
    graph->setPen(QPen(Qt::blue));
    auto graph2 = new QCPCurve(plot->xAxis, plot->yAxis);
    graph2->setData(x, ycut);
    graph2->setName("Cut");
    graph2->setPen(QPen(Qt::red));
    plot->xAxis->setLabel("n [-]");
    plot->yAxis->setLabel("Chi [-]");
    plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    plot->xAxis->setRange(0, x[x.size()-1]);
    plot->yAxis->setRange(0, yorig[yorig.size()-1]);

    plot->rescaleAxes();
    plot->replot();
}
void MainWindow::paintGraph(const QString& filename)
{
    plot->clearPlottables();
    if (filename.mid(filename.size()-4) == ".chi")
    {
	paintChi(filename);
	return;
    }

    std::vector<PreparedResult> results;
    if (!loadPrepared(results,filename))
    {
        plot->replot();
        return;
    }
    QColor c=Qt::blue;
    if (filename.mid(filename.size()-4) == ".acv")
	c = Qt::red;
    paintGraph(results,filename, c);
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

void MainWindow::doClear()
{
    dEpsilon = 0.0002;
    temps.clear();
    report.clear();
    plot->clearPlottables();
    plot->replot();
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

static QColor getMyColor(int nres)
{
    switch(nres)
    {
	case 0: return QColor(0,0,0);
	case 1: return QColor(255,0,0);
	case 2: return QColor(0,0,255);
	case 3: return QColor(0,255,0);
	default:
	    return QColor(nres*5,100+nres*2,255-nres*3);
    }
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
    int ntemp=0;
    plot->clearPlottables();
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
	std::vector<std::vector<PreparedResult> > convs;
	double epmax=0;
	for (auto it=range.first;it!=range.second;it++)
	{
	    ReportEntry & re = it->second;
	    if (!re.valid)
		continue;
	    std::vector<PreparedResult> cur;
	    double epcycle = 0;
	    for (const auto & res: re.results)
	    {
		if (res.cycle != 0 )
		    break;
		double epcur = res.epsilon - res.sigma/sE;
		if (epcur > epcycle)
		     epcycle = epcur;
		cur.push_back(PreparedResult(0,res.sigma, epcur));
	    }
	    convs.push_back(cur);
	    if (epcycle > epmax)
		epmax = epcycle;
	    QString out_name=re.FileName.mid(0,re.FileName.size()-4)+".acv";
	    SaveCutResults(out_name.toStdString(),re.results,re.E,dEpsilon);
	    std::vector<PreparedResult> cutRes;
	    loadPrepared(cutRes,out_name);
	    QString chi_name=re.FileName.mid(0,re.FileName.size()-4)+".chi";
	    saveChi(chi_name.toStdString(),re.results,cutRes,re.E);
	}
	std::vector<size_t> positions(convs.size());
	std::vector<PreparedResult> avResults;
	avResults.push_back(PreparedResult(0,0,0));
	for (int i=0;i<=10;i++)
	{
	    double cur_eps = dEpsilon + i*(epmax*0.999-dEpsilon)/10;
	    double av_sig = 0;
	    size_t n = 0;
	    for (size_t j=0;j<convs.size();j++)
	    {
		auto & cv = convs[j];
		size_t k;
		for (k=positions[j];k<cv.size() && cv[k].epsilon < cur_eps;k++) ;
		if (k == cv.size() || cv[k].epsilon < cur_eps)
		    continue;
		av_sig+=cv[k-1].sigma+(cv[k].sigma - cv[k-1].sigma)/(cv[k].epsilon - cv[k-1].epsilon) * (cur_eps - cv[k-1].epsilon);	
		if (k!=0)
		    positions[j]=k-1;
		n++;    
	    }
	    av_sig/=n;
	    if (i == 1)
	    {
		double e0=avResults[1].epsilon,s0=avResults[1].sigma;
		double e1=cur_eps+av_sig/sE, s1=av_sig;
		for (int k=1;k<=10;k++)
		{
		    double e=e0+((e1-e0)*k)/11.;
		    double s=s0+(s1-s0)/(e1-e0)*(e-e0);
		    avResults.push_back(PreparedResult(0,s,e));
		}
	    }
	    avResults.push_back(PreparedResult(0,av_sig, cur_eps+av_sig/sE));
	}
	out<<"Average curve:"<<'\n';
	for (size_t i=1;i<avResults.size();i++)
	    out<<avResults[i].sigma<<'\t'<<avResults[i].epsilon<<'\n';
	std::vector<PreparedResult> paintRes(10*10+1);
	paintRes.push_back(PreparedResult(0,0,0));
	BEZ<21> b11(avResults);

	for (int i=0;i<=100;i++)
	{
	    BEZ<21>::point p=b11((double)i/100);
	    paintRes.push_back(PreparedResult(0,p.y,p.x));
	}
	paintGraph(paintRes,"Температура: "+QString::number(temp),getMyColor   (ntemp)); 
	ntemp++;
    }
    plot->rescaleAxes();
    plot->replot();
}
