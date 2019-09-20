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
    connect(ui->action_mincycle, SIGNAL(triggered()), this, SLOT(setMinlen()));
    connect(ui->action_clear, SIGNAL(triggered()), this, SLOT(doClear()));
    dEpsilon = 0.0002;
    minlen = 10;
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::openTriggered()
{
    auto fileName = QFileDialog::getOpenFileName(this,
                                                 "Open file with results", "",
                                                 "ALV (*.alv);;Cut ALV(*.acv);;Last cycle(*.acc);;Chi (*.chi);;Reports (*.rpt);;Chi reports (*.rpt2);;All Files (*)");
    if (fileName.isEmpty())
           return;

    paintGraph(fileName);
}

bool MainWindow::loadPrepared(std::vector<PreparedResult> & results,const QString &filename)
{
    results.clear();
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

bool MainWindow::loadReport(std::vector<SavedReport> & rep, const QString & filename)
{
    rep.clear();
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::information(this, "Error",
            "Unable to open file ");
        return false;
    }
    QTextStream in(&file);
    enum {NONE,TEMP,CURVE} state = NONE;
    SavedReport r;
    while (!in.atEnd())
    {
	QString str = in.readLine();
	switch (state)
	{
	    case NONE:
		if (!str.contains("Temperature:"))
		    break;
		if (sscanf(str.mid(str.indexOf(":")+1).toStdString().c_str(),"%lf",&r.temp) != 1)
		    return false;
		state = TEMP;
		break;
	    case TEMP:
		if (str.contains("Temperature:"))
		{
		    if (sscanf(str.mid(str.indexOf(":")+1).toStdString().c_str(),"%lf",&r.temp) != 1)
			return false;
		    break;
		}
		if (str == "Average curve:")
		{
		    state =  CURVE;
		    r.avCurve.push_back(PreparedResult(0,0,0));
		}
		break;
	    case CURVE:
		PreparedResult pr;
		if (sscanf(str.toStdString().c_str(),"%lf%lf",&pr.sigma,&pr.epsilon) == 2)
		{
		    pr.cycle=0;
		    r.avCurve.push_back(pr);
		    break;
		}
		if (r.avCurve.size() != 0)
		{
		    rep.push_back(r);
		    r.avCurve.clear();
		}
		if (str.contains("Temperature:"))
		{
		    if (sscanf(str.mid(str.indexOf(":")+1).toStdString().c_str(),"%lf",&r.temp) != 1)
			return false;
		    state = TEMP;
		    break;
		}
		state = NONE;
	}
    }
    if (in.atEnd() && state == CURVE && r.avCurve.size() != 0)
	rep.push_back(r);

    return true;
}

void MainWindow::paintGraph(std::vector<PreparedResult> & data, const QString& name, QColor color)
{
    plot->xAxis->setLabel("Epsilon [-]");
    plot->yAxis->setLabel("Sigma [GPa]");
    plot->xAxis->setScaleType(QCPAxis::stLinear);
    plot->yAxis->setScaleType(QCPAxis::stLinear);
    QVector<double> x, y;
    double epsmax=data[0].epsilon,epsmin=epsmax;
    double sigmax=data[0].sigma,sigmin=sigmax;
    for (PreparedResult & result: data)
    {
	if (result.epsilon > epsmax )
	    epsmax = result.epsilon;
	else if (result.epsilon < epsmin )
	    epsmin = result.epsilon;
	if (result.sigma > sigmax )
	    sigmax = result.sigma;
	else if (result.sigma < sigmin )
	    sigmin = result.sigma;
        x.append(result.epsilon);
        y.append(result.sigma);
    }
    QSharedPointer<QCPAxisTickerFixed> xtick(new QCPAxisTickerFixed());
    QSharedPointer<QCPAxisTickerFixed> ytick(new QCPAxisTickerFixed());
    xtick->setTickStep(exp(floor(log10(epsmax-epsmin)-0.5)*log(10)));
    ytick->setTickStep(exp(floor(log10(sigmax-sigmin)-0.5)*log(10)));
    plot->xAxis->setTicker(xtick);
    plot->yAxis->setTicker(ytick);
    plot->xAxis->setRange(epsmin,epsmax);
    plot->yAxis->setRange(sigmin, sigmax);

    auto graph = new QCPCurve(plot->xAxis, plot->yAxis);

    graph->setData(x, y);
    graph->setName(name);
    graph->setPen(QPen(color));
}

bool MainWindow::loadChi(const QString & filename, QVector<double> & x, QVector<double> & yorig, QVector<double> & ycut)
{
    x.clear();
    yorig.clear();
    ycut.clear();
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::information(this, "Error",
            "Unable to open file ");
        return false;
    }
    QTextStream in(&file);
    while (!in.atEnd())
    {
	double ncyc,chiorig,chicut;
        if (in.status() != QTextStream::Ok)
        {
            QMessageBox::information(this, "Error",
                "Unable to read file, maybe wrong format");
            return false;
        }
        in >> ncyc >>chiorig>>chicut;
        if (in.atEnd())
	    break;
	x.push_back(ncyc+1);
	yorig.push_back(chiorig);
	ycut.push_back(chicut);
    }
    return true;
}

void MainWindow::paintChi(const QString & filename)
{
    QVector<double> yorig,ycut,x; 
    if (!loadChi(filename,x,yorig,ycut))
	return;
    auto graph = new QCPCurve(plot->xAxis, plot->yAxis);
    graph->setData(x, yorig);
    graph->setName("Original");
    graph->setPen(QPen(Qt::blue));
    auto graph2 = new QCPCurve(plot->xAxis, plot->yAxis);
    graph2->setData(x, ycut);
    graph2->setName("Cut");
    graph2->setPen(QPen(Qt::red));
    plot->xAxis->setLabel("n [-]");
    plot->yAxis->setLabel("Chi [%]");
    plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    QSharedPointer<QCPAxisTickerLog> xtick(new QCPAxisTickerLog);
    xtick->setLogBase(1.05);
    plot->xAxis->setNumberPrecision(0);
    plot->xAxis->setTicker(xtick);
    plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    QSharedPointer<QCPAxisTickerLog> ytick(new QCPAxisTickerLog);
    ytick->setLogBase(1.25892541179416721042);
    plot->yAxis->setNumberPrecision(2);
    plot->yAxis->setTicker(ytick);
    plot->xAxis->setRange(0, x[x.size()-1]);
    plot->yAxis->setRange(0, yorig[yorig.size()-1]);

    plot->rescaleAxes();
    plot->replot();
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


void MainWindow::paintRPT(const QString & filename)
{
    std::vector<SavedReport> rep;
    if (!loadReport(rep,filename))
	return;
    plot->clearPlottables();
    for (size_t ntemp=0;ntemp<rep.size();ntemp++)
    {
	std::vector<PreparedResult> paintRes;
	paintRes.push_back(PreparedResult(0,0,0));
	BEZ<20> b11(rep[ntemp].avCurve);

	for (int i=0;i<=100;i++)
	{
	    BEZ<20>::point p=b11((double)i/100);
	    paintRes.push_back(PreparedResult(0,p.y,p.x));
	}
	paintGraph(paintRes,"Температура: "+QString::number(rep[ntemp].temp), getMyColor(ntemp)); 
    }
    plot->rescaleAxes();
    plot->replot();
}

void MainWindow::paintRPT2(const QString & filename)
{
    QFile file(filename);
    std::vector<RPT2Entry> entries;
    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::information(this, "Error",
            "Unable to open file ");
        return;
    }
    QTextStream in(&file);
    QString str;
    double temp=-1;
    int pos;
    while (in.readLineInto(&str))
    {
	if (str.indexOf("Temperature:")==0)
	{
	    if (sscanf(str.toStdString().c_str(), "Temperature: %lf", &temp)!=1)
	    {
		QMessageBox::information(this, "Error",
		    "invalid file format");
		return;
	    }
	}
	else if ( ( pos = str.indexOf("::") ) != -1)
	{
	    double chi, eps;
	    int ncyc;
	    QString chiname = str.mid(0,pos-3)+"chi";
	    if (sscanf(str.mid(pos+2).toStdString().c_str(),"%d%lf%lf",&ncyc,&eps,&chi)!=3)
	    {
		QMessageBox::information(this, "Error",
		    "invalid file format");
		return;
	    }
	    RPT2Entry rp;
	    rp.temp = temp;
	    rp.cycle=ncyc+1;
	    rp.chi = chi;
	    if (!loadChi(chiname,rp.x,rp.yorig,rp.ycut))
		return;
	    entries.push_back(rp);
	}
	else
	{
	    QMessageBox::information(this, "Error",
		"invalid file format");
	    return;

	}
    }
    int cur_num=0;
    int max_cycle=0;
    double max_chi=0;
    double cur_temp=entries[0].temp;
    QVector<double> x,y;
    char buffer[256];
    
    for (const auto &r: entries)
    {
	if (r.temp !=cur_temp)
	{
	    QCPGraph * gr = new QCPGraph(plot->xAxis, plot->yAxis);
	    gr->setAdaptiveSampling(false);
	    gr->setLineStyle(QCPGraph::lsNone);
	    gr->setScatterStyle(QCPScatterStyle::ssCircle);
	    gr->setPen(QPen(getMyColor(cur_num)));
	    gr->addData(x,y);
	    snprintf(buffer,256,"%G",cur_temp);
	    gr->setName(buffer);
	    cur_num++;
	    cur_temp=r.temp;
	    x.clear();
	    y.clear();
	}
	x.push_back(r.cycle);
	y.push_back(r.chi);
	if (r.cycle > max_cycle)
	    max_cycle = r.cycle;
	if (r.chi > max_chi )
	    max_chi = r.chi;
	auto graph2 = new QCPCurve(plot->xAxis, plot->yAxis);
	graph2->setData(r.x, r.ycut);
	graph2->removeFromLegend();
	graph2->setPen(QPen(getMyColor(cur_num)));
    }
    QCPGraph * gr = new QCPGraph(plot->xAxis, plot->yAxis);
    gr->setAdaptiveSampling(false);
    gr->setLineStyle(QCPGraph::lsNone);
    gr->setScatterStyle(QCPScatterStyle::ssCircle);
    gr->setPen(QPen(getMyColor(cur_num)));
    snprintf(buffer,256,"%G",cur_temp);
    gr->setName(buffer);
    gr->addData(x,y);
    double k,b;
    LQRPT2(entries,k,b);
    x.clear();
    y.clear();
    double dlogx=log(max_cycle)/20;
    for (int i=0;i<=20;i++)
    {
	x.push_back(exp(i*dlogx));
	y.push_back(exp(k*i*dlogx+b));
    }
    auto graph = new QCPCurve(plot->xAxis, plot->yAxis);
    graph->setData(x, y);
    snprintf(buffer,256,"reg: lg y = 1/%lf lg x+ lg(%lf)",1./k,exp(b));
    //snprintf(buffer,256,"reg: y = %lf*x^%lf",exp(b),k);

    graph->setName(buffer);
    graph->setPen(QPen(QBrush(Qt::black),1,Qt::DashLine));
    plot->xAxis->setLabel("n [-]");
    plot->yAxis->setLabel("Chi [%]");
    plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    QSharedPointer<QCPAxisTickerLog> xtick(new QCPAxisTickerLog);
    xtick->setLogBase(1.05);
    plot->xAxis->setNumberPrecision(0);
    plot->xAxis->setTicker(xtick);
    plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    QSharedPointer<QCPAxisTickerLog> ytick(new QCPAxisTickerLog);
    ytick->setLogBase(1.25892541179416721042);
    plot->yAxis->setNumberPrecision(2);
    plot->yAxis->setTicker(ytick);
    plot->xAxis->setRange(0, max_cycle);
    plot->yAxis->setRange(0, max_chi);

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

    if (filename.mid(filename.size()-5) == ".rpt2")
    {
	paintRPT2(filename);
	return;
    }

    if (filename.mid(filename.size()-4) == ".rpt")
    {
	paintRPT(filename);
	return;
    }


    std::vector<PreparedResult> results;
    QString origfilename;
    std::vector<PreparedResult> origresults;
    if (!loadPrepared(results,filename))
    {
        plot->replot();
        return;
    }
    QColor c=Qt::blue;
    bool PaintOrig = false;
    if (filename.mid(filename.size()-4) == ".acv")
    {

	if (QMessageBox::question(this,"Строить ли оригинал","Нужно ли строить исходную кривую?",QMessageBox::Yes|QMessageBox::No,QMessageBox::No) == QMessageBox::Yes)
	{
	    PaintOrig = true; 
	    origfilename = filename.mid(0,filename.size()-4)+".alv";
	    if (!loadPrepared(origresults,origfilename))
	    {
		PaintOrig = false;
	    }
	}
	c = Qt::red;
    }
    else if (filename.mid(filename.size()-4) == ".acc")
	c = Qt::green;
    bool bOk;
    int ini_cycle=results[0].cycle;
    int fin_cycle=results[results.size()-1].cycle;
    int n0=ini_cycle;
    int n1=fin_cycle;
    if (ini_cycle + 1 < fin_cycle)
    {

       n0 = QInputDialog::getInt( this,
                                     "Номер начального полуцикла",
                                     QString::asprintf("Начальный полуцикл (%d-%d)",ini_cycle,fin_cycle),
                                     ini_cycle,
                                     ini_cycle,
                                     fin_cycle,
                                     1,
                                     &bOk
                                    );
	if (!bOk) {
	    n0=ini_cycle;
	}
	n1 = QInputDialog::getInt( this,
					 "Номер конечного полуцикла",
					 "Конечный полуцикл:",
					 fin_cycle,
					 n0,
					 fin_cycle,
					 1,
					 &bOk
					);
	if (!bOk) {
	    n1=fin_cycle;
	}
    }


    size_t c0=0;
    while (c0 < results.size() && results[c0].cycle < n0 )
	c0++;
    size_t c1=c0;
    while (c1 < results.size() && results[c1].cycle <= n1 )
	c1++;
    results.erase(results.begin()+c1,results.end());
    results.erase(results.begin(),results.begin()+c0);
    if (PaintOrig)
    {
	while (c0 < origresults.size() && origresults[c0].cycle < n0 )
	    c0++;
	c1=c0;
	while (c1 < origresults.size() && origresults[c1].cycle <= n1 )
	    c1++;
	origresults.erase(origresults.begin()+c1,origresults.end());
	origresults.erase(origresults.begin(),origresults.begin()+c0);
	paintGraph(origresults,origfilename, Qt::blue);
    }

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
            convertFile(inputFileNames[0].toStdString(), outputFileName.toStdString(), minlen);
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
                   convertFile(inputFileName.toStdString(), outputFileName.toStdString(),minlen);
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
	epsig ep;
	if ((re.valid=findElas(re.results,dEpsilon,re.E,ep,re.approx_good)))
	{
	    re.eps = ep.first;
	    re.S = ep.second;
	}
	report.emplace(temp,re);

	if (re.valid)
	    saveEE(fname.mid(0,fname.size()-4).toStdString(),re.results);
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
                                     100,
                                     3,
                                     &bOk
                                    );
    if (!bOk) {
        // Была нажата кнопка Cancel
        return;
    }
    dEpsilon = delta/100;
}

void MainWindow::setMinlen()
{
    bool bOk;
    int newlen = QInputDialog::getInt( this,
                                     "Введите минимальную длину полуцикла",
                                     "Минимальная длина полуцикла:",
                                     minlen,
                                     0,
                                     100000,
                                     1,
                                     &bOk
                                    );
    if (!bOk || minlen < 0) {
        // Была нажата кнопка Cancel
        return;
    }
    minlen = newlen;
}

void MainWindow::saveReport() {
    QString fileName = QFileDialog::getSaveFileName(this,
                                                 "Save Report", "",
                                                 "RPT (*.rpt);;All Files (*)");
    QString fileName2 = fileName+"2";
    QFile file(fileName),file2(fileName2);


    if (!file.open(QIODevice::WriteOnly) || !file2.open(QIODevice::WriteOnly)) {
        QMessageBox::information(this, "Error",
            "Unable to open file ");
        return;
    }

    QTextStream out(&file);
    QTextStream out2(&file2);
    out<<"Delta:"<<dEpsilon<<"\n";
    int ntemp=0;
    plot->clearPlottables();
    for (double temp: temps)
    {
        out<<"Temperature: "<<temp<<"\n";
        out2<<"Temperature: "<<temp<<"\n";
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
	    QString out_name=re.FileName.mid(0,re.FileName.size()-4);
	    SaveCutResults(out_name.toStdString(),re.results,re.E,dEpsilon);
	    std::vector<PreparedResult> cutRes;
	    loadPrepared(cutRes,out_name+".acv");
	    QString chi_name=out_name+".chi";
	    out2<<re.FileName<<"::";
	    chi ch=saveChi(chi_name.toStdString(),re.results,cutRes,re.E);
	    out2<<ch.cycle<<'\t'<<re.E<<'\t'<<ch.chi_cut*100<<'\n';
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
	std::vector<PreparedResult> paintRes;
	paintRes.push_back(PreparedResult(0,0,0));
	BEZ<20> b11(avResults);

	for (int i=0;i<=100;i++)
	{
	    BEZ<20>::point p=b11((double)i/100);
	    paintRes.push_back(PreparedResult(0,p.y,p.x));
	}
	paintGraph(paintRes,"Температура: "+QString::number(temp),getMyColor   (ntemp)); 
	ntemp++;
    }
    plot->rescaleAxes();
    plot->replot();
}
