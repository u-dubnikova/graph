#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <map>
#include <set>
#include "qcustomplot.h"

#include "result.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void openTriggered();
    void convertTriggered();
    void saveTriggered();
    void setTemperature();
    void saveReport();
    void setDelta();

private:
    void paintGraph(std::vector<PreparedResult> & data, const QString& name, QColor color);
    void paintGraph(const QString& filename);
    bool loadPrepared(std::vector<PreparedResult> & results,const QString &filename);

    Ui::MainWindow *ui;
    QCustomPlot* plot;
    double dEpsilon;
    struct ReportEntry
    {
	QString FileName;
	bool valid;
	bool approx_good;
	double E,S;
	std::vector<PreparedResult> results;
    };
    std::multimap<double,ReportEntry> report;
    std::set<double> temps;
};

#endif // MAINWINDOW_H
