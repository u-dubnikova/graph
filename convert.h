#ifndef CONVERT_H
#define CONVERT_H

#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <QVector>
#include <utility>
#include"result.h"

#define UNUSED(x) ((void)(x))

struct RPT2Entry
{
    double temp;
    int cycle;
    double chi;
    QVector<double> x,yorig,ycut;
};

struct chi
{
    int cycle;
    double chi_orig;
    double chi_cut;
}; 

void LQRPT2(const std::vector<RPT2Entry> & data, double & k, double & b);

void SaveCutResults(
    const std::string & FileName,
    const std::vector<PreparedResult> & src,
    double E2, double delta);

using epsig = std::pair<double,double>;
using esig = std::pair<double,double>;

void convertFile(const std::string& inputFileName, const std::string& outputFileName, unsigned minlen);
bool findElas(std::vector<PreparedResult> & results, double dEps,double & E_0,epsig & epsig_1,bool & approx_good);
double findSigma2(std::vector<PreparedResult> & results, double dEps, double E);
chi saveChi(const std::string & FileName, const std::vector<PreparedResult>& orig,const std::vector<PreparedResult> & cut, double E);
void saveEE(const std::string & FileName, const std::vector<PreparedResult>& results);
double get_EN(std::vector<PreparedResult> & results, size_t & idx);

template<unsigned N> class BEZ
{
    static constexpr int n = N+1;
    static const double bincoeff[n];
    friend class BEZ<N+1>;
public:   
    struct point
    {
	double x,y;
    } points[n],dpoints[n-1];
    double tch,ych;
    BEZ(std::vector<PreparedResult> & avRes)
    {
	double tch1 = 0;
	for (int i=0;i<n;i++)
	{
	    points[i].x=avRes[i].epsilon;
	    points[i].y=avRes[i].sigma;
	}
	for (int i=0;i<n-1;i++)
	{
	    dpoints[i].x=n*(points[i+1].x-points[i].x);
	    dpoints[i].y=n*(points[i+1].y-points[i].y);
	}
	for (int l=0;l<=20;l++)
	    if (derive(l/20.).y < 0)
	    {
	 	tch1=l/20.;
		break;
	    }
	tch=tch1;	
	for (int l=0;l<20;l++)
	    if (derive(tch1-l/400).y <= 0)
		tch=tch1-l/400.;
	if (tch<0.001)
	    tch=1;
	ych=(*this)(tch).y;
    }
    point operator()(double t) const
    {
	point res={0,0};
	double tt=pow(1-t,n-1);
	if (t < 0.00001)
	    return points[0];
	for (int i=0;i<=n-1;i++)
	{
	    res.x+=tt*bincoeff[i]*points[i].x;
	    res.y+=tt*bincoeff[i]*points[i].y;
	    tt*=t/(1-t);
	}
	if (t>tch)
	    res.y=ych;
	return res;
    }
    point derive(double t) const
    {
	point res={0,0};
	double tt=pow(1-t,n-2);
	if (t < 0.00001)
	    return dpoints[0];
	if (t > 1-0.00001)
	    return dpoints[n-2];
	for (int i=0;i<=n-2;i++)
	{
	    res.x+=tt*BEZ<N-1>::bincoeff[i]*dpoints[i].x;
	    res.y+=tt*BEZ<N-1>::bincoeff[i]*dpoints[i].y;
	    tt*=t/(1-t);
	}
	return res;
    }
};

#endif // CONVERT_H
