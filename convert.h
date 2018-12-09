#ifndef CONVERT_H
#define CONVERT_H

#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include"result.h"

void SaveCutResults(
    const std::string & FileName,
    const std::vector<PreparedResult> & src,
    double E2, double delta);

void convertFile(const std::string& inputFileName, const std::string& outputFileName);
bool findElas(std::vector<PreparedResult> & results, double dEps,double & E_0,double & sig_1,bool & approx_good);
double findSigma2(std::vector<PreparedResult> & results, double dEps, double E);

class BEZ11
{
    static constexpr int n = 12;
    static const double bincoeff[n];
public:   
    struct point
    {
	double x,y;
    } points[n];
    BEZ11(std::vector<PreparedResult> & avRes)
    {
	for (int i=0;i<n;i++)
	{
	    points[i].x=avRes[i].epsilon;
	    points[i].y=avRes[i].sigma;
	}
    }
    point operator()(double t) const
    {
	point res={0,0};
	double tt=pow(1-t,n-1);
	if (t < 0.00001)
	    return points[0];
	if (t > 1-0.00001)
	    return points[n-1];
	for (int i=0;i<=n-1;i++)
	{
	    res.x+=tt*bincoeff[i]*points[i].x;
	    res.y+=tt*bincoeff[i]*points[i].y;
	    tt*=t/(1-t);
	}
	return res;
    }
};


#endif // CONVERT_H
