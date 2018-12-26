#include "convert.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <qmath.h>

struct Result
{
    int cycle = 0;
    double sigma = 0.0;
    double epsilon = 0.0;

    Result(int cycle_=0, double sigma_=0, double epsilon_=0):
	cycle(cycle_),sigma(sigma_),epsilon(epsilon_)
    {
    }

    Result(const PreparedResult & r):
	cycle(r.cycle),sigma(r.sigma),epsilon(r.epsilon)
    {
    }

    friend std::istream& operator>>(std::istream& is, Result& result)
    {
        char semicolon;
        double cycleTemp, seconds, stroke;
        is >> seconds >> semicolon
           >> cycleTemp >> semicolon
           >> result.sigma >> semicolon
           >> stroke >> semicolon
           >> result.epsilon >> semicolon;

        result.cycle = (int) std::round(cycleTemp);
        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, const Result& result)
    {
        return os << result.cycle << '\t' << result.sigma << '\t' << result.epsilon;
    }
};

using Results = std::vector<Result>;

inline void skipLine(std::istream& is)
{
    std::string s;
    std::getline(is, s);
}

bool readResults(const std::string& filename, Results& results)
{
    results.clear();
    std::ifstream ifs(filename);
    if (!ifs)
        return false;
    ifs >> std::scientific >> std::uppercase;

    Result result;
    skipLine(ifs);
    while (ifs >> result)
        results.emplace_back(result);

    return true;
}

void cutByCycle(Results& results)
{
    bool firstFound = false;
    for (auto rit = results.rbegin(); rit != results.rend(); ++rit)
    {
        if (rit->cycle == 1)
            firstFound = true;

        if (firstFound && rit->cycle != 1)
        {
            results.erase(results.begin(), rit.base());
            break;
        }
    }
}

void cutNearZero(Results& results)
{
    constexpr double linearEpsilon = 0.01;
    for (auto it = results.begin(); it->cycle == 1 && it < results.end() - 1; ++it)
    {
        auto aS = it->sigma;
        auto bS = (it + 1)->sigma;
        auto aE = it->epsilon;
        auto bE = (it + 1)->epsilon;

        if (aS <= 0 || bS <= 0 || aE <= 0 || bE <= 0)
            continue;

        auto deltaS = bS - aS;
        auto deltaE = bE - aE;
	std::cerr<<"DeltaS="<<deltaS<<",DeltaE="<<deltaE<<std::endl;

        if (deltaS < 0 || deltaE < 0)
            continue;

        if (deltaS > linearEpsilon || deltaE > linearEpsilon)
        {
            results.erase(results.begin(), it);
            break;
        }
    }
}

void transformResults(Results& results)
{
    constexpr double sigmaKoef = 4 / (M_PI * 4.37 * 4.37);
    constexpr int l0 = 13;

    for (auto& res : results)
    {
        res.sigma *= sigmaKoef;
        res.epsilon /= l0;
    }
}

void filterResults(Results & results)
{
    int cutIdx = 0;
    int tmpCycle = results.back().cycle;
    int countCycle = -1;
    double maxEps = 0.0, backEps = 0.0;
    std::vector<double> epsilons;
    for (int i = results.size() - 1; i >= 0; --i)
    {
        if (results[i].cycle != tmpCycle)
        {
            if (countCycle != -1)
                epsilons.push_back(maxEps);
            else
            {
                backEps = maxEps;
                cutIdx = results.size() - i - 1;
            }
            maxEps = 0.0;

            if (++countCycle == 5)
                break;
            tmpCycle = results[i].cycle;
        }

        auto absEps = fabs(results[i].epsilon);
        if (absEps > maxEps)
            maxEps = absEps;
    }

    double avgEps = 0.0;
    for (double eps : epsilons)
        avgEps += eps;
    avgEps /= epsilons.size();

    if (qMax(avgEps, backEps) / qMin(avgEps, backEps) >= 2.0 && cutIdx != 0)
        results.erase(results.end() - cutIdx, results.end());

}

int sign_change(double x,double y)
{
    if ( (x == 0 && y < 0 ) || (x > 0 && y == 0) || (x > 0 && y < 0 ))
	return 1;
    if ( (y == 0 && x < 0 ) || (y > 0 && x == 0) || (y > 0 && x < 0 ))
	return -1;
    return 0;	
}

void splitToHalf(Results & results)
{
    int hcn=0,hcs=1;
    std::vector<int> c_nums;
    int cur_cycle = results[0].cycle;
    c_nums.push_back(cur_cycle);
    for (auto & r: results)
    {
	if (r.cycle != cur_cycle)
	{
	    cur_cycle = r.cycle;
	    c_nums.push_back(cur_cycle);
	}
    }
    for (size_t i=0;i<results.size();)
    {
	size_t j,jmax;
	for (j=jmax=i;j<results.size()-1 && sign_change(results[j].sigma,results[j+1].sigma)!=hcs;j++)
	    if (hcs*results[j].epsilon > hcs*results[jmax].epsilon )
		jmax=j;
	if (j == results.size() - 1)	
	{
	    results.erase(results.begin()+i,results.end());
	    return;
	}
	for (size_t k=i;k<=jmax;k++)
	    results[k].cycle=(c_nums[hcn/2]-1)*2+(hcn%2);
	hcn++;
	hcs=-hcs;
	i=jmax+1;
    }
}

bool printResults(const std::string& filename, const Results& results)
{
    std::ofstream ofs(filename);
    if (!ofs)
        return false;

    ofs << std::scientific << std::uppercase;
    for (const auto& res : results)
        ofs << res << '\n';

    return true;
}

constexpr double bad_E0 = 250;
constexpr double bad_s1 = 0.65;

static double get_sigma(const PreparedResult & p0, const PreparedResult &p1, double E, double dEps)
{
    double e1=p1.epsilon, e0=p0.epsilon;
    double s1=p1.sigma, s0=p0.sigma;
    double ds=s1-s0;
    double de=e1-e0;
    double E1=ds/de;
    return (s0-E1*(e0-dEps))/(1-E1/E);
}

bool findElas(std::vector<PreparedResult> & results, double dEps,double & E_0,double & sig_1, bool & approx_good)
{
    double s1=results[0].sigma*results[0].epsilon;
    double s2=results[0].epsilon*results[0].epsilon;
    double s=s1/s2,sp,d,dp,sig_0,eps_0;
    size_t i;
    sp=s;
    approx_good = false;
    for (i=1;i<results.size() && results[i].cycle == 0 && results[i].epsilon == results[0].epsilon;i++) ;
    s1+=results[i].sigma*results[i].epsilon;
    s2+=results[i].epsilon*results[i].epsilon;
    s=s1/s2;
    if ( results[i].cycle != 0 )
	return false;
    d=(s-sp)/(results[i].epsilon-results[0].epsilon);
    
    for (i++;i<results.size() && results[i].cycle == 0;i++)
    {
	double eps=results[i].epsilon, epsp=results[i-1].epsilon;
	sp=s;
	dp=d;
	s1+=results[i].sigma*eps;
	s2+=eps*eps;
	if (eps == epsp)
	    continue;
	s=s1/s2;
	d=(s-sp)/(eps-epsp);
//	std::cout<<"E="<<s<<std::endl;
//	std::cout<<"E'="<<d<<std::endl;
	if (d<0)
	{
	    eps_0=(epsp*d-eps*dp)/(d-dp);
	    sig_0=results[i-1].sigma+(results[i].sigma-results[i-1].sigma)/(eps-epsp)*(eps_0-epsp);
	    E_0=sig_0/eps_0;
	    break;
	}
    }
   std::cout<<"SIGMA_0="<<sig_0<<",EPSILON_0="<<eps_0<<std::endl;
   std::cout<<"E_0="<<E_0<<std::endl;
    size_t i_s;

    for (i_s=1;i_s<results.size() && results[i_s].cycle == 0 && results[i_s].sigma>= E_0*(results[i_s].epsilon-dEps);i_s++) 
	;
    if (i_s>=results.size() || results[i_s].cycle != 0)
	return false;
    sig_1 = get_sigma(results[i_s-1],results[i_s],E_0,dEps);
    std::cout<<"Sigma_1="<<sig_1<<std::endl;
/* 
 * regression: s=ke+b+u
 * k, b - obvious
 * aEps ... averages
*/
    size_t n=0;
    double aSig=0,aEps=0,aSigEps=0,aEps2=0;
    for (n=0;n<results.size() && results[n].cycle == 0;n++)
    {
	double s=results[n].sigma,e=results[n].epsilon;
	aEps+=e;
	aSig+=s;
	aSigEps+=s*e;
	aEps2+=e*e;
    }
    aEps/=n;
    aSig/=n;
    aSigEps/=n;
    aEps2/=n;
    double k=(aSigEps-aSig*aEps)/(aEps2-aEps*aEps);
    double b=aSig-k*aEps,S2=0;
 /* 
  * Find variance estimate S2 for the regression
  * S2=\sum (s_i-ke_i-b)^2/(n-2)
*/
    for (size_t i=0;i<n;i++)
    {
	double s=results[i].sigma,e=results[i].epsilon;
	double as=k*e+b;
	double ds=s-as;
	S2+=ds*ds;
    }
    S2/=(n-2);
    if (S2 < 0.0025*aSig*aSig ) // sqrt(S2)/|aSig| < 5%
	return false;
    while (results[i_s].cycle == 0)
	i_s++;
    double eps2 = results[i_s].epsilon - dEps;
    double sig2 = results[i_s].sigma;
    for (i_s++;i_s<results.size() && results[i_s].cycle == 1 && results[i_s].sigma<= sig2+E_0*(results[i_s].epsilon-eps2);i_s++) 
	;
    if (i_s<=results.size() && results[i_s].cycle == 1 && results[i_s].sigma <= 0 )
	approx_good = true;
    	
    return true;
}

double findSigma2(std::vector<PreparedResult> & results, double dEps, double E)
{
    size_t i_s;
    for (i_s=0;i_s<results.size() && results[i_s].cycle == 0 && results[i_s].sigma>= E*(results[i_s].epsilon-dEps);i_s++) 
	;
//    if (i_s == results.size())
//	std::cout<<"NO INTERSECTION!!!!"<<std::endl;
    return get_sigma(results[i_s-1],results[i_s],E,dEps);    	
}

void convertFile(const std::string& inputFileName, const std::string& outputFileName)
{
    Results results;
    if (!readResults(inputFileName, results))
        throw std::runtime_error("can't read results");

    cutByCycle(results);
    cutNearZero(results);
    transformResults(results);
    filterResults(results);
    splitToHalf(results);
    if (!printResults(outputFileName, results))
        throw std::runtime_error("can't print results");
}

void SaveCutResults(
    const std::string & FileName,
    const std::vector<PreparedResult> & src,
    double E2, double delta)
{
    std::vector<Result> dst;
    dst.clear();
    dst.push_back(PreparedResult(0,0,0));
    auto it = src.begin();
    while (it->sigma > E2*(it->epsilon - delta) && it->cycle == 0 )
    {
	it++;
	if (it == src.end())
	    return;
    }
    while (it->cycle == 0)
    {
	dst.push_back(*it);
	it++;
	if (it == src.end())
	    return;
    }
    int dir=-1;
    while (it!=src.end())
    {
	int cur_cycle=it->cycle;
	double e0 = it[-1].epsilon;
	double s0 = it[-1].sigma;
	while (it!=src.end() && dir*it->sigma > dir*(s0+E2*(it->epsilon-e0-dir*delta)) && it->cycle == cur_cycle)
	    it++;
	while (it!=src.end() && it->cycle == cur_cycle)
	    dst.push_back(*it++);
	dir=-dir;
    }
    printResults(FileName,dst);
}

//const double BEZ11::bincoeff[n]={1,10,45,120,210,252,210,120,45,10,1};
template<> const double BEZ<11>::bincoeff[]={1,11,55,165,330,462,462,330,165,55,11,1};
template<> const double BEZ<21>::bincoeff[]={1,21,210,1330,5985,20349,54264,116280,203490,293930,352716,352716,293930,203490,116280,54264,20349,5985,1330,210,21,1};
template<> const double BEZ<20>::bincoeff[]={1,20,190,1140,4845,15504,38760,77520,125970,167960,184756,167960,125970,77520,38760,15504,4845,1140,190,20,1};
