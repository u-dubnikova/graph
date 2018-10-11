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
	    results[k].cycle=hcn;
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
//    std::cout<<"SIGMA_0="<<sig_0<<",EPSILON_0="<<eps_0<<std::endl;
//    std::cout<<"E_0="<<E_0<<std::endl;
    size_t i_s;

    for (i_s=1;i_s<results.size() && results[i_s].cycle == 0 && results[i_s].sigma>= E_0*(results[i_s].epsilon-dEps);i_s++) 
	;
    if (i_s>=results.size() || results[i_s].cycle != 0)
	return false;
    double e1=results[i_s].epsilon,e0=results[i_s-1].epsilon;
    s1=results[i_s].sigma;
    double s0=results[i_s-1].epsilon;
    double ds=s1-s0;
    double em=(e0*s1-e1*s0)/ds-dEps;
    sig_1=E_0*em/(1-E_0*(e1-e0)/ds);
//    std::cout<<"Sigma_1="<<sig_1<<std::endl;
//
    if (E_0 > bad_E0 || s1 < bad_s1)
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
