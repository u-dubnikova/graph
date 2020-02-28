#include "convert.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include <qmath.h>
#include <assert.h>

constexpr double BAD_E = -INFINITY;

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
    {
        results.emplace_back(result);
	skipLine(ifs);
    }

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

void splitToHalf(Results & results, unsigned minlen)
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
	for (j=jmax=i;j<results.size()-1 && ( ( jmax-i<=minlen && results[j].cycle != 0) || sign_change(results[j].sigma,results[j+1].sigma)!=hcs);j++)
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

static epsig get_epsigma(const PreparedResult & p0, const PreparedResult &p1, double E, double dEps)
{
    double e1=p1.epsilon, e0=p0.epsilon;
    double s1=p1.sigma, s0=p0.sigma;
    double ds=s1-s0;
    double de=e1-e0;
    double E1=ds/de;
    double ret_sigma=(s0-E1*(e0-dEps))/(1-E1/E);
    double ret_eps=ret_sigma/E+dEps;
    return epsig(ret_eps,ret_sigma);
}

bool get_E0(std::vector<PreparedResult> & results, esig & ret)
{
    double s1,s2,s;    
    size_t i=0,i0;
    if (results.size() == 0)
    {
	ret = esig(BAD_E,BAD_E);
	return false;
    }
    while (results[i].epsilon == 0) 
	i++;
    s1=results[i].sigma*results[i].epsilon;
    s2=results[i].epsilon*results[i].epsilon;
    s=s1/s2;
    i0=i;
    for (i++;i<results.size() && results[i].cycle == 0 && results[i].epsilon == results[i0].epsilon;i++) ;
    if ( i>= results.size() || results[i].cycle != 0 )
    {
	printf("HOLY SHIT\n");
	ret =  esig(BAD_E,BAD_E);
	return false;
    }
    double dmin = 1e-3;
    auto save_i=i;
    while (dmin < 1)
    {
	s1=s2=0;
	s1+=results[save_i].sigma*results[save_i].epsilon;
	s2+=results[save_i].epsilon*results[save_i].epsilon;
	s=s1/s2;
	for (i=save_i+1;i<results.size() && results[i].cycle == 0;i++)
	{
	    s1+=results[i].sigma*results[i].epsilon;
	    s2+=results[i].epsilon*results[i].epsilon;
	    double snew=s1/s2;
	    if (fabs((snew-s)/s) < dmin)
	    {
		ret = esig(snew,results[i].epsilon);
		return true;
	    }
	    s = snew;
	}
	dmin*=10;
    } 

    ret = esig(results[i-1].sigma/results[i-1].epsilon,results[i-1].sigma);
    return false;

}

PreparedResult transform_res(const PreparedResult & x, const PreparedResult & x0, double sgn)
{
    return PreparedResult(0,sgn*(x.sigma-x0.sigma),sgn*(x.epsilon-x0.epsilon));
}

esig get_EN(const std::vector<PreparedResult> & results, size_t & idx)
{
    int cnum = results[idx].cycle;
    std::vector<PreparedResult> v2,v3;
    size_t idx00 = idx;
    if (cnum != 0)
    {
	size_t idx_end;
	for (idx_end=idx;idx_end < results.size() && results[idx_end].cycle == cnum;idx_end++)
	    ;
	int s=2*(cnum%2)-1;
	while (--idx_end >= idx)
	    if (s*results[idx_end].sigma > 0)
		break;
	idx = idx_end+1;
    }
    const PreparedResult & x0 = (cnum == 0 )?PreparedResult():results[idx];
    size_t idx0 = idx;

    double sgn = (cnum%2)?-1:1;
    for (idx++;idx < results.size() && results[idx].cycle == cnum;idx++)
	v2.push_back(transform_res(results[idx],x0,sgn));
    if ( idx0 != 0)
	for (size_t idx2=idx0-1;idx2 > idx00;idx2--)
	    v3.push_back(transform_res(results[idx2],x0,-sgn));
    esig es;
    if (idx0 == 0)
    {
	if (!get_E0(v2,es))
	    es = esig(BAD_E,es.second);
    }
    else if (!get_E0(v3,es))
    {
	esig es2;
	if (!get_E0(v2,es2))
	{
//	    std::cout<<"Failed in both directions, cnum="<<cnum<<"\n";
	    es = esig(BAD_E,es2.second);
	}
	else
	    es = es2;
    }
    if (es.first < 0)
    {
	printf("res=%lf",es.first);
	for (const auto & it:v2)
	    printf(",(%lf;%lf)",it.sigma,it.epsilon);
	printf("\n");
    }
    return es;
}

struct edata 
{
    int cnum;
    double e;
    double sigma;
    edata(int cnum_=0,double e_=BAD_E,double sigma_ = 0):cnum(cnum_),e(e_),sigma(sigma_)
    {
    };
};

template <class URNG> double distort(URNG & re,double mid, double randDelta)
{
#if USE_UNIFORM
    std::uniform_real_distribution<double> unif((1-randDelta)*mid,(1+randDelta)*mid);
    return unif(re);
#elif USE_GAUSS
    std::normal_distribution<double> gauss(mid,randDelta*mid);
    return gauss(re);
#else
    UNUSED(re);
    UNUSED(randDelta);
    return mid;
#endif
}

void FixEE(std::vector<edata> & data)
{
    constexpr double bigDelta = 10;
    constexpr double randDelta = 0.1;
    std::vector<edata> rdata;
    bool modified=false;
    size_t nres = 0;
    double k,b;
    double max_cnum=data[data.size()-1].cnum;

    do
    {
#ifdef DBG_PRINT
	for (const auto & x: data)
	    printf("%d\t%lf\n",x.cnum,x.e);
#endif
        double sumc=0,sumcsq=0,sume=0,sumce=0;
	modified = false;
	double emin = data[0].e, emax=data[0].e;
	for (const auto & x: data)
	{
	    if (x.e == BAD_E)
		continue;
	    if (x.e < emin)
		emin = x.e;
	    if (x.e > emax)
		emax = x.e;
	    double lc=log(x.cnum+max_cnum);
	    nres++;
	    sumc+=lc;
	    sumcsq+=lc*lc;
	    sume+=x.e;
	    sumce+=lc*x.e;
	}
	sumc/=nres;
	sumcsq/=nres;
	sume/=nres;
	sumce/=nres;
	k=(sumce-sumc*sume)/(sumcsq-sumc*sumc);
	b=sume-k*sumc;
#ifdef DBG_PRINT
	printf("emin=%lf, emax=%lf\n",emin, emax);
	printf("k=%lf,b=%lf\n",k,b);
	printf("E[1]=%lf,E[last]=%lf\n",data[1].e,data[data.size()-1].e);
#endif
	for (size_t i = 0;i<data.size();i++)
	{
	    if (data[i].e == BAD_E)
		continue;
	    double calc_e = k*log(data[i].cnum+max_cnum)+b;
	    if (fabs(data[i].e-calc_e)>0.3*fabs(calc_e))
	    {
		modified = true;
		data[i].e = BAD_E;
	    }
	    
	}
    } 
    while (modified);

#ifdef DBG_PRINT
    for (const auto & x: data)
	printf("%d\t%lf\n",x.cnum,x.e);
#endif

    for (auto & x: data)
	if (x.e == BAD_E)
	    x.e = k*log(x.cnum+max_cnum)+b; 
   std::default_random_engine re;
   size_t last = data.size() - 1;
   for (size_t i = 1; i < last ;i++)
	if ( fabs(data[i].e-data[i-1].e) > bigDelta )
	    data[i].e = distort(re,(data[i-1].e+data[i+1].e)/2,randDelta);

   if (fabs(data[last].e-data[last-1].e) > bigDelta)
	data[last].e = distort(re, data[last-2].e+(data[last-1].e-data[last-2].e)/(data[last-1].cnum-data[last-2].cnum)*(data[last].cnum-data[last-2].cnum), randDelta);

}

void saveEE(const std::string & FileName, const std::vector<PreparedResult>& results)
{
    std::ofstream f(FileName+".eee");
    size_t idx=0;
    double eprev=0;
    std::vector<edata> ES;

    do
    {
	int cnum=results[idx].cycle;
	esig es=get_EN(results,idx);
	ES.push_back(edata(cnum,es.first,es.second));
    }
    while (idx<results.size());
    if (ES.size() == 0)
	return;
    FixEE(ES);
    eprev = ES[0].e;
    for (idx=0;idx<ES.size();idx++)
    {
	double EN=ES[idx].e;
	f<<ES[idx].cnum<<"\t"<<EN<<"\t"<<EN-eprev<<std::endl;
	eprev = EN;
    }
}

bool findElas(std::vector<PreparedResult> & results, double dEps,double & E_0,epsig & epsig_1, bool & approx_good)
{
   esig es;
   get_E0(results,es);
   E_0=es.first;
#ifdef DBG_PRINT
   std::cout<<"E_0="<<E_0<<std::endl;
#endif
    if ( E_0 == BAD_E || E_0 < 40. )
	return false;
    size_t i_s;

    for (i_s=1;i_s<results.size() && results[i_s].cycle == 0 && results[i_s].sigma>= E_0*(results[i_s].epsilon-dEps);i_s++) 
	;
    if (i_s>=results.size() || results[i_s].cycle != 0)
	return false;
    epsig_1 = get_epsigma(results[i_s-1],results[i_s],E_0,dEps);
#ifdef DBG_PRINT
    std::cout<<"Sigma_1="<<epsig_1.second<<std::endl;
#endif
    if (epsig_1.second < 0.1)
	return false;
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
    return get_epsigma(results[i_s-1],results[i_s],E,dEps).second;    	
}

void removeLastEven(Results & results)
{
    if (results.size() == 0 || results[results.size()-1].cycle%2 == 1 )
	return;
    int i;
    for (i = results.size()-1;i>=0;i--)
	if (results[i].cycle%2!=0)
	    break;
    results.resize(i+1);
}

void removeDuplicates(Results & results)
{
    Results new_res; 
    new_res.push_back(results[0]);
    size_t k = 0;
    for (size_t n=1;n<results.size();n++)
	if (results[n].epsilon != results[k].epsilon || results[n].sigma != results[k].sigma)
	{
	    k = n;
	    new_res.push_back(results[n]);
	}
    results = new_res;
}


void convertFile(const std::string& inputFileName, const std::string& outputFileName, unsigned minlen)
{
    Results results;
    if (!readResults(inputFileName, results))
        throw std::runtime_error("can't read results");

    cutByCycle(results);
    cutNearZero(results);
    transformResults(results);
    filterResults(results);
    splitToHalf(results,minlen);
    removeLastEven(results);
    removeDuplicates(results);
    if (!printResults(outputFileName, results))
        throw std::runtime_error("can't print results");
}

bool stretchFile(const std::string& inputFileName, const std::string& outputFileName)
{
    std::vector<PreparedResult> results,new_results;
    results.clear();
    std::ifstream ifs(inputFileName);
    if (!ifs)
        return false;

    PreparedResult result;
    while (ifs >> result.cycle >> result.sigma >> result.epsilon)
        results.emplace_back(result);
    size_t idx = 0;
    do
    {
	size_t idx0=idx;
	int cur_cycle = results[idx0].cycle;
	std::vector<double> epss,sigs;
	while (idx < results.size() && results[idx].cycle == cur_cycle)
	{
	    epss.emplace_back(results[idx].epsilon);
	    sigs.emplace_back(results[idx].sigma);
	    idx++;
	}
	if (cur_cycle%2)
	{
	    std::sort(epss.begin(),epss.end(),std::greater<double>());
	    std::sort(sigs.begin(),sigs.end(),std::greater<double>());
	}
	else
	{
	    std::sort(epss.begin(),epss.end());
	    std::sort(sigs.begin(),sigs.end());
	}
	for (size_t i = 0;i < epss.size();i++)
	    new_results.emplace_back(PreparedResult(cur_cycle,sigs[i],epss[i]));
    }
    while (idx<results.size());
    std::ofstream ofs(outputFileName);
    if (!ofs)
        return false;
    ofs<< std::scientific << std::uppercase;
    for (const auto & r: new_results)
	ofs << r.cycle <<"\t" <<r.sigma<<"\t"<<r.epsilon <<std::endl;
    return true;



    return true;

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
    std::ofstream fs(FileName+".sss");
    double sp=0;
    while (it!=src.end())
    {
	int cur_cycle=it->cycle;
	double e0 = it[-1].epsilon;
	double s0 = it[-1].sigma;
	if (dir == 1)
	    dst.push_back(*it++);
	auto it2=it;
	while (it+1!=src.end() && it[1].cycle == cur_cycle)
	    it++;
	while (it != it2 && dir*it->sigma <= dir*(s0+E2*(it->epsilon-e0-dir*delta)))
	    it--;
	it++;
	double sc = fabs(it->sigma-dst[dst.size()-1].sigma);
	fs << it->cycle<<"\t"<<sc<<"\t"<<sc-sp<<"\n";
	sp = sc;
	while (it!=src.end() && it->cycle == cur_cycle)
	    dst.push_back(*it++);
	dir=-dir;
    }
    printResults(FileName+".acv",dst);
}

double intersect(const PreparedResult& r1, const PreparedResult & r2)
{
    double res=(r2.sigma*r1.epsilon-r1.sigma*r2.epsilon)/(r2.sigma-r1.sigma);
    return res;
}
double av(double i, double nstart, double nfinish, double xstart, double xfinish)
{
    double t=(i-nstart)/(nfinish-nstart);
    return (1-t)*xstart+t*xfinish;
}


struct chi1 
{
    int cycle;
    double chi;
};

void makeChi1(const std::vector<PreparedResult> & data, double E, double mul, std::vector<chi1> & res)
{
    res.clear();
    int cnum = 0,save_cycle=-1;
    size_t n=0;
    double save_chi=0;
    double cum_chi=0;;
    while (n<data.size())
    {
	double cur_chi;
	double inter;
	while (n <data.size() && data[n].cycle == cnum)
	    n++;
	if (n >= data.size()-1)
	    break;
	cnum = data[n].cycle;
	while (n+2 < data.size() && sign_change(data[n].sigma,data[n+1].sigma) == 0)
	    n++;
	if (data.size() == n+2)
	    break;
	cur_chi=mul*fabs(inter=intersect(data[n],data[n+1]));
	n++;
	
	if (save_cycle != -1)
	{
	    for (int i=save_cycle+1;i < cnum;i++)
	    {
		cum_chi+=av(i,save_cycle,cnum,save_chi,cur_chi);
		chi1 ch={i-1, cum_chi};
		res.push_back(ch);
	    }
	}
	cum_chi+=cur_chi;
	chi1 ch={cnum-1, cum_chi};
	res.push_back(ch);
	save_cycle=cnum;
	save_chi=cur_chi;
    }
    cum_chi+=mul*fabs(data[data.size()-1].epsilon-data[data.size()-1].sigma/E);
    chi1 ch={cnum,cum_chi};
    res.push_back(ch);
}

chi saveChi(const std::string & FileName, const std::vector<PreparedResult>& orig,const std::vector<PreparedResult> & cut, double E)
{
    const double mul=1;//sqrt(3.)/2.;
    std::vector<chi> Chis;
    std::vector<chi1> res_orig,res_cut;
    makeChi1(orig,E,mul,res_orig);
    makeChi1(cut,E,mul,res_cut);
    size_t norig=0,ncut=0;
    for (;norig<res_orig.size() && ncut <res_cut.size();norig++,ncut++)
    {
	if (res_orig[norig].cycle != res_cut[ncut].cycle )
	{
	    std::cout<<"Cycle mismatch: orig="<<res_orig[norig].cycle<<"res_cut="<<res_cut[ncut].cycle;
	}
	chi ch={res_orig[norig].cycle, res_orig[norig].chi, res_cut[ncut].chi};
	Chis.push_back(ch);
    }
    std::vector<chi> ChiNew(Chis.size());
    std::vector<chi> ChiDelta(Chis.size()-1);
    ChiNew[0]=Chis[0];

    for (size_t i = 1;i< Chis.size(); i++)
    {
	ChiNew[i] = {
	    Chis[i].cycle,
	    Chis[i].chi_orig+Chis[i-1].chi_orig,
	    Chis[i].chi_cut+Chis[i-1].chi_cut
	};
	ChiDelta[i-1]= {
	    ChiNew[i].cycle,
	    ChiNew[i].chi_orig-ChiNew[i-1].chi_orig,
	    ChiNew[i].chi_cut-ChiNew[i-1].chi_cut
	};
    }

    std::ofstream ofs(FileName);
    if (!ofs)
    {
	std::cout<<"Cannot open "<<FileName<<std::endl;
        return ChiNew[ChiNew.size()-1];
    }

    //ofs << std::scientific << std::uppercase;
//    for (const auto& res : Chis)
    for (const auto& res : ChiNew)
        ofs << res.cycle <<'\t' <<res.chi_orig*100 <<'\t' <<res.chi_cut*100 << '\n';

    return ChiNew[ChiNew.size()-1];
}
static inline double sqr(double x)
{
    return x*x;
}

void LQRPT2(const std::vector<RPT2Entry> & data, double & k, double & b)
{
    int n = data.size();
    double cyc_av=0,chi_av=0;
    for (const auto & r: data)
    {
	cyc_av+=log(r.cycle);
	chi_av+=log(r.chi);
    }
    cyc_av/=(double)n;
    chi_av/=(double)n;
    double num=0, denom=0;
    for (const auto & r: data)
    {
	num+=(log(r.cycle)-cyc_av)*(log(r.chi)-chi_av);
	denom+=sqr(log(r.cycle)-cyc_av);
    }
    k=num/denom;
    b=chi_av-k*cyc_av;
}
//const double BEZ11::bincoeff[n]={1,10,45,120,210,252,210,120,45,10,1};
template<> const double BEZ<11>::bincoeff[]={1,11,55,165,330,462,462,330,165,55,11,1};
template<> const double BEZ<21>::bincoeff[]={1,21,210,1330,5985,20349,54264,116280,203490,293930,352716,352716,293930,203490,116280,54264,20349,5985,1330,210,21,1};
template<> const double BEZ<20>::bincoeff[]={1,20,190,1140,4845,15504,38760,77520,125970,167960,184756,167960,125970,77520,38760,15504,4845,1140,190,20,1};
template<> const double BEZ<19>::bincoeff[]={1,19,171,969,3876,11628,27132,50388,75582,92378,92378,75582,50388,27132,11628,3876,969,171,19,1};
