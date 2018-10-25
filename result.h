#ifndef RESULT_H
#define RESULT_H

#include <QTextStream>

struct PreparedResult
{
    int cycle = 0;
    double sigma = 0.0;
    double epsilon = 0.0;
    PreparedResult(int cycle_=0,double sigma_=0, double epsilon_=0):
	cycle(cycle_),sigma(sigma_),epsilon(epsilon_) 
    {
    }

    friend QTextStream& operator>>(QTextStream& is, PreparedResult& result)
    {
        return is >> result.cycle >> result.sigma >> result.epsilon;
    }

    friend QTextStream& operator<<(QTextStream& os, const PreparedResult& result)
    {
        return os << result.cycle << '\t' << result.sigma << '\t' << result.epsilon;
    }
};

#endif // RESULT_H
