#ifndef CONVERT_H
#define CONVERT_H

#include <stdexcept>
#include <string>
#include <vector>
#include"result.h"

void convertFile(const std::string& inputFileName, const std::string& outputFileName);
bool findElas(std::vector<PreparedResult> & results, double dEps,double & E_0,double & sig_1,bool & approx_good);

#endif // CONVERT_H
