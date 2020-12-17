#include<iostream>
#include<iomanip>
#include<math.h>
#include<vector>

using std::vector;
using std::string;
std::vector<double> num_solve(std::string method,std::vector<std::vector<double>> A, std::vector<double> sources,int dimension,int order,double length,double frequency,double K,int max_iterations);
