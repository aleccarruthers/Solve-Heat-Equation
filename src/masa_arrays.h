#include<iostream>
#include<masa.h>
#include<vector>

std::vector<double> masa_sources(int N, int dimension,int order,double length,double frequency,double K,std::string loglevel);
std::vector<double> masa_solution(int N, int dimension,double length,double frequency,double K,std::string loglevel);
