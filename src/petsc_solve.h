#include<iostream>
#include<cmath>

using std::cout;
using std::endl;

PetscErrorCode pet_solve(Mat &A, Vec &src, Vec &exact_sol,int N,std::vector<double> &c_sol, std::vector<double> &masa_sources,std::vector<double> &masa_solution);
