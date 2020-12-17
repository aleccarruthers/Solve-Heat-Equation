#include<iostream>
#include<vector>
#include<iomanip>

using std::cout;
using std::endl;
using std::vector;

void hdf5_prep(std::vector<double> masa_solution, std::vector<double> numeric_sol, double *hd_masa, double *hd_sol,double error_est,double *hd_l2, int dimension, double *coord_x, double *coord_y, double length);
