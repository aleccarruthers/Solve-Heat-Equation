#include<iostream>
#include<vector>

using std::string;
using std::cout;
using std::endl;
using std::vector;

PetscErrorCode petsc_masa_source(int N, int dimension, int finorder, double length, double frequency, double K,std::string loglevel, Vec &masa_src_petsc, Vec &masa_sol_petsc, Vec &num_sol);
