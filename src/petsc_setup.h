#include<iostream>

using std::cout;
using std::vector;

PetscErrorCode petsc_set(int N,int dimension, int finorder, Mat &A_p, std::string mode);
