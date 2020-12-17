#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include "petsc.h"
#include "petscksp.h"
using std::cout;
using std::endl;
using std::vector;

PetscErrorCode pet_solve(Mat &A, Vec &src, Vec &exact_sol,int N,std::vector<double> &c_sol, std::vector<double> &masa_sources,std::vector<double> &masa_solution){
	cout << std::setprecision(12);
	PetscErrorCode ierr;
	// Initialize Solver
	KSP ksp;
	// Solution vector
	Vec sol;
	// Pointers for PETSc Vector Convertion
	PetscScalar *sol_test, *src_test, *exact_test;
	PetscScalar tot_entries = masa_sources.size();
	// Needed for converting PETSc vector into C++ Vector
	PetscInt sol_size, src_size, exact_size;
	PetscInt sol_iter, src_iter, exact_iter;
	
	// Absolute Residual tolerance
	PetscReal norm, tol = 1e-12;
	
	// Initialize Solve
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	ierr = VecDuplicate(src,&sol);
	// Solver specifications
	ierr = KSPSetTolerances(ksp,tol,tol,PETSC_DEFAULT,100000);CHKERRQ(ierr);
	// Initialize Convergence variable: positive=convergence, negative=divergence
	KSPConvergedReason reason;
	
	// Solve
	ierr = KSPSolve(ksp,src,sol);CHKERRQ(ierr);
	// Check Convergence
	KSPGetConvergedReason(ksp,&reason);
	if (reason<0){
		printf("Divergence. \n");}
	else{
		PetscInt its;
		KSPGetIterationNumber(ksp,&its);
		printf("Convergence in %d iterations.\n",(int)its);
	}

	// Convert Numeric solution to C++ vector
	//VecGetLocalSize(sol,&sol_size);
	//VecGetArray(sol,&sol_test);
	PetscInt vec_size;
	VecGetSize(sol,&vec_size);
	PetscScalar holder;
	PetscInt c1;
	int c_i = 0;
	for (c1=0;c1<vec_size;c1+=1){
		VecGetValues(sol,1,&c1,&holder);
		double v = (double)holder;
		c_sol.push_back(v);
		c_i+=1;
	}
	//VecRestoreArray(sol,&sol_test);
	// Convert masa_sources to C++ vector
	PetscInt vec_size2;
	VecGetSize(src,&vec_size2);
	PetscScalar temp;
	PetscInt c2;
	int c_i2 = 0;
	for (c2=0;c2<vec_size2;c2+=1){
		VecGetValues(src,1,&c2,&temp);
		double v2 = (double)temp;
		masa_sources.push_back(v2);
		c_i2+=1;
	}
	
	// Convert masa exact solution to C++ vector
	PetscInt vec_size3;
	VecGetSize(exact_sol,&vec_size3);
	PetscScalar temp3;
	PetscInt c3;
	int c_i3 = 0;
	for (c3=0;c3<vec_size3;c3+=1){
		VecGetValues(exact_sol,1,&c3,&temp3);
		double v3 = (double)temp3;
		masa_solution.push_back(v3);
		c_i3+=1;
	}
	return ierr;
}

