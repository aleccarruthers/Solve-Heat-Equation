#include<iostream>
#include "matrices.h"
#include "masa_arrays.h"
#include "numerical_solve.h"
#include "error_norm.h"
#include "parameter_check.h"
#include "node_val_check.h"
#include "parser.h"
#include "hdf5.h"
#include "hd_prep.h"
#include "hdf5_output.h"
#include<grvy.h>
#include<iomanip>
#include<vector>
#include "petsc_error.h"
#if __has_include("petsc.h")
 #define USE_PETSc
 #include "petsc.h"
 #include "petsc_setup.h"
 #include "petsc_masa.h"
 #include "petsc_solve.h"
#endif

using namespace GRVY;
using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;

int main(int argc, char* argv[]){
	
	// Check for Missing input file
	if (argc==1){
		cout << "Missing input file" << endl;
		exit (EXIT_FAILURE);}

	//////////////////////////////////////////////////
	///////// Initialize important variables /////////
	//////////////////////////////////////////////////
	int dimension;
	int finorder, max_iterations;
	std::string solver, logger, nodes;
	std::string verification;
	double length, frequency;
	double K;
	int N;
	// Masa Sources
	std::vector<double> sources;
	// Masa Solution
	std::vector<double> solution;
	// Numerical Solution
	std::vector<double> num_sol_c;
	
	grvy_timer_init("GRVY Timing");

	//////////////////////////////////////////////////
	//////////////// Input Parsing ///////////////////
	//////////////////////////////////////////////////
	grvy_timer_begin("Input Parsing");
	parse_inputs(dimension,finorder,max_iterations,solver,logger,nodes,verification,length,frequency,K);
	int input_check = params_check(dimension,finorder,solver,length,frequency,logger,max_iterations);
	if (input_check!=0){
		exit(input_check);}
	N = node_check(dimension,nodes);
	grvy_timer_end("Input Parsing");
	cout << std::setprecision(12);


	//////////////////////////////////////////////////
	//////////////// PETSc Solve /////////////////////
	//////////////////////////////////////////////////
	if (solver=="petsc"){
		#ifdef USE_PETSc
		// PETSc Initiation
		PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
		Mat A_p; // A Matrix
		Vec src, masa_sol, num_sol_pet; // Masa Source, Exact Sol, Num. Sol

		// Build A Matrix
		grvy_timer_begin("PETSc: Initialize A Matrix");
		petsc_set(N,dimension,finorder, A_p,logger);
		grvy_timer_end("PETSc: Initialize A Matrix");
		
		// Initialize MASA sources and solution vectors
		grvy_timer_begin("PETSc: Initialize MASA Source and Exact Solution");
		petsc_masa_source(N,dimension,finorder,length,frequency,K,logger,src,masa_sol,num_sol_pet);
		grvy_timer_end("PETSc: Initialize MASA Source and Exact Solution");
		

		// PETSc solve
		grvy_timer_begin("PETSc: Numerical Solve");
		pet_solve(A_p,src,masa_sol,N,num_sol_c,sources,solution);
		grvy_timer_end("PETSc: Numerical Solve");
		PetscFinalize();
		#endif
		#ifndef USE_PETSc
		int e = petsc_error();
		exit(e);
		#endif
	}

	//////////////////////////////////////////////////
	//////////////// Non-PETSc Solve /////////////////
	//////////////////////////////////////////////////
	else{
		// Initialize A Matrix
		grvy_timer_begin("Initialize A Matrix");
		std::vector<std::vector<double>> A = A_matrix(N,dimension,finorder,logger);
		grvy_timer_end("Initialize A Matrix");
		
		// Initialize MASA heat sources
		grvy_timer_begin("Initialize MASA Sources");
		sources = masa_sources(N,dimension,finorder,length,frequency,K,logger);
		grvy_timer_end("Initialize MASA Sources");

		// Initialize MASA solution
		grvy_timer_begin("Initialize MASA Solution");
		solution = masa_solution(N,dimension,length,frequency,K,logger);
		grvy_timer_end("Initialize MASA Solution");
	
		// Solve system of equations
		grvy_timer_begin("Numerical Solver");
		num_sol_c = num_solve(solver,A, sources,dimension,finorder,length,frequency,K,max_iterations);
		grvy_timer_end("Numerical Solver");
	}


	//////////////////////////////////////////////////
	////////////////// L2 Loss ///////////////////////
	//////////////////////////////////////////////////
	grvy_timer_begin("L2 Norm");
	int final_size = num_sol_c.size();
	double error_est = l2Norm(final_size,num_sol_c,solution);
	grvy_timer_end("L2 Norm");
	cout << "L2 Norm: " << error_est << endl;
	
	//////////////////////////////////////////////////
	////////////////// HDF5 Output ///////////////////
	//////////////////////////////////////////////////
	grvy_timer_begin("HDF5 Formatting");
	// Define arrays to hold important output vectors
	double hd_sol[final_size];
	double hd_num_sol[final_size];
	double hd_l2[1];
	double coord_x[final_size];
	double coord_y[final_size];
	// Formatting vectors for hdf5
	hdf5_prep(solution, num_sol_c, hd_sol, hd_num_sol,error_est,hd_l2,dimension,coord_x,coord_y,length);
	// Write HDF5
	hdf5_write(hd_num_sol,hd_sol,hd_l2,final_size,coord_x,coord_y,dimension);
	grvy_timer_end("HDF5 Formatting");
	grvy_timer_finalize();
	grvy_timer_summarize();
	return 0;	
}
