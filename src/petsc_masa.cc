#include<iostream>
#include<vector>
#include<cmath>
#include "petsc.h"
#include<masa.h>
#include<iomanip>
#include<float.h>
using std::cout;
using std::endl;
using std::vector;
using std::string;

PetscErrorCode petsc_masa_source(int N, int dimension, int finorder, double length, double frequency, double K,std::string loglevel, Vec &masa_src_petsc, Vec &masa_sol_petsc, Vec &num_sol){
	PetscErrorCode ierr;
	cout << std::setprecision(12);
	// Masa Source vector for 2nd order FD with 1-D
	if (finorder==2 && dimension==1){
		// Coefficient on FD term
		double coeff = -1;
		double L = length;
		// Set Masa Params
		int model = masa_init("1D Heat Conduction","heateq_1d_steady_const");
		masa_set_param("k_0",K);
		double K_0 = masa_get_param("k_0");
		masa_set_param("A_x",frequency);
		// Spacing between nodes
		double delx = L/(N-1);
		
		// Initialize PETSc source vector
		ierr = VecCreate(PETSC_COMM_WORLD,&masa_src_petsc);CHKERRQ(ierr);
		ierr = VecSetSizes(masa_src_petsc,PETSC_DECIDE,N);CHKERRQ(ierr);
		ierr = VecSetFromOptions(masa_src_petsc);CHKERRQ(ierr);
		ierr = VecDuplicate(masa_src_petsc,&num_sol);CHKERRQ(ierr);
		// Initialize PETSc solution vector
		ierr = VecCreate(PETSC_COMM_WORLD,&masa_sol_petsc);CHKERRQ(ierr);
		ierr = VecSetSizes(masa_sol_petsc,PETSC_DECIDE,N);CHKERRQ(ierr);
		ierr = VecSetFromOptions(masa_sol_petsc);CHKERRQ(ierr);
		for (int i=0;i<N;i++){
			double x = i*delx;
			PetscInt count;
			count = i;
			// Use the exact solution for the boundary nodes
			if (i==0 || i==N-1){
				double ex_val = masa_eval_1d_exact_t(x);
				PetscScalar src[1];
				src[0] = ex_val;
				ierr = VecSetValues(masa_src_petsc,1,&count,src,INSERT_VALUES);CHKERRQ(ierr);
				ierr = VecSetValues(masa_sol_petsc,1,&count,src,INSERT_VALUES);CHKERRQ(ierr);
				
			}
			// Interior nodes with a multiplication of the coefficient on the finite difference
			// term
			else{
				double src_val = coeff*pow(delx,2)*masa_eval_1d_source_t(x)/K_0;
				double sol_val = masa_eval_1d_exact_t(x);
				PetscScalar sol2[1];
				PetscScalar src2[1];
				src2[0] = src_val;
				sol2[0] = sol_val;
				ierr = VecSetValues(masa_src_petsc,1,&count,src2,INSERT_VALUES);CHKERRQ(ierr);
				ierr = VecSetValues(masa_sol_petsc,1,&count,sol2,INSERT_VALUES);CHKERRQ(ierr);
			}
		}
		// Assemble PETSc Vector
		ierr = VecAssemblyBegin(masa_src_petsc);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(masa_src_petsc);CHKERRQ(ierr);
		// Assemble Solution Vector
		ierr = VecAssemblyBegin(masa_sol_petsc);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(masa_sol_petsc);CHKERRQ(ierr);
		
	}
	else if (finorder==2 && dimension==2){
		// Initialize 2D MASA
		int model = masa_init("2D Heat Conduction","heateq_2d_steady_const");
		masa_init_param();
		// The number of nodes along each axis
		int N_x = sqrt(N);
		
		double coeff = -1;
		// Set MASA parameters based off inputs from the input file
		masa_set_param("k_0",K);
		double K_0 = masa_get_param("k_0");
		masa_set_param("A_x",frequency);
		masa_set_param("B_y",frequency);
		PetscInt d2_size = N;
		//Spacing between nodes in both dimensions
		double delx = length/(N_x-1);
		// Increment for each iteration
		PetscInt increment = 0;
		// Initialize PETSc vector
		ierr = VecCreate(PETSC_COMM_WORLD,&masa_src_petsc);CHKERRQ(ierr);
		ierr = VecSetSizes(masa_src_petsc,PETSC_DECIDE,d2_size);CHKERRQ(ierr);
		ierr = VecSetFromOptions(masa_src_petsc);CHKERRQ(ierr);
		ierr = VecDuplicate(masa_src_petsc,&num_sol);CHKERRQ(ierr);
		
		// Initialize PETSc solution vector
		ierr = VecCreate(PETSC_COMM_WORLD,&masa_sol_petsc);CHKERRQ(ierr);
		ierr = VecSetSizes(masa_sol_petsc,PETSC_DECIDE,N);CHKERRQ(ierr);
		ierr = VecSetFromOptions(masa_sol_petsc);CHKERRQ(ierr);
		
		// Treating the (0,0) point as the bottom left node in the grid, so that
		// is what the first term in the vector corresponds to.
		for (int row=0;row<N_x;row++){
			// Current nodes row height
			double y = delx*row;
			for (int col=0;col<N_x;col++){
				// Current x position of node
				double x = delx*col;
				// Set Solution vector at all iteration steps
				double sol_term = masa_eval_2d_exact_t(x,y);
				PetscScalar sol_row;
				sol_row = sol_term;
				ierr = VecSetValues(masa_sol_petsc,1,&increment,&sol_row,INSERT_VALUES);CHKERRQ(ierr);
				// Top and Bottom boundaries of grid set to exact solution
				if (row==0 || row==N_x-1){
					PetscScalar ex_value[1];
					double bottom_exact = masa_eval_2d_exact_t(x,y);
					ex_value[0]=bottom_exact;
					ierr = VecSetValues(masa_src_petsc,1,&increment,ex_value,INSERT_VALUES);CHKERRQ(ierr);}
				// Left and Right boundaries of grid set to exact solution
				else if (col==0 || col==N_x-1){
					PetscScalar ex_value_side[1];
					double side_exact = masa_eval_2d_exact_t(x,y);
					ex_value_side[0]=side_exact;
					ierr = VecSetValues(masa_src_petsc,1,&increment,ex_value_side,INSERT_VALUES);CHKERRQ(ierr);}
				// Interior nodes of grid
				else{
					PetscScalar src_terms[1];
					double inner = masa_eval_2d_source_t(x,y);
					double inner_coeff = (coeff*inner*(pow(delx,2)))/K_0; 
					src_terms[0] = inner_coeff;
					ierr = VecSetValues(masa_src_petsc,1,&increment,src_terms,INSERT_VALUES);CHKERRQ(ierr);}
				increment+=1;
			}
		}
		// Assemble Source Vector
		ierr = VecAssemblyBegin(masa_src_petsc);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(masa_src_petsc);CHKERRQ(ierr);
		// Assemble Solution Vector
		ierr = VecAssemblyBegin(masa_sol_petsc);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(masa_sol_petsc);CHKERRQ(ierr);
	}
				
	else if (finorder==4 && dimension==1){

		// Total nodes increases by two because of the ghost cells on the left and right sides
		int N_tot = N + 2;
		// Coefficient on FD term
		double coeff = -12;
		// Set Masa Params
		int model = masa_init("1D Heat Conduction","heateq_1d_steady_const");
		masa_set_param("k_0",K);
		double K_0 = masa_get_param("k_0");
		masa_set_param("A_x",frequency);
		// Spacing between nodes
		// This distance should stay the same, despite the addition of ghost nodes
		// They are not technically in the grid, so they don't contribute to this distance.
		double delx = length/(N-1);
		
		// Initialize PETSc vector
		ierr = VecCreate(PETSC_COMM_WORLD,&masa_src_petsc);CHKERRQ(ierr);
		ierr = VecSetSizes(masa_src_petsc,PETSC_DECIDE,N_tot);CHKERRQ(ierr);
		ierr = VecSetFromOptions(masa_src_petsc);CHKERRQ(ierr);
		ierr = VecDuplicate(masa_src_petsc,&num_sol);CHKERRQ(ierr);
		// Initialize PETSc solution vector
		ierr = VecCreate(PETSC_COMM_WORLD,&masa_sol_petsc);CHKERRQ(ierr);
		ierr = VecSetSizes(masa_sol_petsc,PETSC_DECIDE,N_tot);CHKERRQ(ierr);
		ierr = VecSetFromOptions(masa_sol_petsc);CHKERRQ(ierr);
		for (int i=-1;i<N+1;i++){
			// Spacing between nodes
			double x = i*delx;
			// Petsc vector index
			PetscInt count;
			count = i+1;
			// Solution vector
			double sol_row = masa_eval_1d_exact_t(x);
			PetscScalar sol_term[1];
			sol_term[0] = sol_row;
			ierr = VecSetValues(masa_sol_petsc,1,&count,sol_term,INSERT_VALUES);CHKERRQ(ierr);

			// Use the exact solution for the boundary nodes
			if (i<=0 || i>=N-1){
				double ex_val = masa_eval_1d_exact_t(x);
				PetscScalar src[1];
				src[0] = ex_val;
				ierr = VecSetValues(masa_src_petsc,1,&count,src,INSERT_VALUES);CHKERRQ(ierr);
			}
			// Interior nodes with a multiplication of the coefficient on the finite difference
			// term
			else{
				double src_val = coeff*pow(delx,2)*masa_eval_1d_source_t(x)/K_0;
				PetscScalar src2[1];
				src2[0] = src_val;
				ierr = VecSetValues(masa_src_petsc,1,&count,src2,INSERT_VALUES);CHKERRQ(ierr);
			}
		}
		// Assemble PETSc Vector
		ierr = VecAssemblyBegin(masa_src_petsc);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(masa_src_petsc);CHKERRQ(ierr);
		ierr = VecAssemblyBegin(masa_sol_petsc);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(masa_sol_petsc);CHKERRQ(ierr);
	}
	else if (finorder==4 && dimension==2){
		// Initialize 2D MASA
		int model = masa_init("2D Heat Conduction","heateq_2d_steady_const");
		masa_init_param();
		// The number of nodes along each axis
		int N_x = sqrt(N);
		int N_x_temp = sqrt(N)+2;
		int N_tot = pow(N_x_temp,2);
		double coeff = -12;
		// Set MASA parameters based off inputs from the input file
		masa_set_param("k_0",K);
		double K_0 = masa_get_param("k_0");
		masa_set_param("A_x",frequency);
		masa_set_param("B_y",frequency);
		
		//Spacing between nodes in both dimensions
		double delx = length/(N_x-1);
		// Increment for each iteration
		PetscInt increment = 0;

		// Initialize PETSc vector
		ierr = VecCreate(PETSC_COMM_WORLD,&masa_src_petsc);CHKERRQ(ierr);
		ierr = VecSetSizes(masa_src_petsc,PETSC_DECIDE,N_tot);CHKERRQ(ierr);
		ierr = VecSetFromOptions(masa_src_petsc);CHKERRQ(ierr);
		ierr = VecDuplicate(masa_src_petsc,&num_sol);CHKERRQ(ierr);
		
		// Initialize PETSc solution vector
		ierr = VecCreate(PETSC_COMM_WORLD,&masa_sol_petsc);CHKERRQ(ierr);
		ierr = VecSetSizes(masa_sol_petsc,PETSC_DECIDE,N_tot);CHKERRQ(ierr);
		ierr = VecSetFromOptions(masa_sol_petsc);CHKERRQ(ierr);
		// Treating the (0,0) point as the bottom left node in the grid, so that
		// is what the first term in the vector corresponds to.
		for (int row=-1;row<N_x+1;row++){
			// Current nodes row height
			double y = delx*row;
			for (int col=-1;col<N_x+1;col++){
				// Current x position of node
				double x = delx*col;

				// Solution Vector fill
				double sol_row = masa_eval_2d_exact_t(x,y);
				PetscScalar sol_term[1];
				sol_term[0] = sol_row;
				ierr = VecSetValues(masa_sol_petsc,1,&increment,sol_term,INSERT_VALUES);CHKERRQ(ierr);
				// Top and Bottom boundaries of grid set to exact solution
				if (row<=0 || row>=N_x-1){
					PetscScalar ex_value[1];
					double bottom_exact = masa_eval_2d_exact_t(x,y);
					ex_value[0]=bottom_exact;
					ierr = VecSetValues(masa_src_petsc,1,&increment,ex_value,INSERT_VALUES);CHKERRQ(ierr);}
				// Left and Right boundaries of grid set to exact solution
				else if (col<=0 || col>=N_x-1){
					PetscScalar ex_value_side[1];
					double side_exact = masa_eval_2d_exact_t(x,y);
					ex_value_side[0]=side_exact;
					ierr = VecSetValues(masa_src_petsc,1,&increment,ex_value_side,INSERT_VALUES);CHKERRQ(ierr);}
				// Interior nodes of grid
				else{
					PetscScalar src_terms[1];
					double inner = masa_eval_2d_source_t(x,y);
					double inner_coeff = (coeff*inner*(pow(delx,2)))/K_0; 
					src_terms[0] = inner_coeff;
					ierr = VecSetValues(masa_src_petsc,1,&increment,src_terms,INSERT_VALUES);CHKERRQ(ierr);}
				increment+=1;
			}
		}
		ierr = VecAssemblyBegin(masa_src_petsc);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(masa_src_petsc);CHKERRQ(ierr);
		ierr = VecAssemblyBegin(masa_sol_petsc);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(masa_sol_petsc);CHKERRQ(ierr);
	}
	return ierr;

}

