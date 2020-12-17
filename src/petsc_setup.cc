#include<iostream>
#include "petsc.h"
#include<cmath>
#include<vector>
#include<iomanip> 
using std::cout;
using std::endl;
using std::vector;
using std::string;
PetscErrorCode petsc_set(int N,int dimension, int finorder, Mat &A_p, std::string mode){
	// Matrix Build //
	PetscErrorCode ierr;
	MPI_Comm comm;
	cout << std::setprecision(12);
	///////////////////////////////////////////////////////////////////////////
	////////////////////////////////  1D -2nd Order  //////////////////////////
	///////////////////////////////////////////////////////////////////////////
	if (dimension==1 && finorder==2){
		PetscInt pet_i, pet_j, pet_dims = N;
		PetscInt col[3];
		PetscInt bc_col[1];
		pet_i = 0; pet_j = 0;
		ierr = MatCreate(PETSC_COMM_WORLD,&A_p);CHKERRQ(ierr);
		ierr = MatSetSizes(A_p,PETSC_DECIDE,PETSC_DECIDE,pet_dims,pet_dims);CHKERRQ(ierr);
		ierr = MatSetUp(A_p);CHKERRQ(ierr);
		PetscScalar bc[1];
		PetscScalar inner[3];
		PetscScalar bc_val = 1;
		bc[0] = 1.0;
		inner[0] = 1.0; inner[1] = -2.0; inner[2] = 1.0;
		for (int i=0;i<N;i++){
			col[0] = pet_i-1; col[1] = pet_i; col[2] = pet_i+1;
			bc_col[0] = pet_i;
			if (i==0){
				ierr = MatSetValues(A_p,1,&pet_i,1,&pet_i,&bc_val,INSERT_VALUES);
				CHKERRQ(ierr);}
			else if (i==N-1){
				ierr = MatSetValues(A_p,1,&pet_i,1,&pet_i,&bc_val,INSERT_VALUES);
				CHKERRQ(ierr);}
			else{
				ierr = MatSetValues(A_p,1,&pet_i,3,col,inner,INSERT_VALUES);
				CHKERRQ(ierr);}
			pet_i+=1;	
		}
		ierr = MatAssemblyBegin(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		/*if (mode=="debug"){
			ierr = MatView(A_p,0);CHKERRQ(ierr);
		}*/
	}

	///////////////////////////////////////////////////////////////////////////
	////////////////////////////////  1D -4th Order  //////////////////////////
	///////////////////////////////////////////////////////////////////////////
	else if (dimension==1 && finorder==4){
		// Increase the dimensions of the matrix to account for ghost cells
		int N_x = N+2;
		PetscInt pet_i, pet_j, pet_dims = N_x;
		PetscInt col[5];
		pet_i = 0; pet_j = 0;
		ierr = MatCreate(PETSC_COMM_WORLD,&A_p);CHKERRQ(ierr);
		ierr = MatSetSizes(A_p,PETSC_DECIDE,PETSC_DECIDE,pet_dims,pet_dims);CHKERRQ(ierr);
		ierr = MatSetUp(A_p);CHKERRQ(ierr);
		PetscScalar inner[5];
		PetscScalar bc_val = 1.0;
		inner[0] = -1.0; inner[1] = 16.0; inner[2] = -30.0; inner[3] = 16.0; inner[4] = -1.0;
		for (int i=0;i<N_x;i++){
			col[0] = pet_i-2; col[1] = pet_i-1; col[2] = pet_i;
			col[3] = pet_i+1; col[4] = pet_i+2;
			if (i<=1){
				ierr = MatSetValues(A_p,1,&pet_i,1,&pet_i,&bc_val,INSERT_VALUES);
				CHKERRQ(ierr);}
			else if (i>=N_x-2){
				ierr = MatSetValues(A_p,1,&pet_i,1,&pet_i,&bc_val,INSERT_VALUES);
				CHKERRQ(ierr);}
			else{
				ierr = MatSetValues(A_p,1,&pet_i,5,col,inner,INSERT_VALUES);
				CHKERRQ(ierr);}
			pet_i+=1;	
		}
		ierr = MatAssemblyBegin(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		/*if (mode=="debug"){
			ierr = MatView(A_p,0);CHKERRQ(ierr);
		}*/

	}
	///////////////////////////////////////////////////////////////////////////
	////////////////////////////////  2D -2nd Order  //////////////////////////
	///////////////////////////////////////////////////////////////////////////
	else if (dimension==2 && finorder==2){
		int N_x = sqrt(N);
		int threshold = N - N_x;
		PetscInt pet_i, pet_j, pet_dims = N;
		pet_i = 0; pet_j = 0;
		ierr = MatCreate(PETSC_COMM_WORLD,&A_p);CHKERRQ(ierr);
		ierr = MatSetSizes(A_p,PETSC_DECIDE,PETSC_DECIDE,pet_dims,pet_dims);CHKERRQ(ierr);
		ierr = MatSetUp(A_p);CHKERRQ(ierr);
		PetscScalar bc_val = 1.0;
		// Top Boundary
		for (int i=0;i<N;i++){
			if (i<N_x || i>=threshold){
				ierr = MatSetValues(A_p,1,&pet_i,1,&pet_i,&bc_val,INSERT_VALUES);
				CHKERRQ(ierr);}
			pet_i+=1;
		}
		
		// Side Boundaries
		for (int bc=N_x;bc<threshold;bc+=N_x){
			PetscInt bc_left = bc;
			PetscInt bc_right = bc + N_x -1;
			ierr = MatSetValues(A_p,1,&bc_left,1,&bc_left,&bc_val,INSERT_VALUES);
			CHKERRQ(ierr);
			ierr = MatSetValues(A_p,1,&bc_right,1,&bc_right,&bc_val,INSERT_VALUES);
			CHKERRQ(ierr);
		}
		// Inner Nodes
		PetscScalar inner[5];
		inner[0] = 1.0; inner[1] = 1.0; inner[2] = -4.0; inner[3] = 1.0; inner[4] = 1.0;
		PetscInt col[5];
		int node_int = N_x + 1;
		int num_non_bc = N_x - 2;
		PetscInt N_per_dim = N_x;
		for (int vals=node_int;vals<threshold;vals+=N_x){
			for (int g = 0;g<num_non_bc;g++){
				PetscInt pos = vals+g;
				col[0] = pos - N_per_dim;col[1] = pos - bc_val; col[2] = pos;
				col[3] = pos + bc_val; col[4] = pos+N_per_dim;
				ierr = MatSetValues(A_p,1,&pos,5,col,inner,INSERT_VALUES);
				CHKERRQ(ierr);
			}
		}
		ierr = MatAssemblyBegin(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		/*if (mode=="debug"){
			ierr = MatView(A_p,0);CHKERRQ(ierr);
		}*/
	}
	///////////////////////////////////////////////////////////////////////////
	////////////////////////////////  2D -4th Order  //////////////////////////
	///////////////////////////////////////////////////////////////////////////
	else if (dimension==2 && finorder==4){
		int N_x_temp = sqrt(N)+2;
		int N_tot = pow(N_x_temp,2);
		int threshold = N_tot - 2*N_x_temp;
		PetscInt pet_i, pet_j, pet_dims = N_tot;
		pet_i = 0; pet_j = 0;
		ierr = MatCreate(PETSC_COMM_WORLD,&A_p);CHKERRQ(ierr);
		ierr = MatSetSizes(A_p,PETSC_DECIDE,PETSC_DECIDE,pet_dims,pet_dims);CHKERRQ(ierr);
		ierr = MatSetUp(A_p);CHKERRQ(ierr);
		PetscScalar bc_val = 1.0;
		// Top rows of A (bottom boundary of grid)
		for (int i=0;i<2*N_x_temp;i++){
			ierr = MatSetValues(A_p,1,&pet_i,1,&pet_i,&bc_val,INSERT_VALUES);
			CHKERRQ(ierr);
			pet_i+=1;
		}
		// Bottom rows of A (top boundary of grid)
		pet_i=threshold;
		for (int i=threshold;i<N_tot;i++){
			ierr = MatSetValues(A_p,1,&pet_i,1,&pet_i,&bc_val,INSERT_VALUES);
			CHKERRQ(ierr);
			pet_i+=1;
		}
		// Side Boundaries
		for (int bc=2*N_x_temp;bc<threshold;bc+=N_x_temp){
			PetscInt bc_left = bc;// Ghost column on left
			PetscInt bc_right = bc + N_x_temp -1; // Ghost column on right
			PetscInt bc_left_real = bc+1;// Ghost column on left
			PetscInt bc_right_real = bc + N_x_temp -2; // Ghost column on right

			// Set Ghost Cells on left Boundary
			ierr = MatSetValues(A_p,1,&bc_left,1,&bc_left,&bc_val,INSERT_VALUES);
			CHKERRQ(ierr);
			// Set Ghost Cells on right Boundary
			ierr = MatSetValues(A_p,1,&bc_right,1,&bc_right,&bc_val,INSERT_VALUES);
			CHKERRQ(ierr);
			// Set True boundary condition on left
			ierr = MatSetValues(A_p,1,&bc_left_real,1,&bc_left_real,&bc_val,INSERT_VALUES);
			CHKERRQ(ierr);
			// Set True boundary condition on right
			ierr = MatSetValues(A_p,1,&bc_right_real,1,&bc_right_real,&bc_val,INSERT_VALUES);
			CHKERRQ(ierr);
		}
		// Inner Nodes
		PetscScalar inner[9];
		inner[0] = -1.0; inner[1] = 16.0; inner[2] = -1.0; inner[3] = 16.0; inner[4] = -60.0;
		inner[5] = 16.0; inner[6] = -1.0; inner[7] = 16.0; inner[8] = -1.0;
		PetscInt col[9];
		PetscInt scale1 = 1;
		PetscInt scale2 = 2;
		// First non-bc node
		int node_int = 2*N_x_temp + 2;
		// Number of nodes before another boundary is reached
		int num_non_bc = N_x_temp - 4;
		PetscInt N_per_dim = N_x_temp;
		for (int vals=node_int;vals<threshold;vals+=N_x_temp){
			for (int g = 0;g<num_non_bc;g++){
				PetscInt pos = vals+g;
				col[0] = pos - scale2*N_per_dim;col[1] = pos - N_per_dim; col[2] = pos-scale2;
				col[3] = pos - scale1; col[4] = pos; col[5] = pos+scale1;
				col[6] = pos+scale2; col[7] = pos+N_per_dim; col[8] = pos+scale2*N_per_dim;
				ierr = MatSetValues(A_p,1,&pos,9,col,inner,INSERT_VALUES);
				CHKERRQ(ierr);
			}
		}
		ierr = MatAssemblyBegin(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		/*if (mode=="debug"){
			ierr = MatView(A_p,0);CHKERRQ(ierr);
		}*/
	}
	return ierr;

}

