#include<iostream>
#include<iomanip>
#include<vector>
#include "numerical_solve.h"
#include<masa.h>
#include<cmath>
using std::cout;
using std::endl;
using std::vector;
using std::string;


std::vector<double> num_solve(std::string method,std::vector<std::vector<double>> A, std::vector<double> sources, int dimension, int order,double length,double frequency,double K,int max_iterations){
	// Setting initial parameters
	cout << std::setprecision(12);
	int vector_length = sources.size();
	double L = length;
	double delx;
	int N_x;
	double freq = frequency;

	// Error tolerance and initial parameters to enter solver while loop
	double tol = 1e-12;
	double diff = 10.0;
	double diff2 = 10.0;

	// Defining step size based on dimension.
	if (dimension==1){
		N_x = vector_length;
		delx = L/(vector_length-1);}
	else{
		N_x = sqrt(vector_length);
		delx = L/(N_x - 1);
	}

	////////////////////////////////////////////////////////////////////////////////////////
	// Jacobi Iteration for 2nd order case (1D or 2D), i.e. no ghost cells
	////////////////////////////////////////////////////////////////////////////////////////
	if (method=="Jacobi" && order==2){
		cout << std::setprecision(12);
		// Initialize solution vector and initial guess/previous iteration vector
		std::vector<double> x_sol(vector_length,0.0);
		std::vector<double> x_guess(vector_length,0.1);

		int counter = 0;
		while (diff>tol){
			diff = 0.0;
			diff2 = 0;
			if (counter >= max_iterations){
				cout << "*** Maximum Number of Iterations Reached ***" << endl;
				fprintf( stderr, "*** Maximum Number of Iterations Reached ***\n");
				exit(1);}
			// Loop over all columns for a particular row of the A matrix
			for (int i=0;i<vector_length;i++){
				// initialize summation parameters
				double val;
				double inner_sum = 0;
				// For all columns in the row, multiply the coefficient in A by the
				// corresponding entry in the guessed vector and sum the result
				for (int j = 0;j<vector_length;j++){
					if (i!=j){
						val = A[i][j]*x_guess[j];
						inner_sum = inner_sum+val;
					}
				}
				// Once all nodes have been accounted for in the row, add the source
				// term for the corresponding row entry to the total 'inner_sum'
				val = -1*inner_sum + sources[i];
				// Divide by the column coefficient in A with index i
				x_sol[i] = (1/A[i][i])*val;
				// Accumulated L2 error
				diff = diff + pow(x_sol[i]-x_guess[i],2);
			}
			// Normalize error
			diff = sqrt(diff/vector_length);
			// Check error between solution and initial vector. If small enough, break
			// the loop. If not, set the initial vector equal to the computed solution 
			// vector and redo the entire computation
			if (diff>tol){
				for (int k=0;k<vector_length;k++){
					x_guess[k] = x_sol[k];
				}
			}
			counter+=1;
		}
		return x_sol;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////
	// Gauss-Seidel for 2nd Order cases (1D/2D)
	////////////////////////////////////////////////////////////////////////////////////////
	else if (method=="GS" && order==2){
		cout << std::setprecision(12);
		// Initialize solution and initial guess vectors
		std::vector<double> x_sol(vector_length,0.1);
		std::vector<double> x_guess(vector_length,0.1);
		// Counter for max iterations and set error tolerance
		int counter = 0;
		tol = 1e-12;

		while (diff>tol){
			diff = 0.0;
			diff2 = 0;
			if (counter >= max_iterations){
				cout << "*** Maximum Number of Iterations Reached ***" << endl;
				fprintf( stderr, "*** Maximum Number of Iterations Reached ***\n");
				exit(1);}
			for (int i=0;i<vector_length;i++){
				// Initialize value and inner sum terms for both the loop that checks
				// the most recently updated terms and the (1) and the loop that 
				// considers nodes ahead of its current position (2)
				double val=0;
				double val2;
				double inner_sum = 0;
				double inner_sum2 = 0;
				// (1) Adding most recently updated terms within the current iteration
				for (int j = 0;j<i;j++){
					val = A[i][j]*x_sol[j];
					inner_sum = inner_sum+val;
				}
				// (2) All variables ahead of the current variable are considered
				for (int k = i+1;k<vector_length;k++){
					val2 = A[i][k]*x_guess[k];
					inner_sum2 = inner_sum2+val2;
				}
				// Combine the summation from both for loops and add the sources.
				// Equivalent to multiplying A by the estimated solution vector,
				// summing all products besides the current node and adding the 
				// corresponding source term
				double val_tot = -1*inner_sum - inner_sum2 + sources[i];
				// Divide by the current nodes coefficient
				x_sol[i] = (1/A[i][i])*val_tot;
			}
			// Compute Loss
			for (int n = 0;n<vector_length;n++){
				diff=diff+pow(x_sol[n] - x_guess[n],2);
			}
			// Normalize and square root the accumulated loss
			diff = sqrt(diff/vector_length);
			// Check if the resulting error is above or below the threshold
			// If greater than the threshold, repeat the loop with the solution 
			// vector of this iteration becoming the initial vector for the next.
			if (diff>tol){
				for (int m=0;m<vector_length;m++){
					x_guess[m] = x_sol[m];}
			}
			counter = counter + 1;
		}
		return x_sol;
	}

	////////////////////////////////////////////////////////////////////////////////////////
	// Gauss-Seidel for 4th Order cases (1D/2D)
	////////////////////////////////////////////////////////////////////////////////////////

	else if (method=="GS" && order==4){
		cout << std::setprecision(12);
		// Initialize solution and initial guess/previous iteration vectors
		std::vector<double> x_sol(vector_length,0.0);
		std::vector<double> x_guess(vector_length,0.1);
		int counter = 0;
		
		// ----------------------------- //
		/////////////// 1-D ///////////////
		// ----------------------------- //

		// Handling ghost cells by querying MASA for the exact temperature
		// at nodes outside the domain. This is done for all nodes with an 
		// incomplete stencil
		if (dimension==1){
			// Initialize MASA and set parameters
			int model = masa_init("1D","heateq_1d_steady_const");
			masa_init_param();
			masa_set_param("A_x",freq);
			masa_set_param("k_0",K);
			while (diff>tol){
				diff = 0.0;
				diff2 = 0;
				// Max iterations check
				if (counter >= max_iterations){
					cout << "*** Maximum Number of Iterations Reached ***" << endl;
					fprintf( stderr, "*** Maximum Number of Iterations Reached ***\n");
					exit(1);}
				for (int i=0;i<vector_length;i++){
					// Using the same variables to keep track of inner product
					// terms
					double val;
					double val2;
					double inner_sum = 0;
					double inner_sum2 = 0;
					// Variables to store the MASA exact solution for nodes
					// with an incomplete stencil. Otherwise they contribute zero
					double heat_ghost1=0.0;
					double heat_ghost2=0.0;

					// Adding most recently updated terms within the
					// current iteration
					for (int j = 0;j<i;j++){
						val = A[i][j]*x_sol[j];
						inner_sum = inner_sum+val;
					}

					// All variables ahead of the current variable are considered
					for (int k = i+1;k<vector_length;k++){
						val2 = A[i][k]*x_guess[k];
						inner_sum2 = inner_sum2+val2;
					}
					
					///////////////////////////////////////////////////////////
					//////////////// Accounting for Ghost Cells ///////////////
					///////////////////////////////////////////////////////////

					// Second node from the start has incomplete stencil in 
					// the negative x direction.
					if (i==1){
						heat_ghost1 = masa_eval_1d_exact_t(-1*i*delx);}

					// Second to last node has an incomplete stencil in the 
					// positive x direction
					if (i==vector_length-2){
						heat_ghost2 = masa_eval_1d_exact_t(vector_length*delx);}

					// Ghost cell terms are added to the inner product terms 
					// because that is what would have occured in the system
					// of equations had the nodes existed
					double val_tot = -1*inner_sum - inner_sum2 + heat_ghost1 + heat_ghost2 + sources[i];
					x_sol[i] = (1/A[i][i])*val_tot;
				}
				// Compute Loss
				for (int n = 0;n<vector_length;n++){
					diff=diff+(pow(x_sol[n] - x_guess[n],2));
				}
				// Normalize Loss
				diff = sqrt(diff/vector_length);
				// Check the error of the solution. If small enough, exit the loop.
				// If too large, set the solution vector of this iteration equal 
				// to the initial vector for the next iteration.
				if (diff>tol){
					for (int m=0;m<vector_length;m++){
						x_guess[m] = x_sol[m];}
				}
				counter = counter + 1;
				}
			return x_sol;
		}

		// ----------------------------- //
		/////////////// 2-D ///////////////
		// ----------------------------- //
		// Similar to what was done in the 1D case, MASA will be used to generate the 
		// exact temperature for nodes outside the domain causing an incomplete stencil.
		// Unlike the 1D case, these incomplete stencils occur at more than just two nodes.
		// The exact solutions will be added to the inner sum of the appropriate nodes as
		// if the node actually existed in the domain. Since the coefficient on the ghost cells
		// from the finite differencing method is negative 1, they add to the inner sum
		else if (dimension==2){
			// Initialize MASA and its parameters
			int model = masa_init("2D","heateq_2d_steady_const");
			masa_init_param();
			masa_set_param("A_x",freq);
			masa_set_param("B_y",freq);
			masa_set_param("k_0",K);

			int counter = 0;
			while (diff>tol){
				diff = 0.0;
				diff2 = 0;
				if (counter >= max_iterations){
					cout << "*** Maximum Number of Iterations Reached ***" << endl;
					fprintf( stderr, "*** Maximum Number of Iterations Reached ***\n");
					exit(1);}
				for (int i=0;i<vector_length;i++){
					// Initialize inner sum variables
					double val;
					double val2;
					double inner_sum=0;	
					double inner_sum2=0;
					// Initialize MASA exact solution variables that only contribute
					// to the solution when a node with an incomplete stencil is
					// encountered
					double heat_ghost1=0.0;
					double heat_ghost2=0.0;

					// Adding most recently updated terms within current i iteration
					for (int j=0;j<i;j++){
						val = A[i][j]*x_sol[j];
						inner_sum = inner_sum+val;}
					// All variables ahead of the current iteration
					for (int k=i+1;k<vector_length;k++){
						val2 = A[i][k]*x_guess[k];
						inner_sum2 = inner_sum2+val2;}

					///////////////////////////////////////////////////////////
					//////////////// Accounting for Ghost Cells ///////////////
					///////////////////////////////////////////////////////////
					
					// Incomplete stencil in the negative y direction for nodes
					// along the second row (from the bottom) of the mesh. Starts
					// from first interior node along the second row and goes until
					// the last interior node of the row.
					if (i>N_x && i<2*N_x-1){
						double x_pos1 = (i-N_x)*delx;
						double y_pos1 = -1*delx;
						heat_ghost1 = masa_eval_2d_exact_t(x_pos1,y_pos1);}
					
					// Incomplete stencils in the positive y direction from nodes
					// that are on the row second from the top in the mesh. Starts
					// from the first interior node along second to last row of 
					// mesh and goes until the last interior node.
					if (i>vector_length-2*N_x && i<vector_length-N_x-1){
						double x_start = vector_length - 2*N_x;
						double x_pos1 = (i-x_start)*delx;
						double y_pos1 = N_x*delx;
						heat_ghost1 = masa_eval_2d_exact_t(x_pos1,y_pos1);}

					// Incomplete Stencils in the negative x directions from nodes
					// along the second column of mesh (left to right). Indexing 
					// the node just after the boundary node on the second column
					// of the mesh.
					if (i>N_x && i<vector_length-N_x && (i-1)%N_x==0){
						double x_pos2 = -1*delx;
						double y_pos2 = ((i-1)/N_x)*delx;
						heat_ghost2 = masa_eval_2d_exact_t(x_pos2,y_pos2);}

					// Incomplete Stencils in the posive x direction from nodes 
					// along the second to last column of the mesh (left to right).
					// Indexing to the node just before the boundary node at the 
					// righthandside of the mesh.
					if (i>N_x && i<vector_length-N_x && (i+2)%N_x==0){
						double x_pos2 = N_x*delx;
						double y_pos2 = (((i+2)/N_x)-1)*delx;
						heat_ghost2 = masa_eval_2d_exact_t(x_pos2,y_pos2);}

					// Evaluating inner sum for this iteration
					double val_tot = -1*inner_sum -inner_sum2 + heat_ghost1 + heat_ghost2 + sources[i];
					x_sol[i] = (1/A[i][i])*val_tot;


				}
				// Compute Loss
				for (int n=0;n<vector_length;n++){
					diff = diff+(pow(x_sol[n]-x_guess[n],2));}
				// Normalize Loss
				diff = sqrt(diff/vector_length);
				// Check error of loss
				if (diff>tol){
					for (int m=0;m<vector_length;m++){
						x_guess[m] = x_sol[m];}
				}
				counter+=1;
			}
			return x_sol;
		} // 2d-4th order end
	} 
}
