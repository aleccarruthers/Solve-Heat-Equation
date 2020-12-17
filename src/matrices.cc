#include<iostream>
#include "matrices.h"
#include<vector>
#include<cmath>

using std::vector;
using std::cout;
using std::endl;

// Creating the A matrix for the system of equations based on the input number of dimensions
// and finite differencing order
std::vector<std::vector<double>> A_matrix(int N, int dimension, int fin_order,std::string loglevel){
	std::vector<std::vector<double>> A_mat;
	// 1-D
	if (dimension==1){
		if (fin_order==2){
			A_mat = D1_2nd(N);
		}
		else if (fin_order==4){
			A_mat = D1_4th(N);
		}
	}
	// 2-D
	else if (dimension==2){
		if (fin_order==2){
			A_mat = D2_2nd(N);
		}
		else if (fin_order==4){
			A_mat = D2_4th(N);
		}
	}
	// Throw an error if the matrix was not constructed to the appropriate number of dimensions
	if (A_mat.size()!=N){
		cout << "*** Allocation failure for the specified number of nodes. Possible causes: ***" << endl;
		cout << "(1) The input number of nodes has been corrupted or modified" << endl;
		cout << "(2) Not enough memory to store a matrix of the requested size" << endl;
		fprintf( stderr, "*** Failed memory allocation for the A matrix in the linear system ***\n");
		exit(1);
	}

	// Check for whether standard or debug output mode was set
	if (loglevel=="debug"){
		cout << " " << endl;
		cout << "DEBUG MODE OUTPUT: " << endl;
		cout << "A: " << endl;
		print_matrix(A_mat,N);}

	return A_mat;
}

////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// 1-Dimensional //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
// 2nd Order
std::vector<std::vector<double>> D1_2nd(int N){
	// Calculate spacing between nodes
	float del_x = N/(N-1);
	std::vector<std::vector<double>> node_vec(N, std::vector<double>(N, 0.0));
	for (int row_c=0;row_c<N; row_c++){
		int ind = row_c-1;
		// Taking care of start and end boundary conditions before looping over interior nodes
		if (row_c==0){node_vec[row_c][0]=1;}
		else if (row_c==N-1){node_vec[row_c][N-1]=1;}
		
		// Loop over all nodes with a complete stencil
		else{
			for (int count=0;count<3;count++){
				if (count==0 || count==2){
					node_vec[row_c][ind]=1.0;}
				else if (count==1){
					node_vec[row_c][ind]=-2.0;
				}
				ind = ind+1;
			}
		}
	}
	return node_vec;	
}	

// 4th order
std::vector<std::vector<double>> D1_4th(int N){
	// Calculate node spacing
	float del_x = N/(N-1);
	std::vector<std::vector<double>> node_vec(N, std::vector<double>(N, 0.0));
	for (int row_c=0;row_c<N; row_c++){
		// Index within each row
		int ind = row_c-1;

		// Resolve nodes on boundary before iterating over interior nodes
		if (row_c==0){node_vec[row_c][0]=1;}
		else if (row_c==N-1){node_vec[row_c][N-1]=1;}
		
		///////////////////// Nodes with incomplete stencils /////////////////////
		// Second Node from start
		else if (row_c==1){
			for (int col_num=0;col_num<4;col_num++)
			{
				if (col_num==0){node_vec[row_c][ind]=16;}
				else if (col_num==1){node_vec[row_c][ind]=-30;}
				else if (col_num==2){node_vec[row_c][ind]=16;}
				else {node_vec[row_c][ind]=-1;}
				ind = ind + 1;
			}

		}
		// Second to last node
		else if (row_c==N-2){
			for (int col_num=0;col_num<4;col_num++)
			{
				if (col_num==0){node_vec[row_c][ind-1]=-1;}
				else if (col_num==1){node_vec[row_c][ind-1]=16;}
				else if (col_num==2){node_vec[row_c][ind-1]=-30;}
				else {node_vec[row_c][ind-1]=16;}
				ind = ind + 1;
			}
		}
		///////////////////// End Nodes with Incomplete Stencil /////////////////////

		// Loop over and fill "A" matrix for all interior nodes with a complete stencil
		else{
			for (int col_num=0;col_num<5;col_num++){
				// First and last components of stencil
				if (col_num==0 || col_num==4){
					node_vec[row_c][ind-1]=-1.0;}

				// Second and second to last components of stencil
				else if (col_num==1 || col_num==3){
					node_vec[row_c][ind-1]=16.0;
				}
				
				// Interior node of stencil
				else {node_vec[row_c][ind-1]=-30;}
				ind = ind+1;
			}
		}
	}
	return node_vec;	
}	

////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// 2-Dimensional //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
// 2nd Order
std::vector<std::vector<double>> D2_2nd(int N){
	// Number of nodes in the x direction (also # of nodes in y direction)
	int N_x = sqrt(N);
	// All nodes above this value are at the top boundary of the domain
	int threshold = N - N_x; 
	// Spacing between nodes
	float del_x = N_x/(N_x-1); 
	std::vector<std::vector<double>> node_vec(N, std::vector<double>(N, 0.0));
	
	// Filling values for top and bottom boundary conditions
	for (int ind=0;ind<N; ind++){
		if (ind<N_x || ind>=threshold){node_vec[ind][ind]=1;}
	}
	// Filling values for left and right boundary conditions
	for (int bc=N_x;bc<threshold;bc+=N_x){
		node_vec[bc][bc]=1;
		node_vec[bc+N_x-1][bc+N_x-1]=1;
	}
	
	// Filling nodes with complete stencils
	int node_int = N_x + 1; //  Location of first complete stencil
	int non_bc = N_x - 2; // Number of interior nodes in a row before boundary
	for (int vals=node_int; vals<threshold; vals+=N_x){
		for (int g=0;g<non_bc;g++){
			node_vec[vals+g][vals+g] = -4;
			node_vec[vals+g][vals-1+g] = 1;
			node_vec[vals+g][vals+1+g] = 1;
			node_vec[vals+g][vals+N_x+g] = 1;
			node_vec[vals+g][vals-N_x+g] = 1;
		}
	}
	return node_vec;

}

// 4th Order
std::vector<std::vector<double>> D2_4th(int N){
	// Number of nodes along the x-axis (also the # of nodes along y-axis)
	int N_x = sqrt(N);
	// All nodes above this value are at the top boundary
	int threshold = N - N_x;
	// Spacing between nodes 
	float del_x = N_x/(N_x-1); 
	std::vector<std::vector<double>> node_vec(N, std::vector<double>(N, 0.0));
	
	// Filling values for top and bottom boundary conditions
	for (int ind=0;ind<N; ind++){
		if (ind<N_x || ind>=threshold){node_vec[ind][ind]=1;}
	}
	// Filling values for left and right boundary conditions
	for (int bc=N_x;bc<threshold;bc+=N_x){
		node_vec[bc][bc]=1;
		node_vec[bc+N_x-1][bc+N_x-1]=1;
	}
	
	/////// Filling in the second and second to last row of nodes (Incomplete stencils) ///////
	// bc2 defines the first interior node with an incomplete stencil
	int bc2 = N_x+1;
	for (int n=0;n<N_x-2;n++){
		// center node of stencil
		node_vec[n+bc2][n+bc2]=-60;
		// nodes above and below center node in stencil
		node_vec[n+bc2][n+bc2+N_x]=16;
		node_vec[n+bc2][n+bc2-N_x]=16;
		// All nodes along the second row of the mesh have a second node in positive y direction
		node_vec[n+bc2][n+bc2+(2*N_x)]=-1;
		
		// Leftmost interior node has an incomplete stencil in both the
		// negative x and y directions
		if (n==0){
			node_vec[n+bc2][n+bc2+1]=16;
			node_vec[n+bc2][n+bc2-1]=16;
			node_vec[n+bc2][n+bc2+2]=-1;
		}
		// Rightmost interior node with incomplete stencil in the 
		// positive x and negative y directions
		else if (n==N_x-3){
			node_vec[n+bc2][n+bc2+1]=16;
			node_vec[n+bc2][n+bc2-1]=16;
			node_vec[n+bc2][n+bc2-2]=-1;
		}
		// Nodes along the row with an incomplete stencil from only the negative y direction
		else{
			node_vec[n+bc2][n+bc2+1]=16;
			node_vec[n+bc2][n+bc2-1]=16;
			node_vec[n+bc2][n+bc2-2]=-1;
			node_vec[n+bc2][n+bc2+2]=-1;
		}

	}
		
	// bcU defines the first interior node along the second highest row of nodes in the mesh
	int bcU = N-(2*N_x)+1;
	for (int n=0;n<N_x-2;n++){
		// center node
		node_vec[n+bcU][n+bcU]=-60;
		// nodes immediately above and below center node
		node_vec[n+bcU][n+bcU+N_x]=16;
		node_vec[n+bcU][n+bcU-N_x]=16;
		// second node in stencil along negative y direction
		node_vec[n+bcU][n+bcU-(2*N_x)]=-1;
		
		// first interior node has two nodes missing. One in the positive y
		// and the other in the negative x
		if (n==0){
			node_vec[n+bcU][n+bcU+1]=16;
			node_vec[n+bcU][n+bcU-1]=16;
			node_vec[n+bcU][n+bcU+2]=-1;
		}
		// The last interior node in this row has missing nodes in the
		// positive y and positive x directions
		else if (n==N_x-3){
			node_vec[n+bcU][n+bcU+1]=16;
			node_vec[n+bcU][n+bcU-1]=16;
			node_vec[n+bcU][n+bcU-2]=-1;
		}
		// All nodes in this row with just a missing node in the positive y direction
		else{
			node_vec[n+bcU][n+bcU+1]=16;
			node_vec[n+bcU][n+bcU-1]=16;
			node_vec[n+bcU][n+bcU-2]=-1;
			node_vec[n+bcU][n+bcU+2]=-1;
		}

	}


	/////// Filling interior rows ///////
	// First non boundary condition node along third row of grid (third from the bottom)
	int node_int = 2*N_x + 1;
	// Index of the last interior node with a complete stencil
	int node_last = N - (3*N_x)+2;
	// Number of nodes in a row on the grid before a boundary
	int non_bc = N_x - 2;
	// Start from the third row and move across all nodes in the row, then increase
	// the index by N_x, which moves the calculations to the next row's starting node.
	for (int vals=node_int; vals<node_last; vals+=N_x){
		for (int g=0;g<non_bc;g++){
			node_vec[vals+g][vals+g] = -60;
			// All nodes in these cases will have complete stencils in the
			// y direction
			node_vec[vals+g][vals+N_x+g] = 16;
			node_vec[vals+g][vals-N_x+g] = 16;
			node_vec[vals+g][vals+(2*N_x)+g] = -1;
			node_vec[vals+g][vals-(2*N_x)+g] = -1;
			// interior node with incomplete stencil in negative x direction
			if (g==0){	
				node_vec[vals+g][vals-1+g] = 16;
				node_vec[vals+g][vals+1+g] = 16;
				node_vec[vals+g][vals+2+g] = -1;
			}
			// second to last interior node in the row and it has an incomplete stencil
			// in the positive x direction
			else if (g==non_bc-1){
				node_vec[vals+g][vals-1+g] = 16;
				node_vec[vals+g][vals+1+g] = 16;
				node_vec[vals+g][vals-2+g] = -1;

			}
			// interior nodes with a complete stencil
			else{
				node_vec[vals+g][vals-1+g] = 16;
				node_vec[vals+g][vals+1+g] = 16;
				node_vec[vals+g][vals-2+g] = -1;
				node_vec[vals+g][vals+2+g] = -1;

			}
		}
	}
	return node_vec;
}
// function to print the resulting A matrix
void print_matrix(std::vector<std::vector<double>> matrix, int n){
	for (int i=0;i<n; i++){
		for (int j = 0;j<n;j++){
			cout << matrix[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;
}	
