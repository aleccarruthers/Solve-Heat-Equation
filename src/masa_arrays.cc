#include<iostream>
#include<vector>
#include "masa_arrays.h"
#include<iomanip>
#include<cmath>

using std::vector;
using std::cout;
using std::endl;

std::vector<double> masa_sources(int N, int dimension, int order,double length,double frequency,double K,std::string loglevel){

	// Initialize parameters
	std::vector<double> source(N,0.0);
	double freq = frequency;

	// Source coefficient that is dependent on the finite differencing order
	double coeff; 
	if (order==2){coeff = -1;}
	else {coeff = -12;}

	// X-axis length in 1D, X and Y-axis length in 2d
	double L = length;

	//////////////////////////// 1D case ////////////////////////////
	if (dimension==1){
		// Initialize Masa
		int model = masa_init("1D Heat Conduction","heateq_1d_steady_const");
		masa_init_param();
		// Set MASA parameters based off inputs from the input file
		masa_set_param("k_0",K);
		double K_0 = masa_get_param("k_0");
		masa_set_param("A_x",freq);
		// Spacing between nodes
		double delx = L/(N-1);
 		// Loop over all nodes and fill a source vector with sources
 		// and exact temperatures (exact temperatures from the boundary condition nodes
		for (int i=0;i<N;i++){
			double x = i*delx;
			// Start and end nodes will get an exact temperature from MASA
			if (i==0 || i==N-1){
				source[i]=masa_eval_1d_exact_t(x); // At ends use Masa exact sol
			}
			// Everywhere else, get the source term and multiply by coeff*delx^2.
			// The multiplication is just a result of the finite differencing expression
			else{ source[i]=coeff*pow(delx,2)*masa_eval_1d_source_t(x)/K_0;}
		}
		// Print the source vector if the loglevel is set to 'debug'
		if (loglevel=="debug"){
			cout << "1-D Source Vector with grid location:" << endl;
			cout << "Source Vector; Grid Position" << endl;
			for (int k=0;k<N;k++){
				double rounded_pos = k*delx;
				cout << std::setprecision(12) << source[k] << "; (";
				cout << std::setprecision(5) << rounded_pos << ")" << endl;}
		}
		return source;
	}
	//////////////////////////// 2D case ////////////////////////////
	else{
		// Initialize MASA
		int model = masa_init("2D Heat Conduction","heateq_2d_steady_const");
		masa_init_param();

		// The number of nodes in each row of the mesh
		double N_x = sqrt(N);

		// Set MASA parameters based off inputs from the input file
		masa_set_param("k_0",K);
		double K_0 = masa_get_param("k_0");
		masa_set_param("A_x",freq);
		masa_set_param("B_y",freq);

		// Spacing between nodes in both x and y direction because
		// a square mesh was assumed.
		double delx = L/(N_x-1);
		// Increment for total nodes in mesh
		int increment = 0;

		// Treating the bottom left node as the first node.
		// Loop over all columns of a row, then move up to the next row and repeat the process
		for (int row=0;row<N_x;row++) 
		{	
			// Current row in grid using actual distance, not the index
			double y = delx*row;
			for (int col=0;col<N_x;col++){
				// Current location along the x axis
				double x = delx*col; 
				// Fill in top and bottom boundary conditions first
				if (row==0 || row==N_x-1){
					source[increment] = masa_eval_2d_exact_t(x,y);
				}
				// Fill in boundary conditions for the sides
				else if (col==0 || col==N_x-1){
					source[increment] = masa_eval_2d_exact_t(x,y);
				}
				// Fill in heat source terms
				else{
					double source_t = masa_eval_2d_source_t(x,y);
					source[increment] = (coeff*source_t*(pow(delx,2)))/K_0;
				}
				increment+=1;
			}
		}
		// Print the source vector along with the spatial coordinates
		if (loglevel=="debug"){
			cout << "2-D Source Vector with grid locations" << endl;
			cout << "Source Vector; Grid Position" << endl;
			int row_count =0;
			int col_count =0;
			for (int g=0;g<N;g++){
				double x_pos = col_count*delx;
				double y_pos = row_count*delx;
				cout << std::setprecision(12) << source[g] << "; (";
				cout << std::setprecision(5) << x_pos << ", " << y_pos << ")" << endl;
				col_count+=1;
				if ((g+1)==N_x+row_count*N_x){
					col_count=0;
					row_count+=1;}
			
			}
				
		}
		return source;
	}
}

std::vector<double> masa_solution(int N, int dimension,double length,double frequency, double K,std::string loglevel){
	double freq = frequency;
	// Initialize parameters
	std::vector<double> solution(N,0.0);
	double L = length;
	// 1D case
	if (dimension==1){
		// Initialize MASA
		int model = masa_init("1D Heat Conduction","heateq_1d_steady_const");
		masa_init_param();
		// Set MASA parameters based off inputs from input file
		masa_set_param("k_0",K);
		masa_set_param("A_x",freq);
		// Spacing between nodes
		double delx = L/(N-1);
		for (int i=0;i<N;i++){ 
			double x = i*delx;
			// Get the exact MASA exact solution at each node position
			solution[i]=masa_eval_1d_exact_t(x); 
		}
		return solution;
	}
	// 2D case
	else{
		// Initialize MASA
		int model = masa_init("2D Heat Conduction","heateq_2d_steady_const");
		masa_init_param();
		// Set MASA parameters based off inputs from the input file
		masa_set_param("k_0",K);
		masa_set_param("A_x",freq);
		masa_set_param("B_y",freq);
		// Number of nodes along each axis
		double N_x = sqrt(N);
		// Spacing between each node in both the x and y direcetions
		double delx = L/(N_x-1);
		// Overall node index counter
		int increment = 0;
		// Considering the bottom left node to be the starting location.
		// Loop over all nodes in a row, before moving to the next row
		for (int row=0;row<N_x;row++)
		{	
			// Current height/row for the nodes using actual distance, not INDICES
			double y = delx*row;
			// Loop over all nodes in a row
			for (int col=0;col<N_x;col++){
				// x position of current node using spatial distance
				double x = delx*col; 
				solution[increment] = masa_eval_2d_exact_t(x,y);
				increment+=1;
			}
		}
		return solution;

	}
}
