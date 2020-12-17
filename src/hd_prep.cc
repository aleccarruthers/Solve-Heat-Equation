#include<iostream>
#include<vector>
#include<cmath>
using std::endl;
using std::cout;
using std::vector;

void hdf5_prep(std::vector<double> masa_solution, std::vector<double> numeric_sol, double *hd_masa, double *hd_sol,double error_est,double *hd_l2, int dimension, double *coord_x, double *coord_y, double length){
	// Convert Solution Vectors into Arrays for HDF5
	int vec_size = masa_solution.size();
	for (int i=0;i<vec_size;i++){
		hd_masa[i] = masa_solution[i];
		hd_sol[i] = numeric_sol[i];
	}
	// Create a 1 element array for L2 Norm
	hd_l2[0] = error_est;
	// Create coordinate array
	if (dimension==1){
		double spacing = length/(vec_size-1);
		for (int s=0;s<vec_size;s++){
			coord_x[s] = s*spacing;
		}
	}
	else if (dimension==2){
		float single_dim_nodes = sqrt(vec_size);
		int nodes_x = static_cast<int>(single_dim_nodes);
		double spacing2 = length/(nodes_x-1);
		int tot_count = 0;
		for (int row=0;row<nodes_x;row++){
			for (int col=0;col<nodes_x;col++){
				double dist_y = row*spacing2;
				double dist_x = col*spacing2;
				coord_x[tot_count] = dist_x;
				coord_y[tot_count] = dist_y;
				tot_count += 1;
			}
		}
	}
}
