#include<iostream>
#include<cmath>
#include<vector>
using std::cout;
using std::endl;
using std::vector;

double l2Norm(int N, std::vector<double> num_solution, std::vector<double> real_solution){
	// Initialized accumulation error term
	double error = 0.0;
	// Compute L2 loss
	for (int i=0;i<N;i++){
		error = error + pow(num_solution[i] - real_solution[i],2);}
	// Normalize error
	error = sqrt(error/N);
	cout << " " << endl;
	return error;
}
