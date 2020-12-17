#include<iostream>
#include<cmath>
#include<string>
using std::cout;
using std::endl;
using std::string;
int node_check(int dimension, std::string nodes){
	int nodes_final;
	char* p;
	double n_val = strtod(nodes.c_str(),&p);

	// Check for special characters in the input number of nodes
	if (*p){
		cout << "*** Special Characters found in the input number of nodes ***" << endl;
		cout << "*** The input number of nodes must be a real, positive integer ***" << endl;
		exit(1);}

	// Check for floating point input for number of nodes
	double double_check = n_val - floor(n_val);
	if (double_check!=0){
		cout << "*** Floating point input for the number of nodes detected ***" << endl;
		cout << "*** The input number of nodes must be a positive integer ***" << endl;
		exit(1);}

	// Check for non-positive number of nodes
	double abs_check = n_val + abs(n_val);
	if (abs_check==0){
		cout << "*** The input number of nodes must be a positive integer ***" << endl;
		exit(1);}

	// If the previous three tests are passed, the input number of n_nodes can be used
	int n_nodes = n_val;
	
	// For the 1D case, 'n_nodes' is the total number of nodes used

	// Check to see if the input number of nodes has a perfect square
	if (dimension==2){
		double node_sq = sqrt(n_nodes);
		// Two methods for setting the number of nodes:
		// (1): If 'n_nodes' has a perfect square, i.e. 25, use 
		// 'n_nodes as the total number of nodes. The
		// domain would have size: sqrt(n_nodes) x sqrt(n_nodes).

		// (2): If 'n_nodes' does not a have a perfect square, i.e. 17,
		// treat 'n_nodes' as the number of nodes
		// along each axis. The domain would have size: n_nodes x n_nodes.

		// Only numbers with a perfect square will pass this condition 
		if ((node_sq - floor(node_sq)) == 0){
			nodes_final = n_nodes;}
		else{
			nodes_final = pow(n_nodes,2);
			cout << " " << endl;
			cout << "Input 'nodes' does not have a perfect square" << endl;
			cout << "Using input 'nodes' as number of nodes in x and y direction" << endl;
			cout << "Domain Size: " << n_nodes << " x " << n_nodes << endl;
			cout << " " << endl;
		}
	}
	else{
		nodes_final = n_nodes;}
	return nodes_final;
}
