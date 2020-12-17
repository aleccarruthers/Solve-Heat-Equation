#include<iostream>
#include<vector>
#include<iomanip>
#include<grvy.h>

using namespace GRVY;
using std::cout;
using std::endl;
using std::vector;
using std::string;

void parse_inputs(int &dimension,int &finorder, int &max_iterations,std::string &solver, std::string &logger, std::string &nodes, std::string &verification, double &length, double &frequency, double &K){
	GRVY_Input_Class iparse;
	if (! iparse.Open("./input.dat")){exit(1);}
	if (iparse.Read_Var("nodes",&nodes)){
	cout << "Nodes: " << nodes << endl;}
	if (iparse.Read_Var("solver",&solver)){
	cout << "Numerical Solver: " << solver << endl;}
	if (iparse.Read_Var("dimension",&dimension)){cout << "Dimension: " << dimension << endl;}
	if (iparse.Read_Var("finorder",&finorder)){
	cout << "Finite Difference Order: " << finorder << endl;}
	if (iparse.Read_Var("max_iterations",&max_iterations)){
	cout << "Maximum Iterations: " << max_iterations << endl;}
	if (iparse.Read_Var("loglevel",&logger)){
	cout << "Output Mode: " << logger << endl;}
	if (iparse.Read_Var("length",&length)){
	cout << "Domain Length: " << length << endl;}
	if (iparse.Read_Var("frequency",&frequency)){
	cout << "MASA Frequency: " << frequency << endl;}
	if (iparse.Read_Var("k",&K)){
	cout << "Thermal Conductivity: " << K << endl;}
	if (iparse.Read_Var("verification",&verification)){
	cout << "Operating Mode: ";}

	// The validity of this input was checked here for convenience. All other 
	// inputs are checked in the following two functions 
	if (verification=="t"){cout << "verification" << endl;}
	else if (verification=="f"){cout << "standard" << endl;}
	else{cout << endl << "*** Not an accepted operating mode ***" << endl << "*** t or f ***" << endl;exit(1);}
}
