#include<iostream>
using std::cout;
using std::endl;

int params_check(int dimension, int order,std::string solver,float length,float freq,std::string loglevel,int max_iterations){
	int return_val = 0;
	// Checking the input dimension for the mesh
	if (dimension!=1 && dimension !=2){
		cout << "*** Supported dimensions are 1 or 2 ***" << endl;
		fprintf( stderr, "*** Not a supported dimension ***\n");
		return_val = 1;}
	// Check Max Iterations
	if (max_iterations <= 0){
		cout << "*** Max Iterations must be a positive integer ***" << endl;
		fprintf( stderr, "*** Max Iterations must be a positive integer ***\n");
		return_val = 1;}
	// Checking the input finite differencing order
	if (order !=2 && order !=4){
		cout << "*** Not a supported finite differencing order ***" << endl;
		cout << "*** Available differencing orders: 2 or 4 ***" << endl;
		//fprintf( stderr, "*** Not a supported finite differencing order ***\n");
		return_val = 1;}

	// Checking the input numerical solver method
	if (solver!="Jacobi" && solver != "GS" && solver != "petsc"){
		cout << "*** Supported Iterative Solvers: 'Jacobi', 'GS' (Gauss-Seidel), and 'petsc' ***" << endl;
		fprintf( stderr, "*** Not a supported numerical solver ***\n");
		return_val = 1;}

	// Checking the output mode
	if (loglevel!="standard" && loglevel != "debug"){
		cout << "*** Not a supported loglevel. Supported loglevels: 'standard' and 'debug' ***" << endl;
		fprintf( stderr, "*** Not a loglevel ***\n");
		return_val = 1;}
	// Checking if the numerically unstable Jacobi 4th order was input
	if (solver=="Jacobi" && order ==4){
		cout << "*** Jacobi 4th order finite differencing is not supported ***" << endl;
		fprintf( stderr, "*** 4th order Jacobi numerical solver is unstable ***\n");
		return_val = 1;}

	// Checking the input length
	if (length<=0){
		cout << "*** Length must be a positive value"<< endl;
		fprintf( stderr, "*** Input length is less than or equal to zero ***\n");
		return_val = 1;}
	
	// Checking the input frequency for MASA
	if (freq<=0){
		cout << "*** MASA frequency must be a positive value ***" << endl;
		fprintf( stderr, "*** Input MASA frequency is less than or equal to zero ***\n");
		return_val = 1;}

	if (return_val==0){
		return 0;}
	else{
		exit(return_val);}
}
