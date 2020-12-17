#include<iostream>

using std::cout;
using std::endl;
using std::string;

int petsc_error(){
	cout << "*** PETSc was not enabled during configuration ***" << endl;
	cout << "To use PETSc: " << endl;
	cout << "module load petsc" << endl;
	cout << "./configure ......... --with-petsc=$TACC_PETSC_DIR" << endl;
	cout << "Make" << endl;
	return 1;
}
