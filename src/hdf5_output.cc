#include<iostream>
#include "hdf5.h"
#include<vector>

using std::cout;
using std::vector;
using std::endl;

void hdf5_write(double *num_sol, double *ana_sol,double *l2, int N, double *coord_x, double *coord_y, int dimension){
	// Create dataset for the numerical solution
	hid_t           file, dataset, dataspace;
	hsize_t         dimsf[1];
	herr_t          status;
	//size_t		precision, offset;
	file = H5Fcreate("data.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	dimsf[0] = N;
	dataspace = H5Screate_simple(1,dimsf,NULL);
	dataset = H5Dcreate2(file, "Numerical Solution", H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, num_sol);
	H5Sclose(dataspace);
	H5Dclose(dataset);

	// Create dataset for analytic solution
	hid_t           dataset2, dataspace2;
	hsize_t         dimsf2[1];
	dimsf2[0] = N;
	dataspace2 = H5Screate_simple(1,dimsf2,NULL);
	dataset2 = H5Dcreate2(file, "Analytic Solution", H5T_NATIVE_DOUBLE, dataspace2,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ana_sol);
	H5Sclose(dataspace2);
	H5Dclose(dataset2);

	// Create Dataset for L2 Norm (This will be used for regression tests)
	hid_t           dataset3, dataspace3;
	hsize_t         dimsf3[1];
	dimsf3[0] = 1;
	dataspace3 = H5Screate_simple(1,dimsf3,NULL);
	dataset3 = H5Dcreate2(file, "L2 Loss", H5T_NATIVE_DOUBLE, dataspace3,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, l2);
	H5Sclose(dataspace3);
	H5Dclose(dataset3);

	// X position
	hid_t           dataset4, dataspace4;
	hsize_t         dimsf4[1];
	dimsf4[0] = N;
	dataspace4 = H5Screate_simple(1,dimsf4,NULL);
	dataset4 = H5Dcreate2(file, "X Coordinates", H5T_NATIVE_DOUBLE, dataspace4,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset4, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord_x);
	H5Sclose(dataspace4);
	H5Dclose(dataset4);

	if (dimension==2){
		
		hid_t           dataset5, dataspace5;
		hsize_t         dimsf5[1];
		dimsf5[0] = N;
		dataspace5 = H5Screate_simple(1,dimsf5,NULL);
		dataset5 = H5Dcreate2(file, "Y Coordinates", H5T_NATIVE_DOUBLE, dataspace5,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset5, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord_y);
		H5Sclose(dataspace5);
		H5Dclose(dataset5);
		
	}

	H5Fclose(file);
}
