#!/usr/bin/env bats
cd ../src/

# Check Accuracy of 1D-2nd Order
@test "Accuracy: PETSc 2nd Order, 1D" {
	sed -i "s/dimension = .*/dimension = 1/g" ../src/input.dat
	sed -i "s/solver = .*/solver = petsc/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 2/g" ../src/input.dat
	sed -i "s/loglevel = .*/loglevel = standard/g" ../src/input.dat
	sed -i "s/nodes = .*/nodes = 60/g" ../src/input.dat
	result=$(../src/heatsolver ../src/input.dat)
	difference=$(h5diff -v --delta=1e-12 ../src/data.h5 ../tests/hdf_tests/pet1dTwoA.h5 /"L2 Loss"/ /"L2 Loss"/ | grep -P "differences" | awk '{print $1}')
	[ "$difference" -eq "0" ]

}

# Check Accuracy of 1D-4th Order
@test "Accuracy: PETSc 4th Order, 1D" {
	sed -i "s/finorder = .*/finorder = 4/g" ../src/input.dat
	sed -i "s/nodes = .*/nodes = 20/g" ../src/input.dat
	result=$(../src/heatsolver ../src/input.dat)
	difference=$(h5diff -v --delta=1e-12 ../src/data.h5 ../tests/hdf_tests/pet1dFourA.h5 /"L2 Loss"/ /"L2 Loss"/ | grep -P "differences" | awk '{print $1}')
	[ "$difference" -eq "0" ]
}

# Check Accuracy of 2D-2nd Order
@test "Accuracy: PETSc 2nd Order, 2D" {
	sed -i "s/dimension = .*/dimension = 2/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 2/g" ../src/input.dat
	sed -i "s/nodes = .*/nodes = 10/g" ../src/input.dat
	result=$(../src/heatsolver ../src/input.dat)
	difference=$(h5diff -v --delta=1e-12 ../src/data.h5 ../tests/hdf_tests/pet2dTwoA.h5 /"L2 Loss"/ /"L2 Loss"/ | grep -P "differences" | awk '{print $1}')
	[ "$difference" -eq "0" ]
}

# Check Accuracy of 2D-4th Order
@test "Accuracy: PETSc 4th Order, 2D" {
	sed -i "s/dimension = .*/dimension = 2/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 4/g" ../src/input.dat
	sed -i "s/nodes = .*/nodes = 10/g" ../src/input.dat
	result=$(../src/heatsolver ../src/input.dat)
	difference=$(h5diff -v --delta=1e-12 ../src/data.h5 ../tests/hdf_tests/pet2dFourA.h5 /"L2 Loss"/ /"L2 Loss"/ | grep -P "differences" | awk '{print $1}')
	[ "$difference" -eq "0" ]

}
