#!/usr/bin/env bats
cd ../src/
# Check return code on running ../src/heatsolver without an ../src/input
@test "Verify ../src/heatsolver exits without input file" {
	run ../src/heatsolver
	[ "$status" -eq 1 ]
}

@test "Verfiy heatsolver exit on an unsupported domain dimension" {
	sed -i "s/dimension = .*/dimension = 3/g" ../src/input.dat
	run ../src/heatsolver ../src/input.dat
	[ "$status" -eq 1 ]
}
	

@test "Verify heatsolver exit on unsupported finite differencing order" {
	sed -i "s/finorder = .*/finorder = 8/g" ../src/input.dat
	sed -i "s/dimension = .*/dimension = 1/g" ../src/input.dat
	sed -i "s/solver = .*/solver = 'GS'/g" ../src/input.dat
	run ../src/heatsolver ../src/input.dat
	[ "$status" -eq "1" ]
}

@test "Verify exit on unsupported numerical solver method" {
	sed -i "s/solver = .*/solver = 'alec'/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 2/g" ../src/input.dat
	run ../src/heatsolver ../src/input.dat
	[ "$status" -eq "1" ]
}

@test "Verify exit on the numerically unstable Jacobi 4th order finite difference" {
	sed -i "s/solver = .*/solver = 'Jacobi'/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 4/g" ../src/input.dat
	run ../src/heatsolver ../src/input.dat
	[ "$status" -eq "1" ]
}


@test "Verify exit on invalid output mode" {
	sed -i "s/solver = .*/solver = 'GS'/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 2/g" ../src/input.dat
	sed -i "s/loglevel = .*/loglevel = foo/g" ../src/input.dat
	run ../src/heatsolver ../src/input.dat
	[ "$status" -eq "1" ]
}
@test "Verify exit on a number of nodes input containing special characters" {
	sed -i "s/solver = .*/solver = 'Jacobi'/g" ../src/input.dat
	sed -i "s/loglevel = .*/loglevel = standard/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 2/g" ../src/input.dat
	sed -i "s/nodes = .*/nodes = 100,000/g" ../src/input.dat
	run ../src/heatsolver ../src/input.dat
	[ "$status" -eq "1" ]
}

@test "Verify exit on floating point input for number of nodes" {
	sed -i "s/nodes = .*/nodes = 100.12/g" ../src/input.dat
	run ../src/heatsolver ../src/input.dat
	[ "$status" -eq "1" ]
}

@test "Verify exit on non-positive input for number of nodes" {
	sed -i "s/nodes = .*/nodes = -25/g" ../src/input.dat
	run ../src/heatsolver ../src/input.dat
	[ "$status" -eq "1" ]
}

@test "Checking exit on negative or zero length" {
	sed -i "s/length = .*/length = -.1/g" ../src/input.dat
	run ../src/heatsolver ../src/input.dat
	[ "$status" -eq "1" ]
}

@test "Checking exit on non-positive MASA frequency" {
	sed -i "s/length = .*/length = 1.0/g" ../src/input.dat
	sed -i "s/frequency = .*/frequency = -1.0/g" ../src/input.dat
	run ../src/heatsolver ../src/input.dat
	[ "$status" -eq "1" ]
}

@test "Verify Numerical Stopping at Max Iterations" {
	sed -i "s/frequency = .*/frequency = 5.0/g" ../src/input.dat
	sed -i "s/max_iterations = .*/max_iterations = 5/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 2/g" ../src/input.dat
	run ../src/heatsolver ../src/input.dat
	[ "$status" -eq "1" ]
}

# Note: I am checking whetehr h5diff says there is a difference or not, not the actual values
# Check accuracy of 2nd order, 1D Jacobi (tol = 1e-9)
@test "Accuracy: Jacobi 2nd Order, 1D" {
	sed -i "s/dimension = .*/dimension = 1/g" ../src/input.dat
	sed -i "s/solver = .*/solver = 'Jacobi'/g" ../src/input.dat
	sed -i "s/nodes = .*/nodes = 60/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 2/g" ../src/input.dat
	sed -i "s/max_iterations = .*/max_iterations = 250000/g" ../src/input.dat
	#result=$(../src/heatsolver ../src/input.dat | grep "L2 Norm: " | awk '{print $3}' | awk '{print $1-0.000592282813522}' | awk '{print sqrt($1*$1)}') 
	#h5diff ../src/data.h5 ../tests/hdf_tests/data_jac1d.h5
	result=$(../src/heatsolver ../src/input.dat)
	difference=$(h5diff -v --delta=1e-12 ../src/data.h5 ../tests/hdf_tests/data_jac1d.h5 /"L2 Loss"/ /"L2 Loss"/ | grep -P "differences" | awk '{print $1}')
	[ "$difference" -eq "0" ]
}

# Check accuracy of 2nd order, 2D Jacobi (tol = 1e-9)
@test "Accuracy: Jacobi 2nd Order, 2D" {
	sed -i "s/dimension = .*/dimension = 2/g" ../src/input.dat
	sed -i "s/solver = .*/solver = 'Jacobi'/g" ../src/input.dat
	sed -i "s/nodes = .*/nodes = 100/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 2/g" ../src/input.dat
	#result=$(../src/heatsolver ../src/input.dat | grep "L2 Norm: " | awk '{print $3}' | awk '{print $1-0.0109285412965}' | awk '{print sqrt($1*$1)}') 
	result=$(../src/heatsolver ../src/input.dat)
	difference=$(h5diff -v --delta=1e-12 ../src/data.h5 ../tests/hdf_tests/data_jac2d.h5 /"L2 Loss"/ /"L2 Loss"/ | grep -P "differences" | awk '{print $1}')
	#(( $(echo "$result < 0.000000001" |bc -l) ))
	[ "$difference" -eq "0" ]
}


# Check accuracy of 2nd order, 1D Gauss-Seidel (tol = 1e-9)
@test "Accuracy: Gauss-Seidel 2nd Order, 1D" {
	sed -i "s/dimension = .*/dimension = 1/g" ../src/input.dat
	sed -i "s/solver = .*/solver = 'GS'/g" ../src/input.dat
	sed -i "s/nodes = .*/nodes = 80/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 2/g" ../src/input.dat
	result=$(../src/heatsolver ../src/input.dat)
	difference=$(h5diff -v --delta=1e-12 ../src/data.h5 ../tests/hdf_tests/data_gs1d_2f.h5 /"L2 Loss"/ /"L2 Loss"/ | grep -P "differences" | awk '{print $1}') 
	[ "$difference" -eq "0" ]
}


# Check accuracy of 2nd order, 2D Gauss-Seidel (tol = 1e-9)
@test "Accuracy: Gauss-Seidel 2nd Order, 2D" {
	sed -i "s/dimension = .*/dimension = 2/g" ../src/input.dat
	sed -i "s/solver = .*/solver = 'GS'/g" ../src/input.dat
	sed -i "s/nodes = .*/nodes = 400/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 2/g" ../src/input.dat
	result=$(../src/heatsolver ../src/input.dat)
	difference=$(h5diff -v --delta=1e-12 ../src/data.h5 ../tests/hdf_tests/data_gs2d_2f.h5 /"L2 Loss"/ /"L2 Loss"/ | grep -P "differences" | awk '{print $1}')
	[ "$difference" -eq "0" ]
}


# Check accuracy of 4th order, 1D Gauss-Seidel (tol = 1e-9)
@test "Accuracy: Gauss-Seidel 4th Order, 1D" {
	sed -i "s/dimension = .*/dimension = 1/g" ../src/input.dat
	sed -i "s/solver = .*/solver = 'GS'/g" ../src/input.dat
	sed -i "s/nodes = .*/nodes = 40/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 4/g" ../src/input.dat
	result=$(../src/heatsolver ../src/input.dat)
	difference=$(h5diff -v --delta=1e-12 ../src/data.h5 ../tests/hdf_tests/data_gs1d_4f.h5 /"L2 Loss"/ /"L2 Loss"/ | grep -P "differences" | awk '{print $1}') 
	[ "$difference" -eq "0" ]
}


# Check accuracy of 4th order, 2D Gauss-Seidel (tol = 1e-9)
@test "Accuracy: Gauss-Seidel 4th Order, 2D" {
	sed -i "s/dimension = .*/dimension = 2/g" ../src/input.dat
	sed -i "s/solver = .*/solver = 'GS'/g" ../src/input.dat
	sed -i "s/nodes = .*/nodes = 400/g" ../src/input.dat
	sed -i "s/finorder = .*/finorder = 4/g" ../src/input.dat
	result=$(../src/heatsolver ../src/input.dat)
	difference=$(h5diff -v --delta=1e-12 ../src/data.h5 ../tests/hdf_tests/data_gs2d_4f.h5 /"L2 Loss"/ /"L2 Loss"/ | grep -P "differences" | awk '{print $1}') 
	[ "$difference" -eq "0" ] 
}


# Check DEBUG mode output
@test "Check Debug Mode Output" {
	sed -i "s/nodes = .*/nodes = 25/g" ../src/input.dat
	sed -i "s/loglevel = .*/loglevel = debug/g" ../src/input.dat
	result=$(../src/heatsolver ../src/input.dat | grep "A:") 
}
