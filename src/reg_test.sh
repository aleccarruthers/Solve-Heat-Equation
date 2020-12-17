#!/work/07043/ac1824/stampede2/bin/bats

# Check return code on running heatsolver without an input
@test "Verify heatsolver exits without input file" {
	run heatsolver
	[ "$status" -eq 1 ]
}

@test "Flag unsupported dimension for finite differencing" {
	sed -i "s/dimension = .*/dimension = 3/g" input.dat
	run heatsolver input.dat
	[ "$status" -eq "3" ]
}
	

@test "Flag unsupported finite differencing" {
	sed -i "s/finorder = .*/finorder = 8/g" input.dat
	sed -i "s/dimension = .*/dimension = 1/g" input.dat
	sed -i "s/solver = .*/solver = 'GS'/g" input.dat
	run heatsolver input.dat
	[ "$status" -eq "4" ]
}

@test "Flag unsupported numerical solver" {
	sed -i "s/solver = .*/solver = 'alec'/g" input.dat
	sed -i "s/finorder = .*/finorder = 2/g" input.dat
	run heatsolver input.dat
	[ "$status" -eq "5" ]
}

@test "Flag Jacobi 4th order finite difference input" {
	sed -i "s/solver = .*/solver = 'Jacobi'/g" input.dat
	sed -i "s/finorder = .*/finorder = 4/g" input.dat
	run heatsolver input.dat
	[ "$status" -eq "6" ]
}

# Check accuracy of 2nd order, 1D Jacobi (tol = 1e-9)
@test "Accuracy: Jacobi 2nd Order, 1D" {
	sed -i "s/dimension = .*/dimension = 1/g" input.dat
	sed -i "s/solver = .*/solver = 'Jacobi'/g" input.dat
	sed -i "s/nodes = .*/nodes = 60/g" input.dat
	sed -i "s/finorder = .*/finorder = 2/g" input.dat
	result=$(heatsolver input.dat | grep "L2 Norm: " | awk '{print $3}' | awk '{print $1-0.000592282813522}' | awk '{print sqrt($1*$1)}') 
	(( $(echo "$result < 0.000000001" |bc -l) )) 

}

# Check accuracy of 2nd order, 2D Jacobi (tol = 1e-9)
@test "Accuracy: Jacobi 2nd Order, 2D" {
	sed -i "s/dimension = .*/dimension = 2/g" input.dat
	sed -i "s/solver = .*/solver = 'Jacobi'/g" input.dat
	sed -i "s/nodes = .*/nodes = 100/g" input.dat
	sed -i "s/finorder = .*/finorder = 2/g" input.dat
	result=$(heatsolver input.dat | grep "L2 Norm: " | awk '{print $3}' | awk '{print $1-0.0109285412965}' | awk '{print sqrt($1*$1)}') 
	(( $(echo "$result < 0.000000001" |bc -l) )) 
}


# Check accuracy of 2nd order, 1D Gauss-Seidel (tol = 1e-9)
@test "Accuracy: Gauss-Seidel 2nd Order, 1D" {
	sed -i "s/dimension = .*/dimension = 1/g" input.dat
	sed -i "s/solver = .*/solver = 'GS'/g" input.dat
	sed -i "s/nodes = .*/nodes = 80/g" input.dat
	sed -i "s/finorder = .*/finorder = 2/g" input.dat
	result=$(heatsolver input.dat | grep "L2 Norm: " | awk '{print $3}' | awk '{print $1-0.000330999836401}' | awk '{print sqrt($1*$1)}') 
	(( $(echo "$result < 0.000000001" |bc -l) )) 
}


# Check accuracy of 2nd order, 2D Gauss-Seidel (tol = 1e-9)
@test "Accuracy: Gauss-Seidel 2nd Order, 2D" {
	sed -i "s/dimension = .*/dimension = 2/g" input.dat
	sed -i "s/solver = .*/solver = 'GS'/g" input.dat
	sed -i "s/nodes = .*/nodes = 400/g" input.dat
	sed -i "s/finorder = .*/finorder = 2/g" input.dat
	result=$(heatsolver input.dat | grep "L2 Norm: " | awk '{print $3}' | awk '{print $1-0.00253791012087}' | awk '{print sqrt($1*$1)}') 
	(( $(echo "$result < 0.000000001" |bc -l) )) 
}


# Check accuracy of 4th order, 1D Gauss-Seidel (tol = 1e-9)
@test "Accuracy: Gauss-Seidel 4th Order, 1D" {
	sed -i "s/dimension = .*/dimension = 1/g" input.dat
	sed -i "s/solver = .*/solver = 'GS'/g" input.dat
	sed -i "s/nodes = .*/nodes = 40/g" input.dat
	sed -i "s/finorder = .*/finorder = 4/g" input.dat
	result=$(heatsolver input.dat | grep "L2 Norm: " | awk '{print $3}' | awk '{print $1-0.00000293932998343}' | awk '{print sqrt($1*$1)}') 
	(( $(echo "$result < 0.000000001" |bc -l) )) 
}


# Check accuracy of 4th order, 2D Gauss-Seidel (tol = 1e-9)
@test "Accuracy: Gauss-Seidel 4th Order, 2D" {
	sed -i "s/dimension = .*/dimension = 2/g" input.dat
	sed -i "s/solver = .*/solver = 'GS'/g" input.dat
	sed -i "s/nodes = .*/nodes = 400/g" input.dat
	sed -i "s/finorder = .*/finorder = 4/g" input.dat
	result=$(heatsolver input.dat | grep "L2 Norm: " | awk '{print $3}' | awk '{print $1-0.00002300939869}' | awk '{print sqrt($1*$1)}') 
	(( $(echo "$result < 0.000000001" |bc -l) )) 
}
