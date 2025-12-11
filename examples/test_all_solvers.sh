#!/bin/bash
echo "===== Testing all examples with solver types 1, 2, 3 ====="
echo ""

for solver in 1 2 3; do
    echo "===== Solver Type $solver ====="
    
    # Ex1
    echo "--- Ex1 (Alpha-Alpha Scattering) ---"
    cd Ex1
    # Modify solver_type in the code dynamically by using sed to replace the line
    sed -i.bak "s/solver_type = [0-9]/solver_type = $solver/" example1_hp.f90 2>/dev/null || \
    sed "s/solver_type = [0-9]/solver_type = $solver/" example1_hp.f90 > tmp && mv tmp example1_hp.f90
    gfortran -O3 -fopenmp -I../../src -o example1_hp example1_hp.f90 ../../libhprmat.a -L/opt/homebrew/opt/openblas/lib -lopenblas -fopenmp 2>/dev/null
    echo "2 30 1 6.0
4 1.0 3.0
0 0 0
-1 0 0 0" | ./example1_hp 2>&1 | grep -E "E \(MeV\)|Total time|Solver type"
    cd ..
    
    # Ex4 (largest matrices - 12 channels)
    echo "--- Ex4 (12C+alpha - 12 channels) ---"
    cd Ex4
    sed -i.bak "s/solver_type = [0-9]/solver_type = $solver/" example4_hp.f90 2>/dev/null || \
    sed "s/solver_type = [0-9]/solver_type = $solver/" example4_hp.f90 > tmp && mv tmp example4_hp.f90
    gfortran -O3 -fopenmp -I../../src -o example4_hp example4_hp.f90 ../../libhprmat.a -L/opt/homebrew/opt/openblas/lib -lopenblas -fopenmp 2>/dev/null
    echo "3 60 1 10.0
3 4.0 8.0
0 0 0
-1 0 0 0" | ./example4_hp 2>&1 | grep -E "E \(MeV\)|Total time|Solver type"
    cd ..
    
    echo ""
done

# Restore original solver_type
cd Ex1 && sed -i.bak "s/solver_type = [0-9]/solver_type = 1/" example1_hp.f90 2>/dev/null || \
sed "s/solver_type = [0-9]/solver_type = 1/" example1_hp.f90 > tmp && mv tmp example1_hp.f90; cd ..
cd Ex4 && sed -i.bak "s/solver_type = [0-9]/solver_type = 1/" example4_hp.f90 2>/dev/null || \
sed "s/solver_type = [0-9]/solver_type = 1/" example4_hp.f90 > tmp && mv tmp example4_hp.f90; cd ..

echo "===== All tests completed ====="
