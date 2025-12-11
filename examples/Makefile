#==============================================================================
# HPRMAT Examples Makefile
#==============================================================================

# Include parent configuration
-include ../make.inc

# Default compilers
FC ?= gfortran
FFLAGS ?= -O3 -fopenmp

# Library path
LIBHPRMAT = ../libhprmat.a

# Default BLAS (can be overridden)
LIBS ?= -llapack -lblas -fopenmp

# Targets
all: Ex0 Ex1 Ex2 Ex3 Ex4 Ex5

Ex0: Ex0/test_solver
Ex1: Ex1/example1_hp
Ex2: Ex2/example2_hp
Ex3: Ex3/example3_hp
Ex4: Ex4/example4_hp
Ex5: Ex5/example5_hp

Ex0/test_solver: Ex0/test_solver.f90 $(LIBHPRMAT)
	$(FC) $(FFLAGS) -I../src -o $@ $< $(LIBHPRMAT) $(LIBS)
	@echo "Built Ex0/test_solver"

Ex1/example1_hp: Ex1/example1_hp.f90 $(LIBHPRMAT)
	$(FC) $(FFLAGS) -I../src -o $@ $< $(LIBHPRMAT) $(LIBS)
	@echo "Built Ex1/example1_hp"

Ex2/example2_hp: Ex2/example2_hp.f90 $(LIBHPRMAT)
	$(FC) $(FFLAGS) -I../src -o $@ $< $(LIBHPRMAT) $(LIBS)
	@echo "Built Ex2/example2_hp"

Ex3/example3_hp: Ex3/example3_hp.f90 $(LIBHPRMAT)
	$(FC) $(FFLAGS) -I../src -o $@ $< $(LIBHPRMAT) $(LIBS)
	@echo "Built Ex3/example3_hp"

Ex4/example4_hp: Ex4/example4_hp.f90 $(LIBHPRMAT)
	$(FC) $(FFLAGS) -I../src -o $@ $< $(LIBHPRMAT) $(LIBS)
	@echo "Built Ex4/example4_hp"

Ex5/example5_hp: Ex5/example5_hp.f90 $(LIBHPRMAT)
	$(FC) $(FFLAGS) -I../src -o $@ $< $(LIBHPRMAT) $(LIBS)
	@echo "Built Ex5/example5_hp"

clean:
	rm -f Ex0/test_solver Ex1/example1_hp Ex2/example2_hp Ex3/example3_hp Ex4/example4_hp Ex5/example5_hp *.o *.mod

.PHONY: all Ex0 Ex1 Ex2 Ex3 Ex4 Ex5 clean
