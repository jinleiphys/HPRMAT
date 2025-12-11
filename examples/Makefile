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
all: Ex0 Ex1

Ex0: Ex0/test_solver
Ex1: Ex1/example1_hp

Ex0/test_solver: Ex0/test_solver.f90 $(LIBHPRMAT)
	$(FC) $(FFLAGS) -I../src -o $@ $< $(LIBHPRMAT) $(LIBS)
	@echo "Built Ex0/test_solver"

Ex1/example1_hp: Ex1/example1_hp.f90 $(LIBHPRMAT)
	$(FC) $(FFLAGS) -I../src -o $@ $< $(LIBHPRMAT) $(LIBS)
	@echo "Built Ex1/example1_hp"

clean:
	rm -f Ex0/test_solver Ex1/example1_hp *.o *.mod

.PHONY: all Ex0 Ex1 clean
