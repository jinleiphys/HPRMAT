# HPRMAT Python Binding Makefile
# Build with: make
# Clean with: make clean

# Compiler settings
FC = gfortran
F2PY = f2py
FFLAGS = -O3 -fPIC -fopenmp

# Paths
SRCDIR = ../../src
BUILDDIR = build

# Libraries
LIBS = -llapack -lblas -lgomp

# Module name
MODULE = hprmat_fortran

# Object files
OBJS = $(BUILDDIR)/precision.o \
       $(BUILDDIR)/rmat_solvers.o \
       $(BUILDDIR)/rmatrix_hp.o \
       $(BUILDDIR)/special_functions.o \
       $(BUILDDIR)/angular_momentum.o

# Default target
all: $(MODULE)
	@echo ""
	@echo "Build complete! To test:"
	@echo "  python -c \"from hprmat import RMatrixSolver; print('OK')\""

# Compile Fortran source files to object files
$(BUILDDIR)/precision.o: $(SRCDIR)/precision.F90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILDDIR)

$(BUILDDIR)/rmat_solvers.o: $(SRCDIR)/rmat_solvers.F90 $(BUILDDIR)/precision.o | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/rmatrix_hp.o: $(SRCDIR)/rmatrix_hp.F90 $(BUILDDIR)/precision.o $(BUILDDIR)/rmat_solvers.o | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

$(BUILDDIR)/special_functions.o: $(SRCDIR)/special_functions.f | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILDDIR)/angular_momentum.o: $(SRCDIR)/angular_momentum.f | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Build the f2py module - only parse hprmat_f2py.F90, link with precompiled objects
$(MODULE): $(OBJS) hprmat_f2py.F90
	$(F2PY) -c hprmat_f2py.F90 $(OBJS) \
		-m $(MODULE) \
		--f90flags="$(FFLAGS)" \
		-I$(BUILDDIR) \
		$(LIBS) \
		--build-dir $(BUILDDIR)

# Create build directory
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Clean
clean:
	rm -rf $(BUILDDIR)
	rm -f $(MODULE)*.so
	rm -f *.pyc
	rm -rf __pycache__

# Install to user site-packages
install: $(MODULE)
	pip install --user .

# Test
test: $(MODULE)
	python -c "import $(MODULE); print('Import successful!')"
	python test_ex1.py

.PHONY: all clean install test
