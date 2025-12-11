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

# Source files needed
FORTRAN_SRCS = $(SRCDIR)/precision.F90 \
               $(SRCDIR)/rmat_solvers.F90 \
               $(SRCDIR)/rmatrix_hp.F90 \
               $(SRCDIR)/special_functions.f \
               $(SRCDIR)/angular_momentum.f

# Libraries
LIBS = -llapack -lblas -lgomp

# Module name
MODULE = hprmat_fortran

# Default target
all: $(MODULE)

# Build the f2py module
$(MODULE): $(BUILDDIR) hprmat_f2py.F90
	$(F2PY) -c hprmat_f2py.F90 $(FORTRAN_SRCS) \
		-m $(MODULE) \
		--f90flags="$(FFLAGS)" \
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
	python example_alpha_alpha.py

.PHONY: all clean install test
