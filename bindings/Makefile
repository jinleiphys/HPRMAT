# HPRMAT Language Bindings Makefile
# Build shared library and language bindings
#
# Usage:
#   make lib       - Build shared library (libhprmat.so/dylib)
#   make python    - Build Python f2py module
#   make all       - Build everything
#   make clean     - Clean build artifacts
#   make install   - Install to system (requires sudo)

# Compiler settings
FC = gfortran
CC = gcc
F2PY = f2py

# Flags
FFLAGS = -O3 -fPIC -fopenmp
CFLAGS = -O3 -fPIC

# Platform detection
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    SHARED_EXT = dylib
    SHARED_FLAGS = -dynamiclib
    RPATH_FLAG = -install_name @rpath/libhprmat.$(SHARED_EXT)
else
    SHARED_EXT = so
    SHARED_FLAGS = -shared
    RPATH_FLAG = -Wl,-soname,libhprmat.$(SHARED_EXT)
endif

# Directories
SRCDIR = ../src
LIBDIR = ../lib
BUILDDIR = build

# Source files
FORTRAN_SRCS = $(SRCDIR)/precision.F90 \
               $(SRCDIR)/rmat_solvers.F90 \
               $(SRCDIR)/rmatrix_hp.F90 \
               $(SRCDIR)/special_functions.f \
               $(SRCDIR)/angular_momentum.f

C_INTERFACE = c/hprmat_c.F90

# Object files
OBJS = $(BUILDDIR)/precision.o \
       $(BUILDDIR)/rmat_solvers.o \
       $(BUILDDIR)/rmatrix_hp.o \
       $(BUILDDIR)/special_functions.o \
       $(BUILDDIR)/angular_momentum.o \
       $(BUILDDIR)/hprmat_c.o

# Libraries
LIBS = -llapack -lblas -lgomp

# Targets
.PHONY: all lib python clean install

all: lib python

# Build shared library
lib: $(LIBDIR)/libhprmat.$(SHARED_EXT)

$(LIBDIR)/libhprmat.$(SHARED_EXT): $(OBJS) | $(LIBDIR)
	$(FC) $(SHARED_FLAGS) $(RPATH_FLAG) -o $@ $(OBJS) $(LIBS)
	@echo "Shared library built: $@"

# Object file rules
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

$(BUILDDIR)/hprmat_c.o: $(C_INTERFACE) $(BUILDDIR)/precision.o $(BUILDDIR)/rmat_solvers.o $(BUILDDIR)/rmatrix_hp.o | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(BUILDDIR) -J$(BUILDDIR)

# Create directories
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

# Build Python f2py module
python: lib
	cd python && $(MAKE)

# Clean
clean:
	rm -rf $(BUILDDIR)
	rm -f $(LIBDIR)/libhprmat.*
	rm -f *.mod *.o
	cd python && $(MAKE) clean 2>/dev/null || true

# Install (requires sudo)
install: lib
	install -d /usr/local/lib
	install -d /usr/local/include
	install -m 644 $(LIBDIR)/libhprmat.$(SHARED_EXT) /usr/local/lib/
	install -m 644 c/hprmat.h /usr/local/include/
	@echo "Installed to /usr/local"
ifeq ($(UNAME_S),Linux)
	ldconfig
endif

# Test
test: lib python
	@echo "Testing C library..."
	cd c && $(MAKE) test 2>/dev/null || echo "C test not available"
	@echo "Testing Python..."
	cd python && python -c "from hprmat import RMatrixSolver; print('Python import OK')"
	@echo "Testing Julia..."
	cd julia && julia -e 'include("HPRMAT.jl"); using .HPRMAT; println("Julia import OK")' 2>/dev/null || echo "Julia test not available"
