#!/bin/bash
#==============================================================================
# HPRMAT Setup Script
# Automatically detects GPU and OpenBLAS, generates make.inc
#==============================================================================

echo "=== HPRMAT Setup ==="
echo ""

# Output file
MAKEINC="make.inc"

# Detect OS
OS=$(uname -s)
echo "Operating System: $OS"

#------------------------------------------------------------------------------
# Detect GPU
#------------------------------------------------------------------------------
GPU_ENABLED="false"
GPU_ARCH="sm_70"
CUDA_PATH=""

# Check if nvcc is in PATH
if command -v nvcc &> /dev/null; then
    CUDA_PATH=$(dirname $(dirname $(which nvcc)))
# Check common CUDA installation paths
elif [ -d "/usr/local/cuda" ]; then
    CUDA_PATH="/usr/local/cuda"
elif [ -d "/usr/local/cuda-12" ]; then
    CUDA_PATH="/usr/local/cuda-12"
elif [ -d "/usr/local/cuda-11" ]; then
    CUDA_PATH="/usr/local/cuda-11"
elif [ -d "/opt/cuda" ]; then
    CUDA_PATH="/opt/cuda"
fi

if [ -n "$CUDA_PATH" ] && [ -f "$CUDA_PATH/bin/nvcc" ]; then
    echo "CUDA found: $CUDA_PATH"

    # Get GPU architecture
    if command -v nvidia-smi &> /dev/null; then
        GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
        echo "GPU detected: $GPU_NAME"

        # Map GPU name to architecture
        case "$GPU_NAME" in
            *"RTX 50"*)
                GPU_ARCH="sm_120"  # Blackwell
                ;;
            *"RTX 40"*)
                GPU_ARCH="sm_89"   # Ada Lovelace
                ;;
            *"RTX 30"*|*"A6000"*)
                GPU_ARCH="sm_86"   # Ampere (consumer)
                ;;
            *"A100"*)
                GPU_ARCH="sm_80"   # Ampere (data center)
                ;;
            *"H100"*)
                GPU_ARCH="sm_90"   # Hopper
                ;;
            *"V100"*)
                GPU_ARCH="sm_70"   # Volta
                ;;
            *"RTX 20"*|*"T4"*)
                GPU_ARCH="sm_75"   # Turing
                ;;
            *)
                GPU_ARCH="sm_70"   # Default fallback
                ;;
        esac
        GPU_ENABLED="true"
        echo "GPU architecture: $GPU_ARCH"
    fi
else
    echo "CUDA not found (GPU support disabled)"
fi

#------------------------------------------------------------------------------
# Detect BLAS library
#------------------------------------------------------------------------------
BLAS_LIBS=""

if [ "$OS" = "Darwin" ]; then
    # macOS
    if [ -d "/opt/homebrew/opt/openblas" ]; then
        OPENBLAS_DIR="/opt/homebrew/opt/openblas"
        BLAS_LIBS="-L${OPENBLAS_DIR}/lib -lopenblas"
        echo "OpenBLAS found: $OPENBLAS_DIR"
    elif [ -d "/usr/local/opt/openblas" ]; then
        OPENBLAS_DIR="/usr/local/opt/openblas"
        BLAS_LIBS="-L${OPENBLAS_DIR}/lib -lopenblas"
        echo "OpenBLAS found: $OPENBLAS_DIR"
    else
        # Use Accelerate framework
        BLAS_LIBS="-framework Accelerate"
        echo "Using Apple Accelerate framework"
    fi
else
    # Linux - search for OpenBLAS in common locations
    OPENBLAS_FOUND="false"
    OPENBLAS_PATHS=(
        "/usr/lib/x86_64-linux-gnu/openblas-openmp"
        "/usr/lib/x86_64-linux-gnu/openblas-pthread"
        "/usr/lib/x86_64-linux-gnu"
        "/usr/lib64/openblas-openmp"
        "/usr/lib64/openblas-pthread"
        "/usr/lib64"
        "/usr/local/lib"
        "/opt/openblas/lib"
    )

    for OPATH in "${OPENBLAS_PATHS[@]}"; do
        if [ -f "$OPATH/libopenblas.so" ] || [ -f "$OPATH/libopenblas.a" ]; then
            BLAS_LIBS="-L$OPATH -lopenblas"
            echo "OpenBLAS found: $OPATH"
            OPENBLAS_FOUND="true"
            break
        fi
    done

    # Also check via ldconfig
    if [ "$OPENBLAS_FOUND" = "false" ]; then
        OPENBLAS_LIB=$(ldconfig -p 2>/dev/null | grep libopenblas | head -1 | awk '{print $NF}')
        if [ -n "$OPENBLAS_LIB" ]; then
            OPENBLAS_DIR=$(dirname "$OPENBLAS_LIB")
            BLAS_LIBS="-L$OPENBLAS_DIR -lopenblas"
            echo "OpenBLAS found via ldconfig: $OPENBLAS_DIR"
            OPENBLAS_FOUND="true"
        fi
    fi

    if [ "$OPENBLAS_FOUND" = "false" ]; then
        echo "OpenBLAS not found. Attempting to install..."
        # Try to install OpenBLAS
        if command -v apt &> /dev/null; then
            sudo apt update && sudo apt install -y libopenblas-openmp-dev
        elif command -v yum &> /dev/null; then
            sudo yum install -y openblas-devel
        elif command -v dnf &> /dev/null; then
            sudo dnf install -y openblas-devel
        fi

        # Search again after installation
        for OPATH in "${OPENBLAS_PATHS[@]}"; do
            if [ -f "$OPATH/libopenblas.so" ] || [ -f "$OPATH/libopenblas.a" ]; then
                BLAS_LIBS="-L$OPATH -lopenblas"
                echo "OpenBLAS installed successfully: $OPATH"
                OPENBLAS_FOUND="true"
                break
            fi
        done

        # Check via ldconfig again
        if [ "$OPENBLAS_FOUND" = "false" ]; then
            sudo ldconfig 2>/dev/null
            OPENBLAS_LIB=$(ldconfig -p 2>/dev/null | grep libopenblas | head -1 | awk '{print $NF}')
            if [ -n "$OPENBLAS_LIB" ]; then
                OPENBLAS_DIR=$(dirname "$OPENBLAS_LIB")
                BLAS_LIBS="-L$OPENBLAS_DIR -lopenblas"
                echo "OpenBLAS installed successfully: $OPENBLAS_DIR"
                OPENBLAS_FOUND="true"
            fi
        fi

        if [ "$OPENBLAS_FOUND" = "false" ]; then
            BLAS_LIBS="-llapack -lblas"
            echo "WARNING: OpenBLAS not available, using system LAPACK/BLAS (slow)"
            echo "Please install OpenBLAS manually for better performance:"
            echo "  Ubuntu/Debian: sudo apt install libopenblas-openmp-dev"
            echo "  CentOS/RHEL:   sudo yum install openblas-devel"
        fi
    fi
fi

#------------------------------------------------------------------------------
# Detect Fortran compiler
#------------------------------------------------------------------------------
FC="gfortran"
if command -v gfortran &> /dev/null; then
    FC="gfortran"
    echo "Fortran compiler: gfortran"
elif command -v ifort &> /dev/null; then
    FC="ifort"
    echo "Fortran compiler: ifort"
fi

#------------------------------------------------------------------------------
# Generate make.inc
#------------------------------------------------------------------------------
echo ""
echo "Generating $MAKEINC..."

cat > $MAKEINC << EOF
#==============================================================================
# HPRMAT Configuration (auto-generated by setup.sh)
# $(date)
#==============================================================================

# Compilers
FC = $FC
CC = gcc

# Fortran flags
FFLAGS = -O3 -fopenmp

# GPU configuration
GPU_ENABLED = $GPU_ENABLED
EOF

if [ "$GPU_ENABLED" = "true" ]; then
cat >> $MAKEINC << EOF
CUDA_PATH = $CUDA_PATH
NVCC = \$(CUDA_PATH)/bin/nvcc
NVCCFLAGS = -O3
GPU_ARCH = $GPU_ARCH
EOF
fi

cat >> $MAKEINC << EOF

# Libraries
LIBS = $BLAS_LIBS -fopenmp
EOF

if [ "$GPU_ENABLED" = "true" ]; then
cat >> $MAKEINC << EOF
LIBS += -L\$(CUDA_PATH)/lib64 -lcudart -lcusolver -lcublas -lstdc++
EOF
fi

echo ""
echo "=== Configuration Summary ==="
echo "Fortran compiler: $FC"
echo "BLAS library: $BLAS_LIBS"
echo "GPU enabled: $GPU_ENABLED"
if [ "$GPU_ENABLED" = "true" ]; then
    echo "GPU architecture: $GPU_ARCH"
    echo "CUDA path: $CUDA_PATH"
fi
echo ""
echo "Configuration saved to $MAKEINC"
echo ""
echo "Next steps:"
echo "  make          # Build library"
echo "  make examples # Build examples"
echo ""
