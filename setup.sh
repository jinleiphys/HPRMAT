#!/bin/bash
#==============================================================================
# HPRMAT Smart Setup Script
# Automatically detects GPU, CUDA, OpenBLAS and offers to install missing deps
#==============================================================================

set -e  # Exit on error

echo "=============================================="
echo "     HPRMAT Smart Setup"
echo "=============================================="
echo ""

# Output file
MAKEINC="make.inc"

# Minimum CUDA version for multi-GPU support
MIN_CUDA_VERSION=12

# Detect OS
OS=$(uname -s)
echo "[1/6] Operating System: $OS"

#------------------------------------------------------------------------------
# Helper functions
#------------------------------------------------------------------------------
ask_yes_no() {
    local prompt="$1"
    local default="${2:-y}"
    local answer

    if [ "$default" = "y" ]; then
        prompt="$prompt [Y/n]: "
    else
        prompt="$prompt [y/N]: "
    fi

    read -p "$prompt" answer
    answer=${answer:-$default}

    case "$answer" in
        [Yy]*) return 0 ;;
        *) return 1 ;;
    esac
}

version_ge() {
    # Returns 0 if $1 >= $2
    [ "$(printf '%s\n' "$2" "$1" | sort -V | head -n1)" = "$2" ]
}

#------------------------------------------------------------------------------
# Step 2: Detect GPU
#------------------------------------------------------------------------------
echo ""
echo "[2/6] Checking GPU..."

GPU_ENABLED="false"
GPU_ARCH="sm_70"
NUM_GPUS=0
MULTI_GPU="false"

if command -v nvidia-smi &> /dev/null; then
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)
    if [ -n "$GPU_NAME" ]; then
        NUM_GPUS=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | wc -l)
        echo "  GPU detected: $GPU_NAME"
        echo "  Number of GPUs: $NUM_GPUS"

        # Map GPU name to architecture
        case "$GPU_NAME" in
            *"RTX 50"*|*"RTX 5"*)
                GPU_ARCH="sm_120"  # Blackwell
                ;;
            *"RTX 40"*|*"RTX 4"*)
                GPU_ARCH="sm_89"   # Ada Lovelace
                ;;
            *"RTX 30"*|*"A6000"*|*"RTX 3"*)
                GPU_ARCH="sm_86"   # Ampere (consumer)
                ;;
            *"A100"*)
                GPU_ARCH="sm_80"   # Ampere (data center)
                ;;
            *"H100"*|*"H200"*)
                GPU_ARCH="sm_90"   # Hopper
                ;;
            *"V100"*)
                GPU_ARCH="sm_70"   # Volta
                ;;
            *"RTX 20"*|*"T4"*|*"RTX 2"*)
                GPU_ARCH="sm_75"   # Turing
                ;;
            *)
                GPU_ARCH="sm_70"   # Default fallback
                ;;
        esac
        echo "  GPU architecture: $GPU_ARCH"
        GPU_ENABLED="true"

        if [ "$NUM_GPUS" -ge 2 ]; then
            MULTI_GPU="true"
            echo "  Multi-GPU: enabled ($NUM_GPUS GPUs)"
        fi
    fi
else
    echo "  No GPU detected (nvidia-smi not found)"
fi

#------------------------------------------------------------------------------
# Step 3: Detect/Install CUDA
#------------------------------------------------------------------------------
echo ""
echo "[3/6] Checking CUDA..."

CUDA_PATH=""
CUDA_VERSION=""
USE_CONDA_CUDA="false"
CONDA_ENV_NAME="cuda12"

# Function to find CUDA
find_cuda() {
    local cuda_path=""
    local cuda_ver=""

    # Check conda environment first
    if [ -n "$CONDA_PREFIX" ] && [ -f "$CONDA_PREFIX/bin/nvcc" ]; then
        cuda_path="$CONDA_PREFIX"
        cuda_ver=$("$cuda_path/bin/nvcc" --version 2>/dev/null | grep "release" | sed 's/.*release \([0-9]*\.[0-9]*\).*/\1/')
        echo "$cuda_path|$cuda_ver|conda"
        return
    fi

    # Check if nvcc is in PATH
    if command -v nvcc &> /dev/null; then
        cuda_path=$(dirname $(dirname $(which nvcc)))
        cuda_ver=$(nvcc --version 2>/dev/null | grep "release" | sed 's/.*release \([0-9]*\.[0-9]*\).*/\1/')
        echo "$cuda_path|$cuda_ver|system"
        return
    fi

    # Check common CUDA installation paths
    for p in "/usr/local/cuda" "/usr/local/cuda-12" "/usr/local/cuda-11" "/opt/cuda"; do
        if [ -f "$p/bin/nvcc" ]; then
            cuda_path="$p"
            cuda_ver=$("$p/bin/nvcc" --version 2>/dev/null | grep "release" | sed 's/.*release \([0-9]*\.[0-9]*\).*/\1/')
            echo "$cuda_path|$cuda_ver|system"
            return
        fi
    done

    echo "||none"
}

# Function to setup conda
setup_conda() {
    # Check if conda exists
    if command -v conda &> /dev/null; then
        return 0
    fi

    # Check common conda locations
    for conda_init in "$HOME/miniconda3/etc/profile.d/conda.sh" \
                      "$HOME/anaconda3/etc/profile.d/conda.sh" \
                      "/opt/conda/etc/profile.d/conda.sh"; do
        if [ -f "$conda_init" ]; then
            source "$conda_init"
            return 0
        fi
    done

    return 1
}

# Function to install miniconda
install_miniconda() {
    echo "  Installing Miniconda..."
    local MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    if [ "$OS" = "Darwin" ]; then
        if [ "$(uname -m)" = "arm64" ]; then
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
        else
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
        fi
    fi

    curl -fsSL "$MINICONDA_URL" -o /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
    rm /tmp/miniconda.sh

    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda init bash
    echo "  Miniconda installed successfully!"
}

# Function to install CUDA via conda
install_cuda_conda() {
    local version="${1:-12.9}"
    echo "  Installing CUDA $version via conda..."

    # Check if environment exists
    if conda env list | grep -q "^$CONDA_ENV_NAME "; then
        echo "  Environment '$CONDA_ENV_NAME' already exists"
        conda activate "$CONDA_ENV_NAME"
    else
        echo "  Creating conda environment '$CONDA_ENV_NAME'..."
        conda create -n "$CONDA_ENV_NAME" -c nvidia cuda-toolkit="$version" -y
        conda activate "$CONDA_ENV_NAME"
    fi

    echo "  CUDA $version installed successfully!"
}

# Main CUDA detection logic
if [ "$GPU_ENABLED" = "true" ]; then
    CUDA_INFO=$(find_cuda)
    CUDA_PATH=$(echo "$CUDA_INFO" | cut -d'|' -f1)
    CUDA_VERSION=$(echo "$CUDA_INFO" | cut -d'|' -f2)
    CUDA_SOURCE=$(echo "$CUDA_INFO" | cut -d'|' -f3)

    if [ -n "$CUDA_PATH" ] && [ -n "$CUDA_VERSION" ]; then
        echo "  CUDA found: $CUDA_PATH (version $CUDA_VERSION, $CUDA_SOURCE)"

        # Check if version is sufficient for multi-GPU
        CUDA_MAJOR=$(echo "$CUDA_VERSION" | cut -d'.' -f1)
        if [ "$MULTI_GPU" = "true" ] && [ "$CUDA_MAJOR" -lt "$MIN_CUDA_VERSION" ]; then
            echo ""
            echo "  WARNING: CUDA $CUDA_VERSION is too old for stable multi-GPU support"
            echo "           Multi-GPU requires CUDA >= $MIN_CUDA_VERSION"
            echo ""

            if ask_yes_no "  Would you like to install CUDA 12.9 via conda?"; then
                if ! setup_conda; then
                    echo ""
                    if ask_yes_no "  Conda not found. Install Miniconda?"; then
                        install_miniconda
                    else
                        echo "  Skipping CUDA upgrade. Multi-GPU will be disabled."
                    fi
                fi

                if command -v conda &> /dev/null; then
                    install_cuda_conda "12.9"
                    USE_CONDA_CUDA="true"
                    CUDA_PATH="$CONDA_PREFIX"
                    CUDA_VERSION="12.9"
                fi
            else
                echo "  Keeping CUDA $CUDA_VERSION. Multi-GPU will be disabled."
                MULTI_GPU="false"
            fi
        fi
    else
        echo "  CUDA not found"
        echo ""

        if ask_yes_no "  Would you like to install CUDA 12.9 via conda?"; then
            if ! setup_conda; then
                echo ""
                if ask_yes_no "  Conda not found. Install Miniconda?"; then
                    install_miniconda
                else
                    echo "  Skipping CUDA installation. GPU support disabled."
                    GPU_ENABLED="false"
                fi
            fi

            if command -v conda &> /dev/null; then
                install_cuda_conda "12.9"
                USE_CONDA_CUDA="true"
                CUDA_PATH="$CONDA_PREFIX"
                CUDA_VERSION="12.9"
            fi
        else
            echo "  GPU support disabled."
            GPU_ENABLED="false"
        fi
    fi
fi

#------------------------------------------------------------------------------
# Step 4: Detect BLAS library
#------------------------------------------------------------------------------
echo ""
echo "[4/6] Checking BLAS library..."

BLAS_LIBS=""

if [ "$OS" = "Darwin" ]; then
    # macOS
    if [ -d "/opt/homebrew/opt/openblas" ]; then
        OPENBLAS_DIR="/opt/homebrew/opt/openblas"
        BLAS_LIBS="-L${OPENBLAS_DIR}/lib -lopenblas"
        echo "  OpenBLAS found: $OPENBLAS_DIR"
    elif [ -d "/usr/local/opt/openblas" ]; then
        OPENBLAS_DIR="/usr/local/opt/openblas"
        BLAS_LIBS="-L${OPENBLAS_DIR}/lib -lopenblas"
        echo "  OpenBLAS found: $OPENBLAS_DIR"
    else
        BLAS_LIBS="-framework Accelerate"
        echo "  Using Apple Accelerate framework"
    fi
else
    # Linux
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
            echo "  OpenBLAS found: $OPATH"
            OPENBLAS_FOUND="true"
            break
        fi
    done

    if [ "$OPENBLAS_FOUND" = "false" ]; then
        # Try ldconfig
        OPENBLAS_LIB=$(ldconfig -p 2>/dev/null | grep libopenblas | head -1 | awk '{print $NF}')
        if [ -n "$OPENBLAS_LIB" ]; then
            OPENBLAS_DIR=$(dirname "$OPENBLAS_LIB")
            BLAS_LIBS="-L$OPENBLAS_DIR -lopenblas"
            echo "  OpenBLAS found via ldconfig: $OPENBLAS_DIR"
            OPENBLAS_FOUND="true"
        fi
    fi

    if [ "$OPENBLAS_FOUND" = "false" ]; then
        echo "  OpenBLAS not found"
        if ask_yes_no "  Would you like to install OpenBLAS?"; then
            if command -v apt &> /dev/null; then
                sudo apt update && sudo apt install -y libopenblas-openmp-dev
            elif command -v yum &> /dev/null; then
                sudo yum install -y openblas-devel
            elif command -v dnf &> /dev/null; then
                sudo dnf install -y openblas-devel
            fi

            # Search again
            for OPATH in "${OPENBLAS_PATHS[@]}"; do
                if [ -f "$OPATH/libopenblas.so" ]; then
                    BLAS_LIBS="-L$OPATH -lopenblas"
                    echo "  OpenBLAS installed: $OPATH"
                    OPENBLAS_FOUND="true"
                    break
                fi
            done
        fi

        if [ "$OPENBLAS_FOUND" = "false" ]; then
            BLAS_LIBS="-llapack -lblas"
            echo "  WARNING: Using system LAPACK/BLAS (slower)"
        fi
    fi
fi

#------------------------------------------------------------------------------
# Step 5: Detect Fortran compiler
#------------------------------------------------------------------------------
echo ""
echo "[5/6] Checking Fortran compiler..."

FC="gfortran"
if command -v gfortran &> /dev/null; then
    FC="gfortran"
    FC_VERSION=$(gfortran --version | head -1)
    echo "  Found: $FC_VERSION"
elif command -v ifort &> /dev/null; then
    FC="ifort"
    echo "  Found: ifort"
else
    echo "  ERROR: No Fortran compiler found!"
    echo "  Please install gfortran: sudo apt install gfortran"
    exit 1
fi

#------------------------------------------------------------------------------
# Step 6: Generate make.inc
#------------------------------------------------------------------------------
echo ""
echo "[6/6] Generating $MAKEINC..."

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
NUM_GPUS = $NUM_GPUS
MULTI_GPU = $MULTI_GPU
EOF
fi

cat >> $MAKEINC << EOF

# Libraries
LIBS = $BLAS_LIBS -fopenmp
EOF

if [ "$GPU_ENABLED" = "true" ]; then
    # Determine CUDA library path (lib64 for system, lib for conda)
    if [ -d "$CUDA_PATH/lib64" ]; then
        CUDA_LIB_PATH="\$(CUDA_PATH)/lib64"
    else
        CUDA_LIB_PATH="\$(CUDA_PATH)/lib"
    fi

    if [ "$MULTI_GPU" = "true" ]; then
cat >> $MAKEINC << EOF
CUDA_LIB_PATH = $CUDA_LIB_PATH
LIBS += -L\$(CUDA_LIB_PATH) -Wl,-rpath,\$(CUDA_LIB_PATH) -lcudart -lcusolver -lcusolverMg -lcublas -lstdc++
EOF
    else
cat >> $MAKEINC << EOF
CUDA_LIB_PATH = $CUDA_LIB_PATH
LIBS += -L\$(CUDA_LIB_PATH) -Wl,-rpath,\$(CUDA_LIB_PATH) -lcudart -lcusolver -lcublas -lstdc++
EOF
    fi
fi

# Generate activation script for conda CUDA
if [ "$USE_CONDA_CUDA" = "true" ]; then
cat > activate_cuda.sh << 'ACTIVATE_EOF'
#!/bin/bash
# Activate CUDA environment for HPRMAT
# Source this file before compiling: source activate_cuda.sh

ACTIVATE_EOF

cat >> activate_cuda.sh << EOF
CONDA_INIT="\$HOME/miniconda3/etc/profile.d/conda.sh"
if [ -f "\$CONDA_INIT" ]; then
    source "\$CONDA_INIT"
elif [ -f "\$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "\$HOME/anaconda3/etc/profile.d/conda.sh"
fi

conda activate $CONDA_ENV_NAME
echo "CUDA environment activated: \$CONDA_PREFIX"
echo "CUDA version: \$(nvcc --version | grep release)"
EOF
    chmod +x activate_cuda.sh
fi

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "     Configuration Summary"
echo "=============================================="
echo "  Fortran compiler: $FC"
echo "  BLAS library:     $BLAS_LIBS"
echo "  GPU enabled:      $GPU_ENABLED"
if [ "$GPU_ENABLED" = "true" ]; then
    echo "  GPU architecture: $GPU_ARCH"
    echo "  CUDA path:        $CUDA_PATH"
    echo "  CUDA version:     $CUDA_VERSION"
    echo "  Number of GPUs:   $NUM_GPUS"
    echo "  Multi-GPU:        $MULTI_GPU"
fi
echo ""
echo "  Configuration saved to: $MAKEINC"

if [ "$USE_CONDA_CUDA" = "true" ]; then
    echo ""
    echo "=============================================="
    echo "     IMPORTANT: Conda CUDA Environment"
    echo "=============================================="
    echo "  Before compiling, activate the CUDA environment:"
    echo ""
    echo "    source activate_cuda.sh"
    echo ""
    echo "  Or manually:"
    echo "    conda activate $CONDA_ENV_NAME"
    echo ""
fi

echo ""
echo "Next steps:"
if [ "$USE_CONDA_CUDA" = "true" ]; then
    echo "  source activate_cuda.sh  # Activate CUDA"
fi
echo "  make                      # Build library"
echo "  make examples             # Build examples"
echo ""
