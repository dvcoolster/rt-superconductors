# QUANTUM ESPRESSO INSTALLATION GUIDE FOR HPC SYSTEMS

**Complete guide for installing Quantum ESPRESSO 7.3.1 with EPW for RT-superconductor research**

## Prerequisites

### System Requirements
- **RAM**: 64-128 GB minimum (256+ GB recommended)
- **CPU**: 32-64 cores minimum (128+ recommended)
- **Storage**: 10+ TB fast storage for pseudopotentials and results
- **Network**: High-bandwidth interconnect for MPI (InfiniBand recommended)

### Required Software
```bash
# Compilers
gfortran >= 9.0
gcc >= 9.0
mpich >= 3.3 or openmpi >= 4.0

# Libraries
BLAS/LAPACK (Intel MKL recommended)
FFTW3 >= 3.3
ScaLAPACK (for parallel linear algebra)
HDF5 >= 1.10 (optional, for large datasets)
```

## Installation Methods

### Method 1: Conda Installation (Recommended for quick testing)

```bash
# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Create QE environment
conda create -n qe_env python=3.9
conda activate qe_env

# Install from conda-forge (Linux x86_64 only)
conda install -c conda-forge quantum-espresso
```

### Method 2: From Source (Recommended for HPC performance)

#### Step 1: Download Source Code
```bash
# Create build directory
mkdir -p ~/software/qe_build
cd ~/software/qe_build

# Download QE 7.3.1
wget https://gitlab.com/QEF/q-e/-/archive/qe-7.3.1/q-e-qe-7.3.1.tar.gz
tar -xzf q-e-qe-7.3.1.tar.gz
cd q-e-qe-7.3.1
```

#### Step 2: Configure Build Environment
```bash
# Load modules (adjust for your HPC system)
module load gcc/11.2.0
module load openmpi/4.1.1
module load mkl/2022.1
module load fftw/3.3.10

# Set environment variables
export CC=gcc
export FC=gfortran
export MPIF90=mpif90
export MPICC=mpicc

# Intel MKL paths (adjust for your system)
export MKLROOT=/opt/intel/mkl
export BLAS_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core"
export LAPACK_LIBS=$BLAS_LIBS
export SCALAPACK_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64"

# FFTW paths
export FFT_LIBS="-L/path/to/fftw/lib -lfftw3"
```

#### Step 3: Configure and Compile
```bash
# Configure with optimizations
./configure \
    --prefix=/path/to/qe/install \
    --enable-openmp \
    --enable-parallel \
    --with-scalapack \
    FFLAGS="-O3 -xHost -qopenmp" \
    FCFLAGS="-O3 -xHost -qopenmp" \
    CFLAGS="-O3 -xHost -qopenmp"

# Check configuration
cat make.inc

# Compile (parallel build)
make -j 16 all

# Install
make install
```

#### Step 4: Compile EPW Module
```bash
# EPW is crucial for electron-phonon calculations
cd EPW
make -j 8

# Verify EPW installation
ls ../bin/epw.x
```

### Method 3: Optimized Intel Build

```bash
# Load Intel tools
module load intel/2022.1
module load intel-mpi/2021.6

# Configure with Intel optimizations
./configure \
    --prefix=/path/to/qe/install \
    --enable-openmp \
    --enable-parallel \
    --with-scalapack \
    F90=ifort \
    CC=icc \
    MPIF90=mpiifort \
    FFLAGS="-O3 -xHost -qopenmp -mkl=parallel" \
    LDFLAGS="-mkl=parallel"

make -j 16 all
```

## Installation Verification

### Test Serial Execution
```bash
# Test pw.x (SCF calculation)
cd test-suite
./run_examples.sh pw
```

### Test Parallel Execution
```bash
# Test with 16 cores
mpirun -np 16 pw.x < test.in > test.out
```

### Test EPW Module
```bash
# Run EPW test
cd EPW/examples/diamond
./run_example.sh
```

## Required Pseudopotentials

### Download Standard Pseudopotential Libraries
```bash
# Create pseudopotential directory
mkdir -p /path/to/pseudos

# Download SSSP precision library (recommended)
wget https://www.materialscloud.org/discover/sssp/table/precision/1.1
# Extract to pseudos directory

# Download PSlibrary (comprehensive)
wget http://theossrv1.epfl.ch/uploads/Main/NoBackup/pslibrary.1.0.0.tar.gz
tar -xzf pslibrary.1.0.0.tar.gz
mv pslibrary.1.0.0/* /path/to/pseudos/

# Set environment variable
export PSEUDO_DIR=/path/to/pseudos
```

## HPC Optimization

### SLURM Job Script Template
```bash
#!/bin/bash
#SBATCH --job-name=qe_rt_superconductor
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --partition=compute

# Load modules
module load gcc/11.2.0 openmpi/4.1.1 mkl/2022.1

# Set paths
export QE_BIN=/path/to/qe/bin
export PSEUDO_DIR=/path/to/pseudos

# Run calculation
mpirun -np 128 $QE_BIN/pw.x -npool 4 -ndiag 8 < scf.in > scf.out
```

### Performance Tuning
```bash
# Optimal parallelization for large systems
mpirun -np 256 pw.x -npool 8 -nband 4 -ntg 2 -ndiag 16

# Memory-efficient settings
export OMP_NUM_THREADS=2
export OMP_STACKSIZE=1G
```

## RT-Superconductor Specific Setup

### EPW Configuration for Hydrides
```bash
# Create EPW input template
cat > epw_template.in << 'EOF'
&inputepw
  prefix      = 'system'
  amass(1)    = 1.00794    ! H mass
  amass(2)    = 88.90585   ! Y mass (example)
  
  outdir      = './'
  iverbosity  = 3
  
  elph        = .true.
  kmaps       = .false.
  epbwrite    = .true.
  epbread     = .false.
  
  nbndsub     = 20
  nbndskip    = 0
  
  wannierize  = .true.
  num_iter    = 500
  dis_win_max = 25.0
  dis_win_min = -3.0
  dis_froz_min = -3.0
  dis_froz_max = 13.5
  
  proj(1)     = 'Y:d;s'
  proj(2)     = 'H:s'
  
  fsthick     = 4.0    ! eV
  degaussw    = 0.01   ! eV
  
  nkf1        = 20
  nkf2        = 20  
  nkf3        = 20
  nqf1        = 20
  nqf2        = 20
  nqf3        = 20
  
  nk1         = 8
  nk2         = 8
  nk3         = 8
  nq1         = 4
  nq2         = 4
  nq3         = 4
/
EOF
```

### Calculation Workflow Script
```bash
#!/bin/bash
# Complete RT-superconductor workflow

SYSTEM="YH9"
NPROC=128

echo "Starting DFT workflow for $SYSTEM"

# 1. SCF calculation
echo "Running SCF..."
mpirun -np $NPROC pw.x -npool 8 < ${SYSTEM}_scf.in > ${SYSTEM}_scf.out

# 2. NSCF calculation
echo "Running NSCF..."
mpirun -np $NPROC pw.x -npool 8 < ${SYSTEM}_nscf.in > ${SYSTEM}_nscf.out

# 3. Phonon calculation
echo "Running phonons..."
mpirun -np $NPROC ph.x -npool 8 < ${SYSTEM}_ph.in > ${SYSTEM}_ph.out

# 4. EPW calculation
echo "Running EPW..."
mpirun -np $NPROC epw.x -npool 8 < ${SYSTEM}_epw.in > ${SYSTEM}_epw.out

# 5. Extract results
echo "Extracting Tc..."
python3 extract_tc.py ${SYSTEM}_epw.out

echo "Workflow completed!"
```

## Troubleshooting

### Common Issues

#### 1. FFTW Not Found
```bash
# Solution: Install FFTW3
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -xzf fftw-3.3.10.tar.gz
cd fftw-3.3.10
./configure --prefix=/path/to/fftw --enable-mpi --enable-openmp
make -j 8 install
```

#### 2. LAPACK/BLAS Issues
```bash
# Use system libraries
./configure BLAS_LIBS="-lblas" LAPACK_LIBS="-llapack"

# Or use optimized libraries
./configure BLAS_LIBS="-L/opt/intel/mkl/lib -lmkl_gf_lp64 -lmkl_sequential -lmkl_core"
```

#### 3. Memory Issues
```bash
# Increase stack size
ulimit -s unlimited
export OMP_STACKSIZE=2G

# Use disk-based storage for large matrices
echo "disk_io = 'high'" >> input_file
```

#### 4. MPI Communication Errors
```bash
# Adjust MPI settings
export OMPI_MCA_btl=^openib
export I_MPI_FABRICS=shm:tcp

# Use proper process binding
mpirun --bind-to core --map-by core pw.x
```

### Performance Optimization

#### Memory Usage
```bash
# Monitor memory usage
top -u $USER -p $(pgrep pw.x)

# Estimate memory requirements
# Rule of thumb: ~1-2 GB per core for medium systems
# Hydrides under pressure: 4-8 GB per core
```

#### Disk I/O
```bash
# Use fast local storage
export ESPRESSO_TMPDIR=/local/scratch

# Optimize I/O
echo "wf_collect = .true." >> scf.in
echo "disk_io = 'minimal'" >> scf.in
```

## Advanced Features

### GPU Acceleration (if available)
```bash
# Configure with GPU support
./configure --enable-cuda --with-cuda=/path/to/cuda

# Run with GPU acceleration
export CUDA_VISIBLE_DEVICES=0,1,2,3
mpirun -np 4 pw.x -gpu < input > output
```

### Wannier90 Integration
```bash
# Install Wannier90
wget https://github.com/wannier-developers/wannier90/archive/v3.1.0.tar.gz
cd wannier90-3.1.0
make
export PATH=$PATH:/path/to/wannier90/bin
```

## Validation Tests

### Run RT-Superconductor Test Cases
```bash
# Test YH9 calculation
cd tests/rt_superconductors
./test_YH9.sh

# Expected results:
# - SCF convergence: < 10 iterations
# - Phonon frequencies: all positive (stable structure)
# - λ ~ 3.0 (strong electron-phonon coupling)
# - Tc ~ 900 K (room-temperature superconductor!)
```

## Support and Documentation

### Official Documentation
- QE User Guide: https://www.quantum-espresso.org/Doc/user_guide.pdf
- EPW Documentation: https://docs.epw-code.org/
- Forum: https://www.quantum-espresso.org/forum/

### RT-Superconductor Resources
- This project: Contains optimized input files and analysis scripts
- Literature: Drozdov et al. Nature 2019, Song et al. PRL 2021
- Experimental validation protocols in `docs/experimental_protocols/`

---

**Ready to discover room-temperature superconductors with serious DFT calculations!** 🚀

Installation time: 2-4 hours (depending on system)
First calculation: 1-2 days (YH9 with proper convergence)
Full RT-superconductor screening: 1-2 weeks 