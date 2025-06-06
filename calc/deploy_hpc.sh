#!/bin/bash
# =================================================================
# HPC DEPLOYMENT SCRIPT FOR RT-SUPERCONDUCTOR DFT CALCULATIONS
# =================================================================
#
# This script deploys our complete RT-superconductor framework on 
# HPC systems with proper Quantum ESPRESSO installation.
#
# Usage: ./deploy_hpc.sh [cluster_name]
# Example: ./deploy_hpc.sh stampede2
#
# Author: RT-Superconductor Research Team
# Date: 2025-06-07

set -e  # Exit on any error

# =================================================================
# CONFIGURATION
# =================================================================

CLUSTER_NAME=${1:-"unknown_cluster"}
QE_VERSION="7.3.1"
PROJECT_NAME="rt-superconductors"
WORK_DIR="/scratch/${USER}/${PROJECT_NAME}"
INSTALL_DIR="${HOME}/software"

echo "🚀 DEPLOYING RT-SUPERCONDUCTOR FRAMEWORK ON ${CLUSTER_NAME}"
echo "================================================================"

# =================================================================
# SYSTEM DETECTION AND MODULE LOADING
# =================================================================

detect_system() {
    echo "🔍 Detecting HPC system..."
    
    if [[ -f /etc/stampede2-release ]]; then
        SYSTEM="stampede2"
        module load gcc/9.1.0 impi/19.0.7 mkl/19.0.7 cmake/3.24.2
        module load fftw3/3.3.8 hdf5/1.10.4
    elif [[ -f /etc/frontera-release ]]; then
        SYSTEM="frontera"
        module load gcc/9.1.0 impi/19.0.7 mkl/19.0.7
        module load cmake/3.24.2 fftw3/3.3.8
    elif [[ -f /etc/bridges2-release ]]; then
        SYSTEM="bridges2"
        module load gcc/10.2.0 openmpi/4.0.5 fftw/3.3.10
        module load cmake/3.18.4 hdf5/1.10.7
    elif command -v srun &> /dev/null; then
        SYSTEM="slurm_generic"
        echo "⚠️  Generic SLURM system detected"
        echo "   Please load appropriate modules manually:"
        echo "   - GCC/GFortran compiler"
        echo "   - MPI implementation"
        echo "   - FFTW3, LAPACK/BLAS"
        echo "   - CMake >= 3.13"
    else
        SYSTEM="generic"
        echo "❌ Unknown system. Please install dependencies manually."
        exit 1
    fi
    
    echo "✅ System detected: ${SYSTEM}"
}

# =================================================================
# QUANTUM ESPRESSO INSTALLATION
# =================================================================

install_quantum_espresso() {
    echo "🔧 Installing Quantum ESPRESSO ${QE_VERSION}..."
    
    mkdir -p "${INSTALL_DIR}"
    cd "${INSTALL_DIR}"
    
    # Download QE
    if [[ ! -f "qe-${QE_VERSION}.tar.gz" ]]; then
        echo "📥 Downloading Quantum ESPRESSO..."
        wget "https://gitlab.com/QEF/q-e/-/archive/qe-${QE_VERSION}/q-e-qe-${QE_VERSION}.tar.gz" \
             -O "qe-${QE_VERSION}.tar.gz"
    fi
    
    # Extract
    echo "📦 Extracting QE source..."
    tar -xzf "qe-${QE_VERSION}.tar.gz"
    cd "q-e-qe-${QE_VERSION}"
    
    # Configure based on system
    echo "⚙️  Configuring QE build..."
    case $SYSTEM in
        "stampede2"|"frontera")
            ./configure --enable-openmp --enable-parallel \
                       --with-scalapack=intel --with-elpa
            ;;
        "bridges2")
            ./configure --enable-openmp --enable-parallel
            ;;
        *)
            ./configure --enable-openmp --enable-parallel
            ;;
    esac
    
    # Compile
    echo "🔨 Compiling QE (this may take 30-60 minutes)..."
    make -j8 pw ph epw pp
    
    # Test installation
    echo "🧪 Testing QE installation..."
    if [[ -f "bin/pw.x" ]]; then
        echo "✅ QE installation successful!"
        export PATH="${INSTALL_DIR}/q-e-qe-${QE_VERSION}/bin:$PATH"
    else
        echo "❌ QE installation failed!"
        exit 1
    fi
}

# =================================================================
# PSEUDOPOTENTIAL SETUP
# =================================================================

setup_pseudopotentials() {
    echo "📚 Setting up pseudopotentials..."
    
    cd "${WORK_DIR}"
    mkdir -p pseudos
    cd pseudos
    
    # Download SSSP library
    echo "📥 Downloading SSSP pseudopotentials..."
    
    # Essential pseudopotentials for our systems
    declare -a PSEUDOS=(
        "Y.pbe-spn-kjpaw_psl.1.0.0.UPF"
        "H.pbe-kjpaw_psl.1.0.0.UPF" 
        "Sc.pbe-spn-kjpaw_psl.1.0.0.UPF"
        "Ca.pbe-spn-kjpaw_psl.1.0.0.UPF"
        "Li.pbe-s-kjpaw_psl.1.0.0.UPF"
        "B.pbe-n-kjpaw_psl.1.0.0.UPF"
        "Be.pbe-s-kjpaw_psl.1.0.0.UPF"
        "C.pbe-n-kjpaw_psl.1.0.0.UPF"
        "Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF"
        "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
    )
    
    for pseudo in "${PSEUDOS[@]}"; do
        if [[ ! -f "$pseudo" ]]; then
            echo "  Downloading $pseudo..."
            wget "https://pseudopotentials.quantum-espresso.org/upf_files/$pseudo"
        fi
    done
    
    echo "✅ Pseudopotentials ready!"
}

# =================================================================
# PROJECT STRUCTURE SETUP
# =================================================================

setup_project_structure() {
    echo "📁 Setting up project structure..."
    
    mkdir -p "${WORK_DIR}"
    cd "${WORK_DIR}"
    
    # Copy our calculation framework
    rsync -av "${HOME}/${PROJECT_NAME}/" "${WORK_DIR}/"
    
    # Create results directories
    mkdir -p {results/hpc_calculations,results/dft_convergence,results/phonon_data,results/epw_output}
    
    # Set up job submission scripts
    generate_slurm_scripts
    
    echo "✅ Project structure ready in ${WORK_DIR}"
}

# =================================================================
# SLURM JOB SCRIPT GENERATION
# =================================================================

generate_slurm_scripts() {
    echo "📝 Generating SLURM job scripts..."
    
    # Master job submission script
    cat > submit_all_calculations.sh << 'EOL'
#!/bin/bash
# =================================================================
# MASTER JOB SUBMISSION FOR RT-SUPERCONDUCTOR CALCULATIONS
# =================================================================

SYSTEMS=("YH9" "ScH9" "CaH10" "Li2BH12" "Be4CH10" "MgFe")

for system in "${SYSTEMS[@]}"; do
    echo "🚀 Submitting calculations for $system..."
    
    # SCF calculation
    job1=$(sbatch --parsable calc/tools/qe_scf.slurm $system)
    echo "  SCF job: $job1"
    
    # NSCF calculation (depends on SCF)
    job2=$(sbatch --parsable --dependency=afterok:$job1 calc/tools/qe_nscf.slurm $system)
    echo "  NSCF job: $job2"
    
    # Phonon calculation (depends on NSCF)  
    job3=$(sbatch --parsable --dependency=afterok:$job2 calc/tools/qe_phonon.slurm $system)
    echo "  Phonon job: $job3"
    
    # EPW calculation (depends on Phonon)
    job4=$(sbatch --parsable --dependency=afterok:$job3 calc/tools/qe_epw.slurm $system)
    echo "  EPW job: $job4"
    
    echo "  All jobs for $system submitted!"
done

echo "🎯 All RT-superconductor calculations submitted!"
EOL

    chmod +x submit_all_calculations.sh
    
    # Individual SLURM scripts
    generate_slurm_scf
    generate_slurm_nscf  
    generate_slurm_phonon
    generate_slurm_epw
}

generate_slurm_scf() {
    cat > calc/tools/qe_scf.slurm << 'EOL'
#!/bin/bash
#SBATCH --job-name=qe_scf
#SBATCH --partition=normal
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=56
#SBATCH --time=02:00:00
#SBATCH --output=results/hpc_calculations/%x_%j.out
#SBATCH --error=results/hpc_calculations/%x_%j.err

# Load modules (system-specific)
source calc/tools/load_modules.sh

SYSTEM=${1:-"YH9"}
echo "🔧 Running SCF calculation for $SYSTEM"

cd calc/$SYSTEM
mpirun -np $SLURM_NTASKS pw.x < scf.in > scf.out

if grep -q "convergence has been achieved" scf.out; then
    echo "✅ SCF calculation for $SYSTEM completed successfully"
else
    echo "❌ SCF calculation for $SYSTEM failed"
    exit 1
fi
EOL
}

generate_slurm_nscf() {
    cat > calc/tools/qe_nscf.slurm << 'EOL'
#!/bin/bash
#SBATCH --job-name=qe_nscf
#SBATCH --partition=normal
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=56
#SBATCH --time=01:00:00
#SBATCH --output=results/hpc_calculations/%x_%j.out
#SBATCH --error=results/hpc_calculations/%x_%j.err

source calc/tools/load_modules.sh

SYSTEM=${1:-"YH9"}
echo "🔧 Running NSCF calculation for $SYSTEM"

cd calc/$SYSTEM
mpirun -np $SLURM_NTASKS pw.x < nscf.in > nscf.out

if grep -q "convergence has been achieved" nscf.out; then
    echo "✅ NSCF calculation for $SYSTEM completed successfully"
else
    echo "❌ NSCF calculation for $SYSTEM failed"
    exit 1
fi
EOL
}

generate_slurm_phonon() {
    cat > calc/tools/qe_phonon.slurm << 'EOL'
#!/bin/bash
#SBATCH --job-name=qe_phonon
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=56
#SBATCH --time=08:00:00
#SBATCH --output=results/hpc_calculations/%x_%j.out
#SBATCH --error=results/hpc_calculations/%x_%j.err

source calc/tools/load_modules.sh

SYSTEM=${1:-"YH9"}
echo "🔧 Running phonon calculation for $SYSTEM"

cd calc/$SYSTEM/ph
mpirun -np $SLURM_NTASKS ph.x < ph.in > ph.out

# Check for completion
if grep -q "PHONON CALCULATION COMPLETED" ph.out; then
    echo "✅ Phonon calculation for $SYSTEM completed successfully"
    
    # Run q2r and matdyn
    mpirun -np $SLURM_NTASKS q2r.x < q2r.in > q2r.out
    mpirun -np $SLURM_NTASKS matdyn.x < matdyn.in > matdyn.out
    
else
    echo "❌ Phonon calculation for $SYSTEM failed"
    exit 1
fi
EOL
}

generate_slurm_epw() {
    cat > calc/tools/qe_epw.slurm << 'EOL'
#!/bin/bash
#SBATCH --job-name=qe_epw
#SBATCH --partition=normal
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=56
#SBATCH --time=12:00:00
#SBATCH --output=results/hpc_calculations/%x_%j.out
#SBATCH --error=results/hpc_calculations/%x_%j.err

source calc/tools/load_modules.sh

SYSTEM=${1:-"YH9"}
echo "🔧 Running EPW calculation for $SYSTEM"

cd calc/$SYSTEM
mpirun -np $SLURM_NTASKS epw.x < epw.in > epw.out

# Extract results
if grep -q "lambda" epw.out; then
    echo "✅ EPW calculation for $SYSTEM completed successfully"
    
    # Extract superconducting properties
    lambda=$(grep "lambda" epw.out | tail -1 | awk '{print $3}')
    omega_log=$(grep "omega_log" epw.out | tail -1 | awk '{print $3}')
    
    # Calculate Tc using our Allen-Dynes implementation
    python3 ../../tools/calculate_tc.py $lambda $omega_log > tc_result.dat
    
    echo "📊 Results for $SYSTEM:"
    echo "  λ = $lambda"
    echo "  ω_log = $omega_log K"
    cat tc_result.dat
    
else
    echo "❌ EPW calculation for $SYSTEM failed"
    exit 1
fi
EOL
}

# =================================================================
# MODULE LOADING SCRIPT
# =================================================================

create_module_script() {
    cat > calc/tools/load_modules.sh << 'EOL'
#!/bin/bash
# =================================================================
# MODULE LOADING FOR DIFFERENT HPC SYSTEMS
# =================================================================

if [[ -f /etc/stampede2-release ]]; then
    # TACC Stampede2
    module load gcc/9.1.0 impi/19.0.7 mkl/19.0.7
    module load cmake/3.24.2 fftw3/3.3.8 hdf5/1.10.4
    
elif [[ -f /etc/frontera-release ]]; then
    # TACC Frontera
    module load gcc/9.1.0 impi/19.0.7 mkl/19.0.7
    module load cmake/3.24.2 fftw3/3.3.8
    
elif [[ -f /etc/bridges2-release ]]; then
    # PSC Bridges2
    module load gcc/10.2.0 openmpi/4.0.5 fftw/3.3.10
    module load cmake/3.18.4 hdf5/1.10.7
    
else
    echo "⚠️  Unknown system - please load modules manually"
fi

# Set QE path
export PATH="${HOME}/software/q-e-qe-7.3.1/bin:$PATH"
EOL

    chmod +x calc/tools/load_modules.sh
}

# =================================================================
# MAIN DEPLOYMENT FUNCTION
# =================================================================

main() {
    echo "🎯 Starting HPC deployment for RT-superconductor calculations..."
    
    detect_system
    install_quantum_espresso
    setup_project_structure
    setup_pseudopotentials
    create_module_script
    
    echo ""
    echo "🚀 DEPLOYMENT COMPLETE!"
    echo "================================================================"
    echo "✅ Quantum ESPRESSO installed: ${INSTALL_DIR}/q-e-qe-${QE_VERSION}"
    echo "✅ Project setup: ${WORK_DIR}"
    echo "✅ SLURM scripts generated"
    echo ""
    echo "🎯 NEXT STEPS:"
    echo "1. cd ${WORK_DIR}"
    echo "2. ./submit_all_calculations.sh"
    echo "3. Monitor with: squeue -u \$USER"
    echo ""
    echo "📊 Expected runtime: 24-48 hours for all systems"
    echo "💾 Expected storage: ~500 GB"
    echo ""
    echo "🔥 BREAKTHROUGH AWAITS! 🔥"
}

# =================================================================
# ERROR HANDLING
# =================================================================

trap 'echo "❌ Deployment failed at line $LINENO"' ERR

# Run main function
main "$@" 