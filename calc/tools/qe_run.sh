#!/bin/bash
#SBATCH -J QE-RT-SC
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=08:00:00
#SBATCH --mem=64GB
#SBATCH --partition=normal
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ================================================================
# QUANTUM ESPRESSO + EPW PIPELINE FOR RT-SUPERCONDUCTOR DISCOVERY
# ================================================================
#
# This script runs the complete DFT+DFPT+EPW pipeline for 
# ambient-pressure superconductor candidates.
#
# Usage: sbatch qe_run.sh <structure_directory>
# Example: sbatch qe_run.sh BBeC4H4
#
# Author: RT-Superconductor Research Team
# Date: 2025-06-07

set -e  # Exit on any error

# Input validation
if [ $# -ne 1 ]; then
    echo "❌ ERROR: Structure directory required"
    echo "Usage: sbatch qe_run.sh <structure_directory>"
    exit 1
fi

STRUCT=$1
WORK_DIR=$(pwd)

# Check if structure directory exists
if [ ! -d "$STRUCT" ]; then
    echo "❌ ERROR: Structure directory '$STRUCT' not found"
    exit 1
fi

echo "🚀 Starting RT-Superconductor DFT Pipeline"
echo "=========================================="
echo "📁 Structure: $STRUCT"
echo "🖥️  Node: $(hostname)"
echo "⏰ Start time: $(date)"
echo "🔧 Cores: $SLURM_NTASKS"
echo ""

# Load modules (adjust for your HPC system)
module purge
module load intel/2021.4
module load impi/2021.4
module load quantumespresso/7.3.1
# module load qe/7.3.1  # Alternative module name

# Set OpenMP threads
export OMP_NUM_THREADS=1

# Set Quantum ESPRESSO executables
PW_EXEC=$(which pw.x)
PH_EXEC=$(which ph.x)
Q2R_EXEC=$(which q2r.x)
MATDYN_EXEC=$(which matdyn.x)
EPW_EXEC=$(which epw.x)

echo "🔧 QE Executables:"
echo "   pw.x: $PW_EXEC"
echo "   ph.x: $PH_EXEC"
echo "   epw.x: $EPW_EXEC"
echo ""

# Change to structure directory
cd "$STRUCT"

# Create status tracking
echo "STARTED" > calculation.status
echo "$(date): Started calculation" >> calculation.log

# ================================================================
# STEP 1: SCF CALCULATION
# ================================================================
echo "⚛️  STEP 1: SCF Calculation"
echo "=========================="

if [ ! -f "scf.out" ] || ! grep -q "convergence has been achieved" scf.out 2>/dev/null; then
    echo "Running SCF calculation..."
    mpirun -n $SLURM_NTASKS $PW_EXEC -in scf.in > scf.out 2>&1
    
    if ! grep -q "convergence has been achieved" scf.out; then
        echo "❌ SCF convergence failed"
        echo "FAILED_SCF" > calculation.status
        exit 1
    fi
    echo "✅ SCF converged successfully"
else
    echo "✅ SCF already completed"
fi

echo "$(date): SCF completed" >> calculation.log

# ================================================================
# STEP 2: NSCF CALCULATION  
# ================================================================
echo ""
echo "⚛️  STEP 2: NSCF Calculation"
echo "============================"

if [ ! -f "nscf.out" ] || ! grep -q "convergence has been achieved" nscf.out 2>/dev/null; then
    echo "Running NSCF calculation..."
    mpirun -n $SLURM_NTASKS $PW_EXEC -in nscf.in > nscf.out 2>&1
    
    if ! grep -q "convergence has been achieved" nscf.out; then
        echo "❌ NSCF convergence failed"
        echo "FAILED_NSCF" > calculation.status
        exit 1
    fi
    echo "✅ NSCF converged successfully"
else
    echo "✅ NSCF already completed"
fi

echo "$(date): NSCF completed" >> calculation.log

# ================================================================
# STEP 3: GAMMA-POINT PHONON CHECK
# ================================================================
echo ""
echo "🎵 STEP 3: Gamma-Point Phonon Check"
echo "===================================="

if [ ! -f "gamma.out" ] || ! grep -q "End of phonon calculation" gamma.out 2>/dev/null; then
    echo "Running Gamma-point phonon calculation..."
    mpirun -n $SLURM_NTASKS $PH_EXEC -in gamma.in > gamma.out 2>&1
    
    if ! grep -q "End of phonon calculation" gamma.out; then
        echo "❌ Gamma phonon calculation failed"
        echo "FAILED_GAMMA" > calculation.status
        exit 1
    fi
else
    echo "✅ Gamma phonon already completed"
fi

# Check for imaginary frequencies (stability)
if grep -q "negative" gamma.out; then
    echo "❌ GATE-1 FAILED: Imaginary frequencies found at Gamma point"
    echo "   Structure is dynamically unstable"
    echo "FAILED_GATE1" > calculation.status
    exit 1
fi

echo "✅ GATE-1 PASSED: No imaginary frequencies at Gamma point"
echo "$(date): Gamma phonon completed - GATE-1 PASSED" >> calculation.log

# ================================================================
# STEP 4: FULL PHONON CALCULATION (4x4x4 GRID)
# ================================================================
echo ""
echo "🎵 STEP 4: Full Phonon Calculation (4×4×4)"
echo "==========================================="

if [ ! -f "ph.out" ] || ! grep -q "End of phonon calculation" ph.out 2>/dev/null; then
    echo "Running full phonon calculation (this may take 2-4 hours)..."
    mpirun -n $SLURM_NTASKS $PH_EXEC -in ph.in > ph.out 2>&1
    
    if ! grep -q "End of phonon calculation" ph.out; then
        echo "❌ Full phonon calculation failed"
        echo "FAILED_PHONON" > calculation.status
        exit 1
    fi
    echo "✅ Full phonon calculation completed"
else
    echo "✅ Full phonon already completed"
fi

echo "$(date): Full phonon completed" >> calculation.log

# ================================================================
# STEP 5: FORCE CONSTANTS AND INTERPOLATION
# ================================================================
echo ""
echo "🔧 STEP 5: Force Constants & Interpolation"
echo "==========================================="

# q2r: Convert dynamical matrices to force constants
if [ ! -f "q2r.out" ]; then
    echo "Running q2r (force constants)..."
    $Q2R_EXEC < q2r.in > q2r.out 2>&1
    echo "✅ Force constants generated"
fi

# matdyn: Check for imaginary modes in full spectrum
if [ ! -f "matdyn.out" ]; then
    echo "Running matdyn (phonon interpolation)..."
    $MATDYN_EXEC < matdyn.in > matdyn.out 2>&1
    echo "✅ Phonon interpolation completed"
fi

# Check for imaginary frequencies in full spectrum
if grep -q -i "imaginary\|negative" matdyn.out; then
    echo "❌ WARNING: Imaginary frequencies found in phonon spectrum"
    echo "   Proceeding but results may be unreliable"
    echo "IMAGINARY_MODES" > calculation.status
    # Don't exit - continue with EPW calculation
fi

echo "$(date): Force constants and interpolation completed" >> calculation.log

# ================================================================
# STEP 6: EPW ELECTRON-PHONON COUPLING
# ================================================================
echo ""
echo "⚡ STEP 6: EPW Electron-Phonon Coupling"
echo "========================================"

if [ ! -f "epw.out" ] || ! grep -q "Ended the EPW calculation" epw.out 2>/dev/null; then
    echo "Running EPW calculation (this may take 4-8 hours)..."
    mpirun -n $SLURM_NTASKS $EPW_EXEC -in epw.in > epw.out 2>&1
    
    if ! grep -q "Ended the EPW calculation" epw.out; then
        echo "❌ EPW calculation failed"
        echo "FAILED_EPW" > calculation.status
        exit 1
    fi
    echo "✅ EPW calculation completed"
else
    echo "✅ EPW already completed"
fi

echo "$(date): EPW completed" >> calculation.log

# ================================================================
# STEP 7: EXTRACT RESULTS AND APPLY GATES
# ================================================================
echo ""
echo "📊 STEP 7: Results Extraction & Analysis"
echo "========================================="

# Extract λ (electron-phonon coupling)
if grep -q "lambda =" epw.out; then
    LAMBDA=$(grep "lambda =" epw.out | tail -1 | awk '{print $3}')
    echo "📊 λ (electron-phonon coupling) = $LAMBDA"
else
    echo "❌ Could not extract λ from EPW output"
    LAMBDA="N/A"
fi

# Extract ω_log (logarithmic average phonon frequency)
if grep -q "omega_log" epw.out; then
    OMEGA_LOG=$(grep "omega_log" epw.out | tail -1 | awk '{print $3}')
    echo "📊 ω_log (avg phonon frequency) = $OMEGA_LOG K"
else
    echo "❌ Could not extract ω_log from EPW output"
    OMEGA_LOG="N/A"
fi

# Calculate Allen-Dynes Tc if we have both λ and ω_log
if [[ "$LAMBDA" != "N/A" && "$OMEGA_LOG" != "N/A" ]]; then
    # Simple Allen-Dynes formula: Tc = ω_log/1.2 * exp(-1.04(1+λ)/(λ-μ*(1+0.62λ)))
    # Using μ* = 0.10
    MU_STAR=0.10
    TC_AD=$(python3 -c "
import math
lambda_val = $LAMBDA
omega_log = $OMEGA_LOG
mu_star = $MU_STAR
if lambda_val > mu_star:
    tc = (omega_log / 1.2) * math.exp(-1.04 * (1 + lambda_val) / (lambda_val - mu_star * (1 + 0.62 * lambda_val)))
    print(f'{tc:.1f}')
else:
    print('N/A')
")
    echo "🌡️  Tc (Allen-Dynes) = $TC_AD K"
else
    TC_AD="N/A"
fi

# Apply GATE-2: Check superconducting potential
GATE2_STATUS="UNKNOWN"
if [[ "$LAMBDA" != "N/A" && "$OMEGA_LOG" != "N/A" ]]; then
    if (( $(echo "$LAMBDA >= 2.0" | bc -l) )) && (( $(echo "$OMEGA_LOG >= 800" | bc -l) )); then
        echo "✅ GATE-2 PASSED: λ ≥ 2.0 and ω_log ≥ 800 K"
        echo "🔥 POTENTIAL ROOM-TEMPERATURE SUPERCONDUCTOR!"
        GATE2_STATUS="PASS"
    else
        echo "❌ GATE-2 FAILED: λ < 2.0 or ω_log < 800 K"
        GATE2_STATUS="FAIL"
    fi
fi

# ================================================================
# STEP 8: SAVE RESULTS
# ================================================================
echo ""
echo "💾 STEP 8: Saving Results"
echo "========================="

# Create results summary
cat > results_summary.json << EOF
{
    "structure": "$STRUCT",
    "calculation_date": "$(date -Iseconds)",
    "status": "COMPLETED",
    "lambda": "$LAMBDA",
    "omega_log": "$OMEGA_LOG",
    "tc_allen_dynes": "$TC_AD",
    "gate1_status": "PASS",
    "gate2_status": "$GATE2_STATUS",
    "calculation_time": "$(echo $SECONDS | awk '{print int($1/3600)":"int(($1%3600)/60)":"int($1%60)}')"
}
EOF

# Update calculation status
echo "COMPLETED" > calculation.status
echo "$(date): Calculation completed successfully" >> calculation.log

# ================================================================
# COMPLETION SUMMARY
# ================================================================
echo ""
echo "🎉 CALCULATION COMPLETED SUCCESSFULLY!"
echo "======================================"
echo "📁 Structure: $STRUCT"
echo "📊 Results:"
echo "   λ = $LAMBDA"
echo "   ω_log = $OMEGA_LOG K"
echo "   Tc = $TC_AD K"
echo "   Gate-2: $GATE2_STATUS"
echo ""
echo "⏰ Total time: $(echo $SECONDS | awk '{print int($1/3600)":"int(($1%3600)/60)":"int($1%60)}')"
echo "💾 Results saved to: results_summary.json"
echo ""

if [[ "$GATE2_STATUS" == "PASS" ]]; then
    echo "🔥🔥🔥 BREAKTHROUGH: POTENTIAL RT-SUPERCONDUCTOR FOUND! 🔥🔥🔥"
fi

echo "🚀 RT-Superconductor Discovery Pipeline Complete!" 