#!/bin/bash
#SBATCH -J QE_RT_SC
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=02:00:00
#SBATCH --mem=64GB
#SBATCH -p compute

# Load modules (adjust for your cluster)
module load quantum-espresso/7.2
module load intel-mpi/2021.4
module load intel/2021.4

# Set environment variables
export OMP_NUM_THREADS=1
export I_MPI_PIN=1

# Create temp directory
mkdir -p tmp

echo "Starting SCF calculation at $(date)"
mpirun -np $SLURM_NTASKS pw.x -npool 8 < scf.in > scf.out

echo "Starting NSCF calculation at $(date)"
mpirun -np $SLURM_NTASKS pw.x -npool 8 < nscf.in > nscf.out

echo "QE calculations completed at $(date)"

# Check for convergence
if grep -q "convergence has been achieved" scf.out; then
    echo "SCF converged successfully"
else
    echo "WARNING: SCF may not have converged properly"
fi

if grep -q "convergence has been achieved" nscf.out; then
    echo "NSCF converged successfully"
else
    echo "WARNING: NSCF may not have converged properly"
fi 