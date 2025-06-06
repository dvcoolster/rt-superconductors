#!/bin/bash
#SBATCH -J EPW_RT_SC
#SBATCH -N 2
#SBATCH -n 64
#SBATCH --time=08:00:00
#SBATCH --mem=128GB
#SBATCH -p compute

# Load modules (adjust for your cluster)
module load quantum-espresso/7.2
module load intel-mpi/2021.4
module load intel/2021.4

# Set environment variables
export OMP_NUM_THREADS=1
export I_MPI_PIN=1

echo "Starting phonon calculation at $(date)"
cd ph/
mpirun -np $SLURM_NTASKS ph.x < ph.in > ph.out
cd ..

echo "Starting q2r and matdyn calculations at $(date)"
mpirun -np 1 q2r.x < q2r.in > q2r.out
mpirun -np 1 matdyn.x < matdyn.in > matdyn.out

echo "Starting EPW calculation at $(date)"
mpirun -np $SLURM_NTASKS epw.x < epw.in > epw.out

echo "EPW calculations completed at $(date)"

# Extract key results
if [ -f epw.out ]; then
    echo "=== EPW RESULTS ==="
    echo "Lambda (electron-phonon coupling):"
    grep "lambda" epw.out | tail -1
    
    echo "Logarithmic average frequency:"
    grep "omega_log" epw.out | tail -1
    
    echo "Critical temperature estimates:"
    grep -A 5 "Superconducting" epw.out
    
    # Calculate Tc using Allen-Dynes formula
    python3 -c "
import re
import math

# Read lambda and omega_log from epw.out
with open('epw.out', 'r') as f:
    content = f.read()

lambda_match = re.search(r'lambda\s*=\s*([\d.]+)', content)
omega_match = re.search(r'omega_log\s*=\s*([\d.]+)', content)

if lambda_match and omega_match:
    λ = float(lambda_match.group(1))
    ω = float(omega_match.group(1))  # in K
    μ_star = 0.10  # typical value
    
    # Allen-Dynes formula
    if λ > μ_star:
        f1 = (1 + (λ/2.46) * (1 + 3.8*μ_star))**(1/3)
        f2 = 1 + ((λ - μ_star*(1+λ))/(λ + 0.15))**2 / 15
        Tc = (f1 * f2 * ω / 1.2) * math.exp(-1.04 * (1+λ) / (λ - μ_star*(1+0.62*λ)))
        print(f'Allen-Dynes Tc = {Tc:.1f} K (λ={λ:.3f}, ω_log={ω:.1f} K)')
    else:
        print('λ too small for superconductivity')
else:
    print('Could not extract lambda and omega_log from epw.out')
"
else
    echo "ERROR: epw.out not found"
fi 