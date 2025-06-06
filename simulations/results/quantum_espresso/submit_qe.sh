#!/bin/bash
#SBATCH --job-name=mgfe_dft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --partition=compute

module load quantum-espresso

for input in *.in; do
    echo "Running $input..."
    mpirun -np 16 pw.x < $input > ${input%.in}.out
done
