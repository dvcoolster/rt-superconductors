#!/usr/bin/env python
"""
High-accuracy DFT calculations using Quantum ESPRESSO for final superconductor screening.

This script performs density functional theory calculations on the most promising
candidates that passed xTB screening. It provides:
- Self-consistent field (SCF) calculations
- Geometry optimization with variable cell
- Electronic structure analysis (DOS, band structure)
- Phonon calculations for superconductivity analysis
- SLURM cluster submission capabilities

Usage:
    # Single structure
    python 04_qe_relax.py structure.cif --ecut 60 --kmesh 6 4 6

    # Batch processing
    python 04_qe_relax.py --batch candidates.csv --n_best 50 --slurm

    # Analysis mode
    python 04_qe_relax.py --analyze results_dir/ --out final_rankings.csv

Reference: Giannozzi, P. et al. J. Phys.: Condens. Matter 2009, 21, 395502
"""

import argparse
import json
import os
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import numpy as np
import pandas as pd
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Materials science libraries
from pymatgen.core import Structure, Composition
from pymatgen.io.pwscf import PWInput, PWOutput
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine

from utils.ledger import rbt_superconductor_score
from utils.io import (
    structure_from_file, structure_to_qe_input, 
    setup_calculation_directory, structure_to_cif
)


class QECalculation:
    """Wrapper for Quantum ESPRESSO calculations."""
    
    def __init__(self, 
                 work_dir: Path,
                 pw_command: str = 'pw.x',
                 ph_command: str = 'ph.x',
                 q2r_command: str = 'q2r.x',
                 matdyn_command: str = 'matdyn.x'):
        self.work_dir = Path(work_dir)
        self.pw_command = pw_command
        self.ph_command = ph_command
        self.q2r_command = q2r_command
        self.matdyn_command = matdyn_command
        
        # Create directory structure
        self.work_dir.mkdir(parents=True, exist_ok=True)
        for subdir in ['input', 'output', 'tmp', 'pseudos']:
            (self.work_dir / subdir).mkdir(exist_ok=True)
    
    def scf_calculation(self,
                       structure: Structure,
                       ecut: float = 60.0,
                       kpoints: List[int] = [6, 6, 6],
                       pseudopotentials: Optional[Dict[str, str]] = None,
                       timeout: int = 3600) -> Tuple[bool, Dict[str, Any]]:
        """
        Perform self-consistent field calculation.
        
        Args:
            structure: Input structure
            ecut: Plane wave cutoff in Ry
            kpoints: k-point mesh [nk1, nk2, nk3]
            pseudopotentials: Dict mapping element to pseudopotential file
            timeout: Calculation timeout in seconds
            
        Returns:
            Tuple of (success, results_dict)
        """
        # Set default pseudopotentials
        if pseudopotentials is None:
            pseudopotentials = self._get_default_pseudopotentials(structure)
        
        # Control parameters
        control = {
            'calculation': 'scf',
            'restart_mode': 'from_scratch',
            'outdir': str(self.work_dir / 'tmp'),
            'pseudo_dir': str(self.work_dir / 'pseudos'),
            'prefix': 'superconductor',
            'verbosity': 'high',
            'disk_io': 'low'
        }
        
        # System parameters  
        system = {
            'ecutwfc': ecut,
            'ecutrho': ecut * 4,  # Charge density cutoff
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02,
            'input_dft': 'PBE'  # Use PBE functional
        }
        
        # Electrons parameters
        electrons = {
            'conv_thr': 1.0e-8,
            'mixing_beta': 0.3,
            'mixing_mode': 'plain',
            'diagonalization': 'david'
        }
        
        # Create input file
        input_file = self.work_dir / 'input' / 'scf.in'
        structure_to_qe_input(
            structure, input_file,
            pseudopotentials=pseudopotentials,
            control=control,
            system=system,
            electrons=electrons,
            kpoints=kpoints
        )
        
        # Run calculation
        output_file = self.work_dir / 'output' / 'scf.out'
        return self._run_pw_calculation(input_file, output_file, timeout)
    
    def relax_calculation(self,
                         structure: Structure,
                         ecut: float = 60.0,
                         kpoints: List[int] = [4, 4, 4],
                         cell_dynamics: str = 'bfgs',
                         timeout: int = 7200) -> Tuple[bool, Dict[str, Any], Optional[Structure]]:
        """
        Perform geometry and cell optimization.
        
        Args:
            structure: Input structure
            ecut: Plane wave cutoff in Ry
            kpoints: k-point mesh
            cell_dynamics: Cell optimization algorithm
            timeout: Calculation timeout in seconds
            
        Returns:
            Tuple of (success, results_dict, optimized_structure)
        """
        # Control parameters for variable cell relaxation
        control = {
            'calculation': 'vc-relax',
            'restart_mode': 'from_scratch',
            'outdir': str(self.work_dir / 'tmp'),
            'pseudo_dir': str(self.work_dir / 'pseudos'),
            'prefix': 'superconductor',
            'verbosity': 'high',
            'forc_conv_thr': 1.0e-4,  # Force convergence threshold
            'etot_conv_thr': 1.0e-6   # Energy convergence threshold
        }
        
        # System parameters
        system = {
            'ecutwfc': ecut,
            'ecutrho': ecut * 4,
            'occupations': 'smearing', 
            'smearing': 'gaussian',
            'degauss': 0.02,
            'input_dft': 'PBE'
        }
        
        # Cell optimization parameters
        cell = {
            'cell_dynamics': cell_dynamics,
            'press_conv_thr': 0.1,  # Pressure convergence (kbar)
            'cell_dofree': 'all'    # Allow all cell degrees of freedom
        }
        
        # Ion optimization parameters  
        ions = {
            'ion_dynamics': 'bfgs'
        }
        
        # Create input
        input_file = self.work_dir / 'input' / 'relax.in'
        
        # Manual creation since we need CELL and IONS sections
        self._write_relax_input(structure, input_file, control, system, cell, ions, kpoints)
        
        # Run calculation
        output_file = self.work_dir / 'output' / 'relax.out'
        success, results = self._run_pw_calculation(input_file, output_file, timeout)
        
        # Extract optimized structure if successful
        optimized_structure = None
        if success:
            try:
                pw_output = PWOutput(str(output_file))
                if pw_output.final_structure:
                    optimized_structure = pw_output.final_structure
            except Exception as e:
                print(f"Warning: Could not extract optimized structure: {e}")
        
        return success, results, optimized_structure
    
    def dos_calculation(self,
                       ecut: float = 60.0,
                       kpoints: List[int] = [8, 8, 8],
                       timeout: int = 1800) -> Tuple[bool, Dict[str, Any]]:
        """
        Calculate density of states (requires prior SCF).
        
        Args:
            ecut: Plane wave cutoff
            kpoints: Dense k-point mesh for DOS
            timeout: Calculation timeout
            
        Returns:
            Tuple of (success, results_dict)
        """
        # NSCF calculation for DOS
        control = {
            'calculation': 'nscf',
            'restart_mode': 'restart',
            'outdir': str(self.work_dir / 'tmp'),
            'pseudo_dir': str(self.work_dir / 'pseudos'),
            'prefix': 'superconductor',
            'verbosity': 'high'
        }
        
        system = {
            'ecutwfc': ecut,
            'ecutrho': ecut * 4,
            'occupations': 'tetrahedra',  # Better for DOS
            'input_dft': 'PBE'
        }
        
        # Create dummy structure (will use from restart)
        dummy_structure = Structure.from_spacegroup("Pm-3m", [[1, 0, 0], [0, 1, 0], [0, 0, 1]], ["H"], [[0, 0, 0]])
        
        input_file = self.work_dir / 'input' / 'nscf.in'
        structure_to_qe_input(
            dummy_structure, input_file,
            control=control,
            system=system,
            kpoints=kpoints
        )
        
        output_file = self.work_dir / 'output' / 'nscf.out'
        return self._run_pw_calculation(input_file, output_file, timeout)
    
    def _get_default_pseudopotentials(self, structure: Structure) -> Dict[str, str]:
        """Get default pseudopotential mapping."""
        pseudos = {}
        for element in structure.composition.elements:
            # Use SSSP precision pseudopotentials (would need to be downloaded)
            pseudos[element.symbol] = f"{element.symbol}_ONCV_PBE-1.0.upf"
        return pseudos
    
    def _write_relax_input(self, structure, input_file, control, system, cell, ions, kpoints):
        """Write QE input file for relaxation with all necessary sections."""
        with open(input_file, 'w') as f:
            # Control section
            f.write("&CONTROL\n")
            for key, value in control.items():
                if isinstance(value, str):
                    f.write(f"  {key} = '{value}'\n")
                else:
                    f.write(f"  {key} = {value}\n")
            f.write("/\n\n")
            
            # System section
            f.write("&SYSTEM\n")
            f.write(f"  nat = {structure.num_sites}\n")
            f.write(f"  ntyp = {len(structure.composition.elements)}\n")
            f.write(f"  ibrav = 0\n")  # Use CELL_PARAMETERS
            
            for key, value in system.items():
                if isinstance(value, str):
                    f.write(f"  {key} = '{value}'\n")
                else:
                    f.write(f"  {key} = {value}\n")
            f.write("/\n\n")
            
            # Electrons section
            f.write("&ELECTRONS\n")
            f.write("  conv_thr = 1.0e-8\n")
            f.write("  mixing_beta = 0.3\n")
            f.write("/\n\n")
            
            # Ions section
            f.write("&IONS\n")
            for key, value in ions.items():
                f.write(f"  {key} = '{value}'\n")
            f.write("/\n\n")
            
            # Cell section
            f.write("&CELL\n")
            for key, value in cell.items():
                if isinstance(value, str):
                    f.write(f"  {key} = '{value}'\n")
                else:
                    f.write(f"  {key} = {value}\n")
            f.write("/\n\n")
            
            # Atomic species
            f.write("ATOMIC_SPECIES\n")
            for element in structure.composition.elements:
                mass = element.atomic_mass
                pseudo = f"{element.symbol}_ONCV_PBE-1.0.upf"
                f.write(f"{element.symbol} {mass:.3f} {pseudo}\n")
            f.write("\n")
            
            # Cell parameters
            f.write("CELL_PARAMETERS angstrom\n")
            for row in structure.lattice.matrix:
                f.write(f"{row[0]:12.8f} {row[1]:12.8f} {row[2]:12.8f}\n")
            f.write("\n")
            
            # Atomic positions
            f.write("ATOMIC_POSITIONS crystal\n")
            for site in structure:
                f.write(f"{site.specie.symbol} {site.frac_coords[0]:12.8f} {site.frac_coords[1]:12.8f} {site.frac_coords[2]:12.8f}\n")
            f.write("\n")
            
            # K-points
            f.write("K_POINTS automatic\n")
            f.write(f"{kpoints[0]} {kpoints[1]} {kpoints[2]} 0 0 0\n")
    
    def _run_pw_calculation(self, input_file: Path, output_file: Path, timeout: int) -> Tuple[bool, Dict[str, Any]]:
        """Run pw.x calculation."""
        cmd = [self.pw_command, '-input', str(input_file)]
        
        try:
            with open(output_file, 'w') as f:
                result = subprocess.run(
                    cmd, 
                    stdout=f, 
                    stderr=subprocess.STDOUT,
                    timeout=timeout,
                    cwd=self.work_dir
                )
            
            # Parse output
            success = result.returncode == 0
            results = self._parse_pw_output(output_file) if success else {}
            
            if not success:
                results['error'] = f"pw.x returned code {result.returncode}"
            
            return success, results
            
        except subprocess.TimeoutExpired:
            return False, {'error': 'Calculation timed out'}
        except Exception as e:
            return False, {'error': str(e)}
    
    def _parse_pw_output(self, output_file: Path) -> Dict[str, Any]:
        """Parse pw.x output file."""
        results = {}
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Look for final energy
            for line in content.split('\n'):
                if 'Final energy' in line or 'total energy' in line:
                    try:
                        energy_ry = float(line.split('=')[-1].split()[0])
                        results['total_energy_ry'] = energy_ry
                        results['total_energy_eV'] = energy_ry * 13.605693  # Ry to eV
                    except (IndexError, ValueError):
                        pass
                
                elif 'highest occupied level' in line:
                    try:
                        ef_eV = float(line.split('=')[-1].split()[0])
                        results['fermi_energy_eV'] = ef_eV
                    except (IndexError, ValueError):
                        pass
            
            # Check convergence
            if 'convergence has been achieved' in content:
                results['converged'] = True
            else:
                results['converged'] = False
        
        except Exception as e:
            print(f"Warning: Could not parse QE output: {e}")
        
        return results


def create_slurm_script(
    calculation_dir: Path,
    job_name: str,
    qe_commands: List[str],
    nodes: int = 1,
    ntasks: int = 24,
    time_hours: int = 24,
    partition: str = 'normal'
) -> Path:
    """Create SLURM submission script for QE calculations."""
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --time={time_hours}:00:00
#SBATCH --partition={partition}
#SBATCH --output={calculation_dir}/slurm_%j.out
#SBATCH --error={calculation_dir}/slurm_%j.err

# Load modules (adjust for your cluster)
module load quantum-espresso
module load intel-mpi

# Set up environment
export OMP_NUM_THREADS=1
cd {calculation_dir}

# Run calculations
"""
    
    for cmd in qe_commands:
        script_content += f"mpirun -np $SLURM_NTASKS {cmd}\n"
    
    script_file = calculation_dir / 'submit.sh'
    with open(script_file, 'w') as f:
        f.write(script_content)
    
    # Make executable
    script_file.chmod(0o755)
    
    return script_file


def process_single_structure(
    structure_file: str,
    output_dir: Path,
    ecut: float = 60.0,
    kpoints: List[int] = [6, 6, 6],
    use_slurm: bool = False
) -> Dict[str, Any]:
    """Process a single structure with QE."""
    
    # Load structure
    structure = structure_from_file(structure_file)
    formula = structure.composition.reduced_formula
    
    # Set up calculation directory
    calc_dir = output_dir / f"qe_{formula}"
    calc_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize QE calculation
    qe_calc = QECalculation(calc_dir)
    
    results = {
        'formula': formula,
        'structure_file': str(structure_file),
        'calculation_dir': str(calc_dir)
    }
    
    if use_slurm:
        # Create SLURM script
        script_file = create_slurm_script(
            calc_dir,
            f"qe_{formula}",
            ['pw.x -input input/scf.in > output/scf.out']
        )
        
        # Submit job
        try:
            result = subprocess.run(['sbatch', str(script_file)], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                job_id = result.stdout.strip().split()[-1]
                results['slurm_job_id'] = job_id
                results['status'] = 'submitted'
            else:
                results['status'] = 'submission_failed'
                results['error'] = result.stderr
        except Exception as e:
            results['status'] = 'submission_failed' 
            results['error'] = str(e)
    
    else:
        # Run locally
        print(f"Running SCF calculation for {formula}...")
        success, scf_results = qe_calc.scf_calculation(structure, ecut, kpoints)
        
        if success:
            results.update(scf_results)
            results['status'] = 'completed'
            
            # Run relaxation if SCF succeeded
            print(f"Running relaxation for {formula}...")
            relax_success, relax_results, opt_structure = qe_calc.relax_calculation(
                structure, ecut, [max(1, k-2) for k in kpoints]  # Coarser k-mesh for relax
            )
            
            if relax_success and opt_structure:
                results.update({f"relax_{k}": v for k, v in relax_results.items()})
                results['relaxation_successful'] = True
                
                # Save optimized structure
                opt_file = calc_dir / 'optimized_structure.cif'
                structure_to_cif(opt_structure, opt_file)
                results['optimized_structure_file'] = str(opt_file)
            else:
                results['relaxation_successful'] = False
        
        else:
            results['status'] = 'failed'
            results.update(scf_results)
    
    return results


def batch_process_candidates(
    candidates_file: str,
    output_dir: Path,
    n_best: int = 50,
    ecut: float = 60.0,
    use_slurm: bool = False
) -> List[Dict[str, Any]]:
    """Process batch of candidates with QE."""
    
    # Load candidates
    df = pd.read_csv(candidates_file)
    
    # Sort by best metric (RBT score or xTB stability)
    if 'xtb_stability_score' in df.columns:
        df = df.sort_values('xtb_stability_score', ascending=False)
    elif 'best_rbt_score' in df.columns:
        df = df.sort_values('best_rbt_score', ascending=False)
    else:
        print("Warning: No suitable ranking column found")
    
    # Take top candidates
    top_candidates = df.head(n_best)
    
    print(f"Processing {len(top_candidates)} candidates with QE DFT...")
    
    all_results = []
    
    for i, candidate in tqdm(top_candidates.iterrows(), total=len(top_candidates), desc="QE calculations"):
        try:
            formula = candidate.get('reduced_formula', f'candidate_{i}')
            
            # Generate structure if not provided
            if 'optimized_structure_file' in candidate and pd.notna(candidate['optimized_structure_file']):
                structure_file = candidate['optimized_structure_file']
            else:
                # Create structure from composition
                composition = Composition(formula)
                structure = structure_from_composition(composition)
                
                # Save temporary structure file
                temp_dir = output_dir / 'temp_structures'
                temp_dir.mkdir(exist_ok=True)
                structure_file = temp_dir / f"{formula}.cif"
                structure_to_cif(structure, structure_file)
            
            # Determine k-points based on structure size
            try:
                structure = structure_from_file(structure_file)
                abc = structure.lattice.abc
                kpoints = [max(1, int(25 / a)) for a in abc]  # Adaptive k-mesh
            except:
                kpoints = [6, 6, 6]  # Default
            
            # Process structure
            results = process_single_structure(
                str(structure_file),
                output_dir,
                ecut=ecut,
                kpoints=kpoints,
                use_slurm=use_slurm
            )
            
            # Combine with candidate data
            combined_results = {**candidate.to_dict(), **results}
            all_results.append(combined_results)
            
        except Exception as e:
            print(f"Error processing candidate {i}: {e}")
            error_results = candidate.to_dict()
            error_results.update({
                'status': 'error',
                'error': str(e)
            })
            all_results.append(error_results)
    
    return all_results


def main():
    parser = argparse.ArgumentParser(
        description="High-accuracy DFT calculations with Quantum ESPRESSO"
    )
    
    # Input modes
    parser.add_argument(
        "input", nargs='?',
        help="Input structure file (CIF) or candidates CSV"
    )
    parser.add_argument(
        "--batch", type=str,
        help="Batch mode: process candidates from CSV file"
    )
    parser.add_argument(
        "--analyze", type=str,
        help="Analyze completed calculations in directory"
    )
    
    # Calculation parameters
    parser.add_argument(
        "--ecut", type=float, default=60.0,
        help="Plane wave cutoff in Ry"
    )
    parser.add_argument(
        "--kmesh", type=int, nargs=3, default=[6, 6, 6],
        help="k-point mesh (nk1 nk2 nk3)"
    )
    parser.add_argument(
        "--functional", choices=['PBE', 'PBEsol', 'LDA'], default='PBE',
        help="Exchange-correlation functional"
    )
    
    # Job control
    parser.add_argument(
        "--n_best", type=int, default=50,
        help="Number of best candidates to process (batch mode)"
    )
    parser.add_argument(
        "--slurm", action="store_true",
        help="Submit jobs to SLURM cluster"
    )
    parser.add_argument(
        "--nodes", type=int, default=1,
        help="Number of SLURM nodes"
    )
    parser.add_argument(
        "--ntasks", type=int, default=24,
        help="Number of SLURM tasks"
    )
    parser.add_argument(
        "--time", type=int, default=24,
        help="SLURM job time limit (hours)"
    )
    
    # Output
    parser.add_argument(
        "--out", type=str, default="results/qe_results",
        help="Output directory or file"
    )
    
    args = parser.parse_args()
    
    output_path = Path(args.out)
    
    if args.batch:
        # Batch processing mode
        output_path.mkdir(parents=True, exist_ok=True)
        
        results = batch_process_candidates(
            args.batch,
            output_path,
            n_best=args.n_best,
            ecut=args.ecut,
            use_slurm=args.slurm
        )
        
        # Save results
        results_file = output_path / 'qe_batch_results.csv'
        pd.DataFrame(results).to_csv(results_file, index=False)
        print(f"Saved batch results to {results_file}")
        
    elif args.analyze:
        # Analysis mode - collect results from completed calculations
        print(f"Analyzing completed calculations in {args.analyze}...")
        
        analysis_results = []
        calc_dirs = list(Path(args.analyze).glob("qe_*"))
        
        for calc_dir in calc_dirs:
            try:
                # Look for output files
                scf_out = calc_dir / 'output' / 'scf.out'
                relax_out = calc_dir / 'output' / 'relax.out'
                
                result = {'calculation_dir': str(calc_dir)}
                
                if scf_out.exists():
                    qe_calc = QECalculation(calc_dir)
                    scf_results = qe_calc._parse_pw_output(scf_out)
                    result.update({f"scf_{k}": v for k, v in scf_results.items()})
                
                if relax_out.exists():
                    relax_results = qe_calc._parse_pw_output(relax_out)
                    result.update({f"relax_{k}": v for k, v in relax_results.items()})
                
                analysis_results.append(result)
                
            except Exception as e:
                print(f"Error analyzing {calc_dir}: {e}")
        
        # Save analysis
        df_analysis = pd.DataFrame(analysis_results)
        df_analysis.to_csv(args.out, index=False)
        print(f"Saved analysis to {args.out}")
        
    elif args.input:
        # Single structure mode
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        results = process_single_structure(
            args.input,
            output_path.parent,
            ecut=args.ecut,
            kpoints=args.kmesh,
            use_slurm=args.slurm
        )
        
        print(f"Results: {results}")
        
        # Save results
        results_file = output_path.parent / f"qe_results_{Path(args.input).stem}.json"
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
    else:
        print("Must specify input file, --batch, or --analyze mode")


if __name__ == "__main__":
    main() 