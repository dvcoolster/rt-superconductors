#!/usr/bin/env python
"""
Fast geometry optimization using xTB (extended tight-binding) for superconductor screening.

This script performs semi-empirical geometry optimizations on candidate structures
to filter out unstable configurations before expensive DFT calculations. Uses the
GFN2-xTB method for speed and reasonable accuracy.

Usage:
    python 03_xtb_relax.py --input candidates.csv --n_candidates 100 --out results/xtb_optimized.csv

The xTB method provides:
- Fast geometry optimization (~seconds per structure)
- Reasonable energetics for stability assessment
- Vibrational frequency analysis
- Electronic structure properties

Reference: Bannwarth, C. et al. J. Chem. Theory Comput. 2019, 15, 1652–1671
"""

import argparse
import json
import os
import tempfile
import subprocess
import shutil
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import numpy as np
import pandas as pd
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# ASE for structure manipulation
from ase import Atoms
from ase.io import write as ase_write, read as ase_read
from ase.optimize import BFGS
from ase.constraints import FixAtoms

# xTB integration
try:
    from xtb.ase.calculator import XTB
    XTB_AVAILABLE = True
except ImportError:
    print("Warning: xtb-python not available. Will use external xTB binary.")
    XTB_AVAILABLE = False

# Pymatgen for structure handling
from pymatgen.core import Structure, Composition
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer

from utils.ledger import rbt_superconductor_score
from utils.io import structure_from_file, structure_to_cif, structure_from_composition


class XTBCalculator:
    """Wrapper for xTB calculations using external binary or Python interface."""
    
    def __init__(self, method: str = 'GFN2-xTB', charge: int = 0, uhf: int = 0):
        self.method = method
        self.charge = charge
        self.uhf = uhf  # Number of unpaired electrons
        self.use_python_xtb = XTB_AVAILABLE
        
        if not self.use_python_xtb:
            # Check for external xTB binary
            self.xtb_command = self._find_xtb_binary()
            if not self.xtb_command:
                raise RuntimeError("Neither xtb-python nor external xTB binary found")
    
    def _find_xtb_binary(self) -> Optional[str]:
        """Find xTB binary in PATH."""
        for cmd in ['xtb', 'xtb.exe']:
            if shutil.which(cmd):
                return cmd
        return None
    
    def optimize_structure(self, structure: Structure, max_steps: int = 200) -> Tuple[Structure, Dict[str, Any]]:
        """
        Optimize structure geometry using xTB.
        
        Args:
            structure: Input structure
            max_steps: Maximum optimization steps
            
        Returns:
            Tuple of (optimized_structure, results_dict)
        """
        if self.use_python_xtb:
            return self._optimize_with_python_xtb(structure, max_steps)
        else:
            return self._optimize_with_external_xtb(structure, max_steps)
    
    def _optimize_with_python_xtb(self, structure: Structure, max_steps: int) -> Tuple[Structure, Dict[str, Any]]:
        """Optimize using xtb-python interface."""
        # Convert to ASE atoms
        atoms = AseAtomsAdaptor.get_atoms(structure)
        
        # Set up xTB calculator
        calc = XTB(method=self.method, charge=self.charge, uhf=self.uhf)
        atoms.set_calculator(calc)
        
        # Initial energy
        initial_energy = atoms.get_potential_energy()
        
        # Optimize geometry
        optimizer = BFGS(atoms, logfile=None)
        try:
            optimizer.run(fmax=0.05, steps=max_steps)  # 0.05 eV/Å convergence
            converged = optimizer.converged()
        except Exception as e:
            print(f"Optimization failed: {e}")
            converged = False
        
        # Final energy and properties
        final_energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        max_force = np.max(np.linalg.norm(forces, axis=1))
        
        # Convert back to pymatgen structure
        optimized_structure = AseAtomsAdaptor.get_structure(atoms)
        
        results = {
            'converged': converged,
            'initial_energy_eV': initial_energy,
            'final_energy_eV': final_energy,
            'energy_change_eV': final_energy - initial_energy,
            'max_force_eV_A': max_force,
            'n_atoms': len(atoms)
        }
        
        return optimized_structure, results
    
    def _optimize_with_external_xtb(self, structure: Structure, max_steps: int) -> Tuple[Structure, Dict[str, Any]]:
        """Optimize using external xTB binary."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Write input structure
            input_file = tmpdir / "input.xyz"
            atoms = AseAtomsAdaptor.get_atoms(structure)
            ase_write(str(input_file), atoms)
            
            # Run xTB optimization
            output_file = tmpdir / "xtbopt.xyz"
            log_file = tmpdir / "xtb.out"
            
            cmd = [
                self.xtb_command,
                str(input_file),
                f"--{self.method.lower()}",
                "--opt",
                f"--chrg {self.charge}",
                f"--uhf {self.uhf}",
                f"--cycles {max_steps}",
                "--parallel", "1"
            ]
            
            try:
                with open(log_file, 'w') as f:
                    subprocess.run(cmd, cwd=tmpdir, stdout=f, stderr=subprocess.STDOUT, 
                                 check=True, timeout=300)  # 5 minute timeout
                
                # Read optimized structure
                if output_file.exists():
                    opt_atoms = ase_read(str(output_file))
                    optimized_structure = AseAtomsAdaptor.get_structure(opt_atoms)
                else:
                    # Fallback to input if optimization failed
                    optimized_structure = structure
                
                # Parse results from log
                results = self._parse_xtb_output(log_file)
                
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
                print(f"xTB calculation failed: {e}")
                optimized_structure = structure
                results = {
                    'converged': False,
                    'error': str(e)
                }
        
        return optimized_structure, results
    
    def _parse_xtb_output(self, log_file: Path) -> Dict[str, Any]:
        """Parse xTB output file for results."""
        results = {
            'converged': False,
            'final_energy_eV': None,
            'max_force_eV_A': None,
            'homo_eV': None,
            'lumo_eV': None,
            'gap_eV': None
        }
        
        try:
            with open(log_file, 'r') as f:
                content = f.read()
            
            # Look for convergence
            if "GEOMETRY OPTIMIZATION CONVERGED" in content:
                results['converged'] = True
            
            # Extract final energy (in Hartree, convert to eV)
            for line in content.split('\n'):
                if "TOTAL ENERGY" in line:
                    try:
                        energy_hartree = float(line.split()[-2])
                        results['final_energy_eV'] = energy_hartree * 27.211386  # Hartree to eV
                    except (IndexError, ValueError):
                        pass
                
                elif "HOMO-LUMO GAP" in line:
                    try:
                        gap_eV = float(line.split()[-2])
                        results['gap_eV'] = gap_eV
                    except (IndexError, ValueError):
                        pass
        
        except Exception as e:
            print(f"Warning: Could not parse xTB output: {e}")
        
        return results


def assess_structure_stability(structure: Structure, xtb_results: Dict[str, Any]) -> Dict[str, float]:
    """
    Assess structural stability from xTB results and geometric criteria.
    
    Args:
        structure: Optimized structure
        xtb_results: Results from xTB calculation
        
    Returns:
        Dictionary of stability metrics
    """
    stability = {}
    
    # Convergence check
    stability['optimization_converged'] = float(xtb_results.get('converged', False))
    
    # Energy per atom (rough stability indicator)
    if xtb_results.get('final_energy_eV') is not None:
        energy_per_atom = xtb_results['final_energy_eV'] / structure.num_sites
        stability['energy_per_atom_eV'] = energy_per_atom
        
        # Heuristic: stable materials typically have -5 to -15 eV/atom
        if -20 < energy_per_atom < -2:
            stability['energy_reasonable'] = 1.0
        else:
            stability['energy_reasonable'] = 0.0
    else:
        stability['energy_per_atom_eV'] = 0.0
        stability['energy_reasonable'] = 0.0
    
    # Electronic gap check
    gap_eV = xtb_results.get('gap_eV', 0)
    stability['electronic_gap_eV'] = gap_eV
    
    # For superconductors, we often want small gaps or metallic behavior
    if gap_eV < 0.5:  # Metallic or small gap
        stability['gap_favorable'] = 1.0
    elif gap_eV < 2.0:  # Small semiconductor gap
        stability['gap_favorable'] = 0.7
    else:  # Large gap (less favorable)
        stability['gap_favorable'] = 0.3
    
    # Structural quality checks
    try:
        # Check for reasonable bond distances
        distances = structure.distance_matrix
        min_distance = np.min(distances[distances > 0])
        
        if min_distance > 1.0:  # Reasonable minimum bond distance
            stability['bond_distances_ok'] = 1.0
        else:
            stability['bond_distances_ok'] = 0.0
            
        stability['min_bond_distance_A'] = min_distance
        
    except Exception:
        stability['bond_distances_ok'] = 0.5
        stability['min_bond_distance_A'] = 2.0
    
    # Density check
    try:
        density = structure.density
        if 1.0 < density < 25.0:  # Reasonable density range
            stability['density_reasonable'] = 1.0
        else:
            stability['density_reasonable'] = 0.0
        stability['density_g_cm3'] = density
    except Exception:
        stability['density_reasonable'] = 0.5
        stability['density_g_cm3'] = 5.0
    
    # Combined stability score
    stability['overall_stability'] = np.mean([
        stability['optimization_converged'],
        stability['energy_reasonable'], 
        stability['gap_favorable'],
        stability['bond_distances_ok'],
        stability['density_reasonable']
    ])
    
    return stability


def process_candidate_batch(
    candidates: List[Dict[str, Any]],
    calculator: XTBCalculator,
    output_dir: Path,
    save_structures: bool = True
) -> List[Dict[str, Any]]:
    """
    Process a batch of candidates with xTB optimization.
    
    Args:
        candidates: List of candidate dictionaries
        calculator: XTB calculator instance
        output_dir: Directory to save results
        save_structures: Whether to save optimized structures
        
    Returns:
        List of processed candidates with xTB results
    """
    processed = []
    
    for i, candidate in enumerate(tqdm(candidates, desc="xTB optimization")):
        try:
            # Get structure
            formula = candidate.get('reduced_formula', candidate.get('formula', ''))
            
            if 'structure' in candidate:
                # Use provided structure
                if isinstance(candidate['structure'], dict):
                    structure = Structure.from_dict(candidate['structure'])
                else:
                    structure = candidate['structure']
            else:
                # Generate prototype structure
                composition = Composition(formula)
                prototype = candidate.get('best_prototype', 'rocksalt')
                structure = structure_from_composition(composition, prototype=prototype)
            
            # Run xTB optimization
            opt_structure, xtb_results = calculator.optimize_structure(structure)
            
            # Assess stability
            stability = assess_structure_stability(opt_structure, xtb_results)
            
            # Update candidate with results
            enhanced_candidate = {
                **candidate,
                'xtb_converged': xtb_results.get('converged', False),
                'xtb_energy_eV': xtb_results.get('final_energy_eV', 0),
                'xtb_energy_per_atom': stability.get('energy_per_atom_eV', 0),
                'xtb_gap_eV': xtb_results.get('gap_eV', 0),
                'xtb_stability_score': stability.get('overall_stability', 0),
                'min_bond_distance': stability.get('min_bond_distance_A', 0),
                'optimized_density': stability.get('density_g_cm3', 0),
                'optimization_successful': True
            }
            
            # Save optimized structure if requested
            if save_structures and xtb_results.get('converged', False):
                struct_file = output_dir / f"structures/opt_{formula}_{i:04d}.cif"
                struct_file.parent.mkdir(exist_ok=True)
                try:
                    structure_to_cif(opt_structure, struct_file)
                    enhanced_candidate['optimized_structure_file'] = str(struct_file)
                except Exception as e:
                    print(f"Warning: Could not save structure for {formula}: {e}")
            
            processed.append(enhanced_candidate)
            
        except Exception as e:
            print(f"Error processing {candidate.get('formula', 'unknown')}: {e}")
            # Add failed candidate with error info
            failed_candidate = {
                **candidate,
                'xtb_converged': False,
                'xtb_stability_score': 0.0,
                'optimization_successful': False,
                'error': str(e)
            }
            processed.append(failed_candidate)
    
    return processed


def filter_by_xtb_criteria(
    candidates: List[Dict[str, Any]],
    min_stability: float = 0.6,
    require_convergence: bool = True,
    max_gap_eV: float = 3.0
) -> List[Dict[str, Any]]:
    """
    Filter candidates based on xTB results.
    
    Args:
        candidates: List of candidates with xTB results
        min_stability: Minimum stability score required
        require_convergence: Whether to require optimization convergence
        max_gap_eV: Maximum allowed electronic gap
        
    Returns:
        Filtered list of candidates
    """
    filtered = []
    
    for candidate in candidates:
        # Check convergence
        if require_convergence and not candidate.get('xtb_converged', False):
            continue
        
        # Check stability score
        stability = candidate.get('xtb_stability_score', 0)
        if stability < min_stability:
            continue
        
        # Check electronic gap
        gap = candidate.get('xtb_gap_eV', float('inf'))
        if gap > max_gap_eV:
            continue
        
        filtered.append(candidate)
    
    print(f"xTB filtering: {len(filtered)}/{len(candidates)} candidates passed")
    return filtered


def main():
    parser = argparse.ArgumentParser(
        description="Fast xTB geometry optimization for superconductor screening"
    )
    parser.add_argument(
        "--input", required=True,
        help="Input CSV file with candidates"
    )
    parser.add_argument(
        "--n_candidates", type=int, default=100,
        help="Number of top candidates to process"
    )
    parser.add_argument(
        "--method", choices=['GFN2-xTB', 'GFN1-xTB', 'GFN0-xTB'], default='GFN2-xTB',
        help="xTB method to use"
    )
    parser.add_argument(
        "--charge", type=int, default=0,
        help="Total charge of the system"
    )
    parser.add_argument(
        "--max_steps", type=int, default=200,
        help="Maximum optimization steps"
    )
    parser.add_argument(
        "--min_stability", type=float, default=0.6,
        help="Minimum stability score for filtering"
    )
    parser.add_argument(
        "--max_gap", type=float, default=3.0,
        help="Maximum allowed electronic gap (eV)"
    )
    parser.add_argument(
        "--out", required=True,
        help="Output CSV file path"
    )
    parser.add_argument(
        "--save_structures", action="store_true",
        help="Save optimized structures as CIF files"
    )
    parser.add_argument(
        "--parallel", type=int, default=1,
        help="Number of parallel processes (not implemented yet)"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    output_path = Path(args.out)
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load candidates
    print(f"Loading candidates from {args.input}...")
    candidates_df = pd.read_csv(args.input)
    candidates = candidates_df.head(args.n_candidates).to_dict('records')
    print(f"Processing {len(candidates)} candidates with xTB")
    
    # Initialize xTB calculator
    calculator = XTBCalculator(
        method=args.method,
        charge=args.charge
    )
    
    # Process candidates
    processed_candidates = process_candidate_batch(
        candidates,
        calculator,
        output_dir,
        save_structures=args.save_structures
    )
    
    # Filter by xTB criteria
    filtered_candidates = filter_by_xtb_criteria(
        processed_candidates,
        min_stability=args.min_stability,
        require_convergence=True,
        max_gap_eV=args.max_gap
    )
    
    # Sort by stability score
    filtered_candidates.sort(
        key=lambda x: x.get('xtb_stability_score', 0), 
        reverse=True
    )
    
    # Save results
    df_out = pd.DataFrame(filtered_candidates)
    
    # Select key columns for output
    output_columns = [
        'reduced_formula', 'xtb_stability_score', 'xtb_converged',
        'xtb_energy_per_atom', 'xtb_gap_eV', 'min_bond_distance',
        'optimized_density', 'best_rbt_score', 'best_tc_estimate'
    ]
    
    available_columns = [col for col in output_columns if col in df_out.columns]
    df_output = df_out[available_columns]
    
    df_output.to_csv(output_path, index=False, float_format='%.4f')
    
    print(f"\nSaved {len(filtered_candidates)} optimized candidates to {output_path}")
    
    # Print summary statistics
    print("\nxTB Optimization Summary:")
    print(f"Candidates processed: {len(processed_candidates)}")
    print(f"Optimizations converged: {sum(1 for c in processed_candidates if c.get('xtb_converged', False))}")
    print(f"Passed stability filter: {len(filtered_candidates)}")
    
    if filtered_candidates:
        stability_scores = [c.get('xtb_stability_score', 0) for c in filtered_candidates]
        gap_values = [c.get('xtb_gap_eV', 0) for c in filtered_candidates]
        
        print(f"Stability score range: {np.min(stability_scores):.3f} - {np.max(stability_scores):.3f}")
        print(f"Electronic gap range: {np.min(gap_values):.3f} - {np.max(gap_values):.3f} eV")
        
        # Top 10 candidates
        print("\nTop 10 candidates by stability:")
        print("Rank | Formula | Stability | Gap (eV) | Energy/atom | RBT Score")
        print("-" * 70)
        
        for i, candidate in enumerate(filtered_candidates[:10]):
            print(f"{i+1:4d} | {candidate.get('reduced_formula', 'N/A'):12s} | "
                  f"{candidate.get('xtb_stability_score', 0):8.3f} | "
                  f"{candidate.get('xtb_gap_eV', 0):7.3f} | "
                  f"{candidate.get('xtb_energy_per_atom', 0):10.2f} | "
                  f"{candidate.get('best_rbt_score', 0):8.4f}")


if __name__ == "__main__":
    main() 