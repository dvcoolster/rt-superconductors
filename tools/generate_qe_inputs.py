#!/usr/bin/env python3
"""
QUANTUM ESPRESSO INPUT GENERATOR FOR AMBIENT-PRESSURE SUPERCONDUCTORS

Generates complete QE input files for SCF, NSCF, phonon, and EPW calculations
from xTB-optimized structures.

Author: RT-Superconductor Research Team
Date: 2025-06-07
"""

import os
import re
import argparse
import numpy as np
from pathlib import Path

class QEInputGenerator:
    """Generate Quantum ESPRESSO input files for superconductor calculations"""
    
    def __init__(self, ecut=80, ecut_rho=640):
        self.ecut = ecut  # Plane wave cutoff (Ry)
        self.ecut_rho = ecut_rho  # Charge density cutoff (Ry)
        
        # Pseudopotential files (assuming SSSP library)
        self.pseudopotentials = {
            'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
            'Li': 'Li.pbe-n-rrkjus_psl.1.0.0.UPF', 
            'Be': 'Be.pbe-n-rrkjus_psl.1.0.0.UPF',
            'B': 'B.pbe-n-rrkjus_psl.1.0.0.UPF',
            'C': 'C.pbe-n-rrkjus_psl.1.0.0.UPF',
            'N': 'N.pbe-n-rrkjus_psl.1.0.0.UPF',
            'O': 'O.pbe-n-rrkjus_psl.1.0.0.UPF',
            'Ca': 'Ca.pbe-nsp-rrkjus_psl.1.0.0.UPF',
            'Y': 'Y.pbe-spn-rrkjus_psl.1.0.0.UPF',
            'La': 'La.pbe-spfn-rrkjus_psl.1.0.0.UPF'
        }
        
        # Atomic masses (amu)
        self.atomic_masses = {
            'H': 1.008, 'Li': 6.94, 'Be': 9.01, 'B': 10.81, 'C': 12.01,
            'N': 14.01, 'O': 16.00, 'Ca': 40.08, 'Y': 88.91, 'La': 138.91
        }
    
    def parse_xyz_file(self, xyz_file):
        """Parse XYZ file to extract atomic positions"""
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        
        n_atoms = int(lines[0].strip())
        comment = lines[1].strip()
        
        atoms = []
        positions = []
        
        for i in range(2, 2 + n_atoms):
            parts = lines[i].strip().split()
            element = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            
            atoms.append(element)
            positions.append([x, y, z])
        
        return atoms, positions, comment
    
    def estimate_lattice_parameters(self, atoms, positions):
        """Estimate reasonable lattice parameters from molecular coordinates"""
        positions = np.array(positions)
        
        # Find molecular bounding box
        min_coords = np.min(positions, axis=0)
        max_coords = np.max(positions, axis=0)
        
        # Add padding for isolated molecule calculation
        padding = 8.0  # Angstrom
        
        a = max_coords[0] - min_coords[0] + 2 * padding
        b = max_coords[1] - min_coords[1] + 2 * padding  
        c = max_coords[2] - min_coords[2] + 2 * padding
        
        # Center molecule in box
        center_offset = np.array([a/2, b/2, c/2]) - (max_coords + min_coords) / 2
        centered_positions = positions + center_offset
        
        return a, b, c, centered_positions
    
    def generate_scf_input(self, atoms, positions, formula):
        """Generate SCF input file"""
        
        # Convert to lattice coordinates  
        a, b, c, cart_positions = self.estimate_lattice_parameters(atoms, positions)
        
        # Get unique species
        unique_species = list(set(atoms))
        n_species = len(unique_species)
        
        # Calculate total electrons for charge neutrality
        total_electrons = sum(self.get_valence_electrons(atom) for atom in atoms)
        
        scf_input = f"""&CONTROL
    calculation = 'scf',
    restart_mode = 'from_scratch',
    prefix = '{formula}',
    outdir = './tmp/',
    pseudo_dir = '../../../pseudos/',
    verbosity = 'high',
    tprnfor = .true.,
    tstress = .true.,
/

&SYSTEM
    ibrav = 1,
    A = {a:.6f},
    nat = {len(atoms)},
    ntyp = {n_species},
    ecutwfc = {self.ecut},
    ecutrho = {self.ecut_rho},
    occupations = 'smearing',
    smearing = 'marzari-vanderbilt',
    degauss = 0.02,
    input_dft = 'PBE',
/

&ELECTRONS
    diagonalization = 'david',
    mixing_mode = 'plain',
    mixing_beta = 0.7,
    conv_thr = 1.0d-8,
    electron_maxstep = 200,
/

ATOMIC_SPECIES
"""
        
        # Add atomic species
        for species in unique_species:
            mass = self.atomic_masses[species]
            pseudo = self.pseudopotentials[species]
            scf_input += f"{species:3s} {mass:8.3f} {pseudo}\n"
        
        scf_input += f"\nATOMIC_POSITIONS angstrom\n"
        
        # Add atomic positions
        for atom, pos in zip(atoms, cart_positions):
            scf_input += f"{atom:3s} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n"
        
        # K-points for isolated molecule
        scf_input += f"\nK_POINTS gamma\n"
        
        return scf_input
    
    def generate_nscf_input(self, atoms, positions, formula):
        """Generate NSCF input file for denser k-point grid"""
        
        scf_input = self.generate_scf_input(atoms, positions, formula)
        
        # Modify for NSCF
        nscf_input = scf_input.replace("calculation = 'scf'", "calculation = 'nscf'")
        nscf_input = nscf_input.replace("K_POINTS gamma", "K_POINTS automatic\n4 4 4 0 0 0")
        
        # Add more bands for unoccupied states
        total_electrons = sum(self.get_valence_electrons(atom) for atom in atoms)
        nbnd = max(20, int(total_electrons * 1.5))
        
        # Insert nbnd after degauss line
        nscf_input = nscf_input.replace(
            "degauss = 0.02,", 
            f"degauss = 0.02,\n    nbnd = {nbnd},"
        )
        
        return nscf_input
    
    def generate_ph_gamma_input(self, formula):
        """Generate phonon input for Gamma-point check"""
        
        ph_input = f"""Gamma-point phonon calculation for {formula}
&INPUTPH
    tr2_ph = 1.0d-14,
    prefix = '{formula}',
    outdir = './tmp/',
    alpha_mix(1) = 0.7,
    fildyn = '{formula}.dyn',
/
0.0 0.0 0.0
"""
        return ph_input
    
    def generate_ph_grid_input(self, formula):
        """Generate phonon input for 4x4x4 grid"""
        
        ph_input = f"""Full phonon calculation for {formula}
&INPUTPH
    tr2_ph = 1.0d-14,
    prefix = '{formula}',
    outdir = './tmp/',
    ldisp = .true.,
    nq1 = 4, nq2 = 4, nq3 = 4,
    alpha_mix(1) = 0.7,
    fildyn = '{formula}.dyn',
/
"""
        return ph_input
    
    def generate_q2r_input(self, formula):
        """Generate q2r input file"""
        
        q2r_input = f"""&INPUT
    fildyn = '{formula}.dyn',
    flfrc = '{formula}.fc',
/
"""
        return q2r_input
    
    def generate_matdyn_input(self, formula):
        """Generate matdyn input file"""
        
        matdyn_input = f"""&INPUT
    asr = 'crystal',
    flfrc = '{formula}.fc',
    flfrq = '{formula}.freq',
    q_in_band_form = .true.,
/
4
0.0 0.0 0.0 100
0.5 0.0 0.0 100  
0.5 0.5 0.0 100
0.0 0.0 0.0 1
"""
        return matdyn_input
    
    def generate_epw_input(self, atoms, formula):
        """Generate EPW input file"""
        
        # Count electrons for Fermi level
        total_electrons = sum(self.get_valence_electrons(atom) for atom in atoms)
        
        epw_input = f"""--
&INPUTEPW
    prefix = '{formula}',
    outdir = './tmp/',
    
    elph = .true.,
    kmaps = .false.,
    epbwrite = .true.,
    epbread = .false.,
    
    epwwrite = .true.,
    epwread = .false.,
    
    wannierize = .true.,
    num_iter = 1000,
    dis_win_max = 15.0,
    dis_win_min = -3.0,
    dis_froz_min = -3.0,
    dis_froz_max = 10.0,
    
    proj(1) = 's',
    proj(2) = 'p',
    
    iverbosity = 2,
    
    eps_acustic = 2.0,
    ephwrite = .true.,
    
    fsthick = 0.4,
    degaussw = 0.05,
    
    nkf1 = 8, nkf2 = 8, nkf3 = 8,
    nqf1 = 4, nqf2 = 4, nqf3 = 4,
    
    nk1 = 4, nk2 = 4, nk3 = 4,
    nq1 = 4, nq2 = 4, nq3 = 4,
/
"""
        return epw_input
    
    def get_valence_electrons(self, element):
        """Get number of valence electrons for pseudopotential"""
        valence = {
            'H': 1, 'Li': 1, 'Be': 2, 'B': 3, 'C': 4,
            'N': 5, 'O': 6, 'Ca': 10, 'Y': 11, 'La': 11
        }
        return valence.get(element, 4)
    
    def create_directory_structure(self, output_dir, formula):
        """Create directory structure for QE calculations"""
        
        calc_dir = Path(output_dir) / formula
        calc_dir.mkdir(parents=True, exist_ok=True)
        
        # Create subdirectories
        (calc_dir / 'tmp').mkdir(exist_ok=True)
        
        return calc_dir
    
    def generate_all_inputs(self, xyz_file, output_dir):
        """Generate all QE input files for a structure"""
        
        # Parse XYZ file
        atoms, positions, comment = self.parse_xyz_file(xyz_file)
        formula = Path(xyz_file).stem
        
        # Create directory structure
        calc_dir = self.create_directory_structure(output_dir, formula)
        
        print(f"  📁 Creating inputs for {formula}")
        
        # Generate all input files
        inputs = {
            'scf.in': self.generate_scf_input(atoms, positions, formula),
            'nscf.in': self.generate_nscf_input(atoms, positions, formula),
            'gamma.in': self.generate_ph_gamma_input(formula),
            'ph.in': self.generate_ph_grid_input(formula),
            'q2r.in': self.generate_q2r_input(formula),
            'matdyn.in': self.generate_matdyn_input(formula),
            'epw.in': self.generate_epw_input(atoms, formula)
        }
        
        # Write input files
        for filename, content in inputs.items():
            with open(calc_dir / filename, 'w') as f:
                f.write(content)
        
        # Copy optimized coordinates
        optimized_xyz = calc_dir / f"{formula}_opt.xyz"
        with open(optimized_xyz, 'w') as f:
            f.write(f"{len(atoms)}\n")
            f.write(f"{comment}\n")
            for atom, pos in zip(atoms, positions):
                f.write(f"{atom:3s} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
        
        return calc_dir

def main():
    parser = argparse.ArgumentParser(description='Generate QE input files for superconductor candidates')
    parser.add_argument('--list', type=str, required=True,
                       help='File containing list of survivor structures')
    parser.add_argument('--template', type=str, default='scf',
                       help='Template type (scf, nscf, ph, epw)')
    parser.add_argument('--ecut', type=int, default=80,
                       help='Plane wave cutoff (Ry)')
    parser.add_argument('--ecut-rho', dest='ecut_rho', type=int, default=640,
                       help='Charge density cutoff (Ry)')
    parser.add_argument('--out', type=str, required=True,
                       help='Output directory for QE calculations')
    
    args = parser.parse_args()
    
    print("⚛️  QUANTUM ESPRESSO INPUT GENERATOR")
    print("=" * 50)
    print(f"📋 Survivors list: {args.list}")
    print(f"🔧 Cutoffs: ecutwfc={args.ecut} Ry, ecutrho={args.ecut_rho} Ry")
    print(f"📂 Output directory: {args.out}")
    print("")
    
    # Read survivors list
    with open(args.list, 'r') as f:
        survivors = [line.strip() for line in f if line.strip()]
    
    print(f"🎯 Processing {len(survivors)} stable structures...")
    
    # Initialize generator
    generator = QEInputGenerator(ecut=args.ecut, ecut_rho=args.ecut_rho)
    
    # Create output directory
    output_dir = Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each survivor
    generated = 0
    for survivor in survivors:
        try:
            # Look for corresponding XYZ file
            xyz_file = f"{survivor}.xyz"
            if not Path(xyz_file).exists():
                print(f"  ❌ XYZ file not found: {xyz_file}")
                continue
            
            # Generate QE inputs
            calc_dir = generator.generate_all_inputs(xyz_file, output_dir)
            generated += 1
            
        except Exception as e:
            print(f"  ❌ Error processing {survivor}: {e}")
            continue
    
    print("")
    print("✅ QE INPUT GENERATION COMPLETE!")
    print("=" * 40)
    print(f"📊 Structures processed: {generated}")
    print(f"📂 Output directory: {output_dir}")
    print("")
    print("🚀 Ready for HPC submission!")
    print("   Next: Submit SLURM jobs for SCF+NSCF+Phonon calculations")

if __name__ == "__main__":
    main() 