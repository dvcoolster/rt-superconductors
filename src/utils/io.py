"""
I/O utilities for rt-superconductors: CIF ↔ POSCAR ↔ pw.in conversions.

Handles file format conversions needed for the quantum chemistry pipeline:
- CIF files (crystallographic information)
- VASP POSCAR files 
- Quantum ESPRESSO pw.in input files
- JSON serialization for intermediate results
"""

import json
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
import numpy as np
from pymatgen.core import Structure, Composition
from pymatgen.io.cif import CifWriter, CifParser
from pymatgen.io.vasp import Poscar
from pymatgen.io.pwscf import PWInput
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher
from ase import Atoms
from ase.io import write as ase_write, read as ase_read


def structure_from_file(filepath: Union[str, Path]) -> Structure:
    """
    Load structure from various file formats.
    
    Args:
        filepath: Path to structure file (CIF, POSCAR, etc.)
        
    Returns:
        Pymatgen Structure object
    """
    filepath = Path(filepath)
    
    if filepath.suffix.lower() in ['.cif']:
        parser = CifParser(str(filepath))
        return parser.get_structures()[0]
    elif filepath.name.upper() in ['POSCAR', 'CONTCAR'] or filepath.suffix.lower() == '.vasp':
        poscar = Poscar.from_file(str(filepath))
        return poscar.structure
    elif filepath.suffix.lower() in ['.xyz', '.pdb']:
        atoms = ase_read(str(filepath))
        return AseAtomsAdaptor.get_structure(atoms)
    else:
        # Try generic ASE reader
        try:
            atoms = ase_read(str(filepath))
            return AseAtomsAdaptor.get_structure(atoms)
        except Exception as e:
            raise ValueError(f"Unsupported file format for {filepath}: {e}")


def structure_to_cif(structure: Structure, filepath: Union[str, Path]) -> None:
    """
    Write structure to CIF file.
    
    Args:
        structure: Pymatgen Structure
        filepath: Output CIF file path
    """
    writer = CifWriter(structure)
    writer.write_file(str(filepath))


def structure_to_poscar(structure: Structure, filepath: Union[str, Path]) -> None:
    """
    Write structure to VASP POSCAR format.
    
    Args:
        structure: Pymatgen Structure
        filepath: Output POSCAR file path
    """
    poscar = Poscar(structure)
    poscar.write_file(str(filepath))


def structure_to_qe_input(
    structure: Structure,
    filepath: Union[str, Path],
    pseudopotentials: Optional[Dict[str, str]] = None,
    control: Optional[Dict[str, Any]] = None,
    system: Optional[Dict[str, Any]] = None,
    electrons: Optional[Dict[str, Any]] = None,
    kpoints: Optional[List[List[int]]] = None
) -> None:
    """
    Write Quantum ESPRESSO pw.x input file.
    
    Args:
        structure: Pymatgen Structure
        filepath: Output pw.in file path
        pseudopotentials: Dict mapping element to pseudopotential filename
        control: Control section parameters
        system: System section parameters  
        electrons: Electrons section parameters
        kpoints: K-points mesh as [nk1, nk2, nk3]
    """
    # Default pseudopotentials (assuming norm-conserving ONCV)
    if pseudopotentials is None:
        pseudopotentials = {el.symbol: f"{el.symbol}_ONCV_PBE-1.0.upf" 
                           for el in structure.composition.elements}
    
    # Default control parameters
    if control is None:
        control = {
            'calculation': 'scf',
            'restart_mode': 'from_scratch',
            'outdir': './tmp/',
            'pseudo_dir': './pseudos/',
            'prefix': 'pwscf',
            'verbosity': 'high'
        }
    
    # Default system parameters
    if system is None:
        system = {
            'ecutwfc': 60.0,  # Ry
            'ecutrho': 240.0,  # Ry
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.02
        }
    
    # Default electrons parameters
    if electrons is None:
        electrons = {
            'conv_thr': 1.0e-8,
            'mixing_beta': 0.3,
            'mixing_mode': 'plain'
        }
    
    # Default k-points
    if kpoints is None:
        # Simple heuristic based on lattice parameters
        abc = structure.lattice.abc
        kpoints = [max(1, int(20 / a)) for a in abc]
    
    # Create PWInput object
    pw_input = PWInput(
        structure=structure,
        pseudo=pseudopotentials,
        control=control,
        system=system,
        electrons=electrons,
        kpoints_grid=kpoints
    )
    
    # Write to file
    pw_input.write_file(str(filepath))


def load_supercon_data(filepath: Union[str, Path]) -> List[Dict[str, Any]]:
    """
    Load SuperCon database CSV file.
    
    Expected columns: formula, Tc, pressure, etc.
    
    Args:
        filepath: Path to SuperCon CSV file
        
    Returns:
        List of superconductor data dictionaries
    """
    import pandas as pd
    
    df = pd.read_csv(filepath)
    
    # Standard column name mappings
    column_map = {
        'Formula': 'formula',
        'Tc': 'tc_k',
        'Pressure': 'pressure_gpa',
        'Class': 'material_class'
    }
    
    # Rename columns if they exist
    for old_name, new_name in column_map.items():
        if old_name in df.columns:
            df = df.rename(columns={old_name: new_name})
    
    return df.to_dict('records')


def save_candidates_json(
    candidates: List[Dict[str, Any]], 
    filepath: Union[str, Path]
) -> None:
    """
    Save candidate structures/compositions to JSON.
    
    Args:
        candidates: List of candidate dictionaries
        filepath: Output JSON file path
    """
    # Convert Pymatgen objects to serializable format
    serializable_candidates = []
    
    for candidate in candidates:
        serializable = {}
        for key, value in candidate.items():
            if isinstance(value, (Structure, Composition)):
                serializable[key] = value.as_dict()
            elif isinstance(value, np.ndarray):
                serializable[key] = value.tolist()
            elif isinstance(value, (np.integer, np.floating)):
                serializable[key] = value.item()
            else:
                serializable[key] = value
        serializable_candidates.append(serializable)
    
    with open(filepath, 'w') as f:
        json.dump(serializable_candidates, f, indent=2)


def load_candidates_json(filepath: Union[str, Path]) -> List[Dict[str, Any]]:
    """
    Load candidate structures/compositions from JSON.
    
    Args:
        filepath: Input JSON file path
        
    Returns:
        List of candidate dictionaries with Pymatgen objects restored
    """
    with open(filepath, 'r') as f:
        data = json.load(f)
    
    candidates = []
    for item in data:
        candidate = {}
        for key, value in item.items():
            if key == 'structure' and isinstance(value, dict):
                candidate[key] = Structure.from_dict(value)
            elif key == 'composition' and isinstance(value, dict):
                candidate[key] = Composition.from_dict(value)
            else:
                candidate[key] = value
        candidates.append(candidate)
    
    return candidates


def structure_from_composition(
    composition: Composition,
    prototype: str = 'rocksalt',
    a: float = 5.0
) -> Structure:
    """
    Generate prototype structure from composition.
    
    Args:
        composition: Target composition
        prototype: Structure prototype ('rocksalt', 'perovskite', 'fluorite', etc.)
        a: Lattice parameter in Angstroms
        
    Returns:
        Prototype structure with target composition
    """
    elements = list(composition.elements)
    amounts = [composition[el] for el in elements]
    
    if prototype == 'rocksalt' and len(elements) == 2:
        # Simple cubic NaCl-type structure
        lattice = [[a, 0, 0], [0, a, 0], [0, 0, a]]
        species = [elements[0], elements[1]]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        
        return Structure(lattice, species, coords)
    
    elif prototype == 'perovskite' and len(elements) == 3:
        # Simple cubic ABX3 perovskite
        lattice = [[a, 0, 0], [0, a, 0], [0, 0, a]]
        species = [elements[0], elements[1], elements[2], elements[2], elements[2]]
        coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]]
        
        return Structure(lattice, species, coords)
    
    elif prototype == 'fluorite' and len(elements) == 2:
        # CaF2-type structure
        lattice = [[a, 0, 0], [0, a, 0], [0, 0, a]]
        species = [elements[0], elements[1], elements[1], elements[1], elements[1]]
        coords = [[0, 0, 0], [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], 
                  [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]]
        
        return Structure(lattice, species, coords)
    
    else:
        # Default to simple cubic with random placement
        n_atoms = sum(amounts)
        lattice = [[a, 0, 0], [0, a, 0], [0, 0, a]]
        
        species = []
        for el, amt in zip(elements, amounts):
            species.extend([el] * int(amt))
        
        # Generate random coordinates
        np.random.seed(42)  # Reproducible
        coords = np.random.random((len(species), 3))
        
        return Structure(lattice, species, coords)


def batch_convert_structures(
    input_dir: Union[str, Path],
    output_dir: Union[str, Path],
    input_format: str = 'cif',
    output_format: str = 'poscar'
) -> None:
    """
    Batch convert structures between formats.
    
    Args:
        input_dir: Directory containing input files
        output_dir: Directory for output files
        input_format: Input format ('cif', 'poscar', 'xyz')
        output_format: Output format ('cif', 'poscar', 'qe', 'xyz')
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find input files
    if input_format == 'cif':
        pattern = '*.cif'
    elif input_format == 'poscar':
        pattern = 'POSCAR*'
    elif input_format == 'xyz':
        pattern = '*.xyz'
    else:
        pattern = f'*.{input_format}'
    
    input_files = list(input_dir.glob(pattern))
    
    for input_file in input_files:
        try:
            # Load structure
            structure = structure_from_file(input_file)
            
            # Generate output filename
            stem = input_file.stem
            if output_format == 'cif':
                output_file = output_dir / f"{stem}.cif"
                structure_to_cif(structure, output_file)
            elif output_format == 'poscar':
                output_file = output_dir / f"POSCAR_{stem}"
                structure_to_poscar(structure, output_file)
            elif output_format == 'qe':
                output_file = output_dir / f"{stem}.in"
                structure_to_qe_input(structure, output_file)
            elif output_format == 'xyz':
                output_file = output_dir / f"{stem}.xyz"
                atoms = AseAtomsAdaptor.get_atoms(structure)
                ase_write(str(output_file), atoms)
            
            print(f"Converted {input_file.name} → {output_file.name}")
            
        except Exception as e:
            print(f"Failed to convert {input_file.name}: {e}")


# Utility for creating directory structure
def setup_calculation_directory(
    base_dir: Union[str, Path],
    structure_name: str,
    subdirs: Optional[List[str]] = None
) -> Path:
    """
    Set up directory structure for calculations.
    
    Args:
        base_dir: Base directory for calculations
        structure_name: Name of the structure/calculation
        subdirs: List of subdirectories to create
        
    Returns:
        Path to the calculation directory
    """
    if subdirs is None:
        subdirs = ['input', 'output', 'tmp', 'pseudos']
    
    calc_dir = Path(base_dir) / structure_name
    calc_dir.mkdir(parents=True, exist_ok=True)
    
    for subdir in subdirs:
        (calc_dir / subdir).mkdir(exist_ok=True)
    
    return calc_dir


# Export main functions
__all__ = [
    'structure_from_file',
    'structure_to_cif',
    'structure_to_poscar', 
    'structure_to_qe_input',
    'load_supercon_data',
    'save_candidates_json',
    'load_candidates_json',
    'structure_from_composition',
    'batch_convert_structures',
    'setup_calculation_directory'
] 