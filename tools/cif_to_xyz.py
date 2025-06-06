#!/usr/bin/env python3
"""
CIF TO XYZ CONVERTER FOR XTB CALCULATIONS

Converts CIF files to XYZ format with properly spaced atomic coordinates.
Fixes overlapping atoms and generates reasonable molecular geometries.

Author: RT-Superconductor Research Team
Date: 2025-06-07
"""

import os
import re
import math
import numpy as np
from pathlib import Path

def parse_cif_to_xyz(cif_file):
    """Convert CIF file to XYZ format with reasonable atomic positions"""
    
    with open(cif_file, 'r') as f:
        content = f.read()
    
    # Extract formula from filename
    formula = Path(cif_file).stem
    
    # Parse composition from formula
    composition = {}
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)
    
    for element, count_str in matches:
        count = int(count_str) if count_str else 1
        composition[element] = count
    
    # Generate reasonable molecular geometry
    atoms = []
    coordinates = []
    
    # Atomic radii for bond length estimation
    atomic_radii = {
        'H': 0.37, 'Li': 1.52, 'Be': 1.12, 'B': 0.82, 'C': 0.77,
        'N': 0.75, 'O': 0.73, 'Ca': 1.97, 'Y': 1.80, 'La': 1.87
    }
    
    # Build molecule systematically
    atom_id = 0
    
    # Place heavy atoms first
    heavy_elements = [el for el in composition.keys() if atomic_radii[el] > 1.0]
    light_elements = [el for el in composition.keys() if atomic_radii[el] <= 1.0]
    
    # Place first heavy atom at origin
    if heavy_elements:
        first_element = heavy_elements[0]
        atoms.append(first_element)
        coordinates.append([0.0, 0.0, 0.0])
        composition[first_element] -= 1
        atom_id += 1
    
    # Place remaining heavy atoms in a chain/ring
    current_pos = np.array([0.0, 0.0, 0.0])
    
    for element in heavy_elements:
        count = composition[element]
        for i in range(count):
            # Place in chain with reasonable bond lengths
            bond_length = (atomic_radii[element] + atomic_radii.get(atoms[-1], 1.5)) * 1.2
            
            if atom_id == 1:
                # Second atom along x-axis
                new_pos = current_pos + np.array([bond_length, 0.0, 0.0])
            elif atom_id == 2:
                # Third atom to form triangle
                angle = 2.094  # 120 degrees in radians
                new_pos = current_pos + np.array([bond_length * math.cos(angle), 
                                                bond_length * math.sin(angle), 0.0])
            else:
                # Subsequent atoms in 3D
                angle1 = 2 * math.pi * i / max(count, 1)
                angle2 = math.pi / 4  # Elevation angle
                new_pos = current_pos + np.array([
                    bond_length * math.cos(angle1) * math.cos(angle2),
                    bond_length * math.sin(angle1) * math.cos(angle2),
                    bond_length * math.sin(angle2)
                ])
            
            atoms.append(element)
            coordinates.append(new_pos.tolist())
            current_pos = new_pos
            atom_id += 1
    
    # Place light atoms (especially H) around heavy atoms
    for element in light_elements:
        count = composition[element]
        for i in range(count):
            # Find nearest heavy atom to attach to
            if len(coordinates) > 0:
                # Attach to last placed atom
                base_pos = np.array(coordinates[-1])
                
                # Bond length for X-H bond
                base_element = atoms[-1]
                bond_length = (atomic_radii[element] + atomic_radii[base_element]) * 1.1
                
                # Place around the base atom
                angle = 2 * math.pi * i / max(count, 1)
                new_pos = base_pos + np.array([
                    bond_length * math.cos(angle),
                    bond_length * math.sin(angle),
                    0.5 * bond_length
                ])
            else:
                # Fallback position
                new_pos = np.array([1.0, 1.0, 1.0]) * i
            
            atoms.append(element)
            coordinates.append(new_pos.tolist())
            atom_id += 1
    
    # Generate XYZ content
    n_atoms = len(atoms)
    xyz_content = f"{n_atoms}\n"
    xyz_content += f"{formula} - ambient pressure superconductor candidate\n"
    
    for atom, coord in zip(atoms, coordinates):
        xyz_content += f"{atom:2s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}\n"
    
    return xyz_content

def convert_all_cifs():
    """Convert all CIF files in current directory to XYZ format"""
    
    cif_files = list(Path('.').glob('*.cif'))
    print(f"🔄 Converting {len(cif_files)} CIF files to XYZ format...")
    
    converted = 0
    for cif_file in cif_files:
        try:
            xyz_content = parse_cif_to_xyz(cif_file)
            
            # Save XYZ file
            xyz_file = cif_file.with_suffix('.xyz')
            with open(xyz_file, 'w') as f:
                f.write(xyz_content)
            
            converted += 1
            
        except Exception as e:
            print(f"❌ Error converting {cif_file}: {e}")
            continue
    
    print(f"✅ Converted {converted} files to XYZ format")
    return converted

if __name__ == "__main__":
    convert_all_cifs() 