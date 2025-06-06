#!/usr/bin/env python3
"""
CIF EXPORT TOOL FOR AMBIENT-PRESSURE SUPERCONDUCTOR CANDIDATES

Converts chemical pressure analysis results to CIF format for xTB calculations.
Generates reasonable crystal structures from composition formulas.

Author: RT-Superconductor Research Team
Date: 2025-06-07
"""

import json
import argparse
import os
import numpy as np
from pathlib import Path
import math

class CIFGenerator:
    """Generate CIF files from composition formulas"""
    
    def __init__(self):
        # Atomic radii (Å) for structure generation
        self.atomic_radii = {
            'H': 0.37, 'Li': 1.52, 'Be': 1.12, 'B': 0.82, 'C': 0.77,
            'N': 0.75, 'O': 0.73, 'Ca': 1.97, 'Y': 1.80, 'La': 1.87
        }
        
        # Common crystal systems for each material category
        self.structure_templates = {
            'Hydride': {
                'space_group': 'Fm-3m',  # Face-centered cubic
                'lattice_type': 'cubic',
                'coordination': 8
            },
            'Electride': {
                'space_group': 'I4/mmm', # Body-centered tetragonal
                'lattice_type': 'tetragonal',
                'coordination': 6
            },
            'Clathrate': {
                'space_group': 'Pm-3n',  # Clathrate structure
                'lattice_type': 'cubic',
                'coordination': 20
            },
            'Li-B-C-N': {
                'space_group': 'P6/mmm', # Hexagonal layered
                'lattice_type': 'hexagonal',
                'coordination': 3
            },
            'Other': {
                'space_group': 'P1',     # Triclinic (general)
                'lattice_type': 'triclinic',
                'coordination': 6
            }
        }
    
    def estimate_lattice_parameters(self, composition: dict, category: str) -> dict:
        """Estimate lattice parameters from composition and category"""
        total_atoms = sum(composition.values())
        
        # Estimate volume per atom
        volume_per_atom = 0
        for element, count in composition.items():
            radius = self.atomic_radii[element]
            atomic_volume = (4/3) * math.pi * radius**3
            volume_per_atom += count * atomic_volume / total_atoms
        
        # Apply packing efficiency
        packing_efficiency = 0.74  # Close packing
        cell_volume = total_atoms * volume_per_atom / packing_efficiency
        
        # Estimate lattice parameters based on crystal system
        template = self.structure_templates[category]
        
        if template['lattice_type'] == 'cubic':
            a = b = c = cell_volume**(1/3)
            alpha = beta = gamma = 90.0
        elif template['lattice_type'] == 'tetragonal':
            a = b = (cell_volume / 1.2)**(1/3)  # Slightly compressed c
            c = 1.2 * a
            alpha = beta = gamma = 90.0
        elif template['lattice_type'] == 'hexagonal':
            a = b = (cell_volume / 1.633)**(1/3)  # c/a = 1.633 for ideal hex
            c = 1.633 * a
            alpha = beta = 90.0
            gamma = 120.0
        else:  # triclinic
            a = cell_volume**(1/3)
            b = a * 1.1
            c = a * 0.9
            alpha = 85.0
            beta = 95.0
            gamma = 100.0
        
        return {
            'a': a, 'b': b, 'c': c,
            'alpha': alpha, 'beta': beta, 'gamma': gamma,
            'space_group': template['space_group']
        }
    
    def generate_atomic_positions(self, composition: dict, lattice_params: dict) -> list:
        """Generate reasonable atomic positions"""
        positions = []
        
        # Simple position generation based on space group
        space_group = lattice_params['space_group']
        
        if space_group in ['Fm-3m', 'Pm-3n']:  # Cubic structures
            # Place atoms at special positions
            elements = list(composition.keys())
            
            # Heavy atoms at origin and corners
            heavy_elements = [el for el in elements if self.atomic_radii[el] > 1.0]
            light_elements = [el for el in elements if self.atomic_radii[el] <= 1.0]
            
            atom_id = 1
            for element in heavy_elements:
                count = composition[element]
                for i in range(count):
                    if i == 0:
                        # Origin
                        positions.append({
                            'label': f'{element}{atom_id}',
                            'element': element,
                            'x': 0.0, 'y': 0.0, 'z': 0.0,
                            'occupancy': 1.0
                        })
                    elif i == 1 and count > 1:
                        # Face center
                        positions.append({
                            'label': f'{element}{atom_id}',
                            'element': element,
                            'x': 0.5, 'y': 0.5, 'z': 0.0,
                            'occupancy': 1.0
                        })
                    atom_id += 1
            
            # Light atoms at interstitial sites
            for element in light_elements:
                count = composition[element]
                for i in range(count):
                    frac = i / max(count, 1)
                    positions.append({
                        'label': f'{element}{atom_id}',
                        'element': element,
                        'x': 0.25 + frac * 0.5,
                        'y': 0.25,
                        'z': 0.25,
                        'occupancy': 1.0
                    })
                    atom_id += 1
        
        elif space_group == 'I4/mmm':  # Tetragonal (electrides)
            # Layered structure with cavity spaces
            elements = list(composition.keys())
            z_layers = [0.0, 0.25, 0.5, 0.75]
            
            atom_id = 1
            layer_idx = 0
            for element in elements:
                count = composition[element]
                for i in range(count):
                    z = z_layers[layer_idx % len(z_layers)]
                    x = 0.0 if element in ['Ca', 'Y', 'La'] else 0.5
                    y = 0.0 if element in ['Ca', 'Y', 'La'] else 0.5
                    
                    positions.append({
                        'label': f'{element}{atom_id}',
                        'element': element,
                        'x': x, 'y': y, 'z': z,
                        'occupancy': 1.0
                    })
                    atom_id += 1
                    layer_idx += 1
        
        elif space_group == 'P6/mmm':  # Hexagonal layers
            # Layered structure for Li-B-C-N networks
            elements = list(composition.keys())
            
            atom_id = 1
            for element in elements:
                count = composition[element]
                for i in range(count):
                    # Distribute in hexagonal layers
                    angle = 2 * math.pi * i / count
                    x = 0.333 + 0.2 * math.cos(angle)
                    y = 0.667 + 0.2 * math.sin(angle)
                    z = 0.0 if element in ['B', 'C'] else 0.5
                    
                    positions.append({
                        'label': f'{element}{atom_id}',
                        'element': element,
                        'x': x, 'y': y, 'z': z,
                        'occupancy': 1.0
                    })
                    atom_id += 1
        
        else:  # General positions for triclinic
            elements = list(composition.keys())
            atom_id = 1
            total_atoms = sum(composition.values())
            
            for element in elements:
                count = composition[element]
                for i in range(count):
                    # Random but reasonable positions
                    x = 0.1 + 0.8 * (atom_id - 1) / total_atoms
                    y = 0.2 + 0.6 * ((atom_id - 1) * 2 % total_atoms) / total_atoms
                    z = 0.3 + 0.4 * ((atom_id - 1) * 3 % total_atoms) / total_atoms
                    
                    positions.append({
                        'label': f'{element}{atom_id}',
                        'element': element,
                        'x': x, 'y': y, 'z': z,
                        'occupancy': 1.0
                    })
                    atom_id += 1
        
        return positions
    
    def generate_cif_content(self, formula: str, composition: dict, 
                           category: str, lattice_params: dict, 
                           positions: list) -> str:
        """Generate CIF file content"""
        
        cif_content = f"""#######################################################################
#
# CIF file for {formula} (ambient-pressure superconductor candidate)
# Generated by RT-Superconductor Discovery Pipeline
# Category: {category}
# Expected self-compression: >200 GPa
#
#######################################################################

data_{formula}

_chemical_name_common                  '{formula}'
_chemical_formula_sum                  '{formula}'
_chemical_formula_weight               {self.calculate_formula_weight(composition):.2f}

_space_group_name_H-M_alt              '{lattice_params["space_group"]}'
_space_group_IT_number                 1

_cell_length_a                         {lattice_params["a"]:.6f}
_cell_length_b                         {lattice_params["b"]:.6f}
_cell_length_c                         {lattice_params["c"]:.6f}
_cell_angle_alpha                      {lattice_params["alpha"]:.2f}
_cell_angle_beta                       {lattice_params["beta"]:.2f}
_cell_angle_gamma                      {lattice_params["gamma"]:.2f}
_cell_volume                           {self.calculate_cell_volume(lattice_params):.2f}

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
"""
        
        for pos in positions:
            cif_content += f"{pos['label']:8} {pos['element']:2} {pos['x']:8.5f} {pos['y']:8.5f} {pos['z']:8.5f} {pos['occupancy']:6.3f}\n"
        
        return cif_content
    
    def calculate_formula_weight(self, composition: dict) -> float:
        """Calculate molecular weight from composition"""
        atomic_weights = {
            'H': 1.008, 'Li': 6.94, 'Be': 9.01, 'B': 10.81, 'C': 12.01,
            'N': 14.01, 'O': 16.00, 'Ca': 40.08, 'Y': 88.91, 'La': 138.91
        }
        
        weight = 0
        for element, count in composition.items():
            weight += count * atomic_weights.get(element, 12.0)
        
        return weight
    
    def calculate_cell_volume(self, lattice_params: dict) -> float:
        """Calculate unit cell volume"""
        a, b, c = lattice_params['a'], lattice_params['b'], lattice_params['c']
        alpha = math.radians(lattice_params['alpha'])
        beta = math.radians(lattice_params['beta'])
        gamma = math.radians(lattice_params['gamma'])
        
        volume = a * b * c * math.sqrt(
            1 + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma) -
            math.cos(alpha)**2 - math.cos(beta)**2 - math.cos(gamma)**2
        )
        
        return volume
    
    def parse_formula_to_composition(self, formula: str) -> dict:
        """Parse chemical formula string to composition dictionary"""
        import re
        
        composition = {}
        
        # Find all element-count pairs
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)
        
        for element, count_str in matches:
            count = int(count_str) if count_str else 1
            composition[element] = count
        
        return composition

def main():
    parser = argparse.ArgumentParser(description='Export superconductor candidates to CIF format')
    parser.add_argument('--input', type=str, required=True,
                       help='Input JSON file with chemical pressure analysis')
    parser.add_argument('--top-n', type=int, default=50,
                       help='Number of top candidates to export')
    parser.add_argument('--out', type=str, required=True,
                       help='Output directory for CIF files')
    
    args = parser.parse_args()
    
    print("🔬 CIF EXPORT TOOL FOR AMBIENT-PRESSURE SUPERCONDUCTORS")
    print("=" * 60)
    
    # Load chemical pressure analysis
    with open(args.input, 'r') as f:
        data = json.load(f)
    
    candidates = data['chemical_pressure_analysis'][:args.top_n]
    
    # Create output directory
    output_dir = Path(args.out)
    output_dir.mkdir(exist_ok=True)
    
    print(f"📊 Exporting {len(candidates)} candidates to {output_dir}")
    
    # Initialize CIF generator
    generator = CIFGenerator()
    
    exported_count = 0
    for i, candidate in enumerate(candidates):
        try:
            formula = candidate['formula']
            category = candidate['category']
            
            print(f"  {i+1:2d}. {formula:15} | {category:10} | P={candidate['chemical_pressure']:.0f} GPa")
            
            # Parse composition
            composition = generator.parse_formula_to_composition(formula)
            
            # Generate lattice parameters
            lattice_params = generator.estimate_lattice_parameters(composition, category)
            
            # Generate atomic positions
            positions = generator.generate_atomic_positions(composition, lattice_params)
            
            # Generate CIF content
            cif_content = generator.generate_cif_content(
                formula, composition, category, lattice_params, positions
            )
            
            # Save CIF file
            cif_filename = output_dir / f"{formula}.cif"
            with open(cif_filename, 'w') as f:
                f.write(cif_content)
            
            exported_count += 1
            
        except Exception as e:
            print(f"    ❌ Error exporting {formula}: {e}")
            continue
    
    print("")
    print("✅ CIF EXPORT COMPLETE!")
    print("=" * 30)
    print(f"📁 Files exported: {exported_count}")
    print(f"📂 Output directory: {output_dir}")
    print("")
    print("🚀 Ready for xTB stability screening!")
    print("   Next command: cd calc/00_xtb_cifs && xtb *.cif --gfn 2 --opt --freq")

if __name__ == "__main__":
    main() 