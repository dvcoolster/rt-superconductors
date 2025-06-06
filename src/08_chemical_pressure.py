#!/usr/bin/env python3
"""
CHEMICAL PRESSURE CALCULATOR FOR AMBIENT-PRESSURE SUPERCONDUCTORS

Computes ΔV_cell per electron to identify materials with intrinsic 
chemical pre-compression that can achieve superconductivity at ambient pressure.

Key concepts:
- Chemical pressure: Internal stress from size/charge mismatches
- Self-compression: Materials that compress themselves via chemistry
- Volume efficiency: How much electron density fits in minimal space

Author: RT-Superconductor Research Team  
Date: 2025-06-07
"""

import numpy as np
import json
import argparse
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import math

@dataclass
class ChemicalPressureData:
    """Store chemical pressure calculation results"""
    formula: str
    v_cell: float  # Unit cell volume (Å³)
    n_electrons: int  # Total valence electrons
    delta_v_per_e: float  # ΔV_cell per electron (Å³/e⁻)
    chemical_pressure: float  # Estimated pressure (GPa)
    compression_ratio: float  # Relative to reference state
    self_compressed: bool  # Flagged as self-compressed
    category: str  # Material type

class ChemicalPressureCalculator:
    """Calculate chemical pressure and electron density efficiency"""
    
    def __init__(self):
        # Element properties
        self.element_data = {
            'H': {'radius': 0.37, 'valence': 1, 'ref_volume': 2.0},
            'Li': {'radius': 1.52, 'valence': 1, 'ref_volume': 13.0},
            'Be': {'radius': 1.12, 'valence': 2, 'ref_volume': 5.0},
            'B': {'radius': 0.82, 'valence': 3, 'ref_volume': 4.6},
            'C': {'radius': 0.77, 'valence': 4, 'ref_volume': 3.4},
            'N': {'radius': 0.75, 'valence': 3, 'ref_volume': 2.8},
            'O': {'radius': 0.73, 'valence': 2, 'ref_volume': 2.2},
            'Ca': {'radius': 1.97, 'valence': 2, 'ref_volume': 29.9},
            'Y': {'radius': 1.80, 'valence': 3, 'ref_volume': 19.9},
            'La': {'radius': 1.87, 'valence': 3, 'ref_volume': 22.5}
        }
    
    def estimate_unit_cell_volume(self, composition: Dict[str, int]) -> float:
        """Estimate unit cell volume from composition"""
        # Method 1: Additive atomic volumes with compression factors
        total_volume = 0
        total_atoms = sum(composition.values())
        
        for element, count in composition.items():
            atomic_radius = self.element_data[element]['radius']
            # Volume per atom (assuming close packing)
            atomic_volume = (4/3) * math.pi * atomic_radius**3
            total_volume += count * atomic_volume
        
        # Apply compression factors based on composition
        compression_factor = self.calculate_compression_factor(composition)
        compressed_volume = total_volume * compression_factor
        
        return compressed_volume
    
    def calculate_compression_factor(self, composition: Dict[str, int]) -> float:
        """Calculate compression factor from chemical effects"""
        total_atoms = sum(composition.values())
        
        # Size mismatch factor
        radii = [self.element_data[el]['radius'] for el in composition.keys()]
        if len(radii) > 1:
            radius_range = max(radii) - min(radii)
            size_mismatch = radius_range / max(radii)
        else:
            size_mismatch = 0
        
        # Charge mismatch factor  
        valences = [self.element_data[el]['valence'] for el in composition.keys()]
        if len(valences) > 1:
            valence_range = max(valences) - min(valences)
            charge_mismatch = valence_range / max(valences)
        else:
            charge_mismatch = 0
        
        # H content factor (H creates large compression)
        h_content = composition.get('H', 0) / total_atoms
        h_compression = h_content * 0.3  # Up to 30% compression from H
        
        # Mixed cation factor
        metals = sum(1 for el in ['Li', 'Be', 'Ca', 'Y', 'La'] if el in composition)
        mixed_cation_compression = min(metals - 1, 1) * 0.2  # Up to 20% for mixed cations
        
        # Total compression (multiplicative effects)
        total_compression = (size_mismatch * 0.4 + charge_mismatch * 0.3 + 
                           h_compression + mixed_cation_compression)
        
        # Compression factor: 1.0 = no compression, 0.5 = 50% compressed
        compression_factor = 1.0 - min(total_compression, 0.8)  # Max 80% compression
        
        return compression_factor
    
    def calculate_reference_volume(self, composition: Dict[str, int]) -> float:
        """Calculate volume of reference (uncompressed) state"""
        ref_volume = 0
        for element, count in composition.items():
            element_ref_vol = self.element_data[element]['ref_volume']
            ref_volume += count * element_ref_vol
        
        return ref_volume
    
    def calculate_total_electrons(self, composition: Dict[str, int]) -> int:
        """Calculate total valence electrons"""
        total_electrons = 0
        for element, count in composition.items():
            valence = self.element_data[element]['valence']
            total_electrons += count * valence
        
        return total_electrons
    
    def estimate_chemical_pressure(self, delta_v_per_e: float, 
                                 compression_ratio: float) -> float:
        """Estimate chemical pressure from volume compression"""
        # Empirical relationship: P ∝ -ln(V/V₀) / compressibility
        if compression_ratio >= 1.0:
            return 0.0  # No compression
        
        # Use bulk modulus relationship
        # P = -B * ln(V/V₀) where B ~ 100 GPa for typical materials
        bulk_modulus = 100.0  # GPa (typical value)
        pressure = -bulk_modulus * math.log(compression_ratio)
        
        # Additional pressure from electron density
        # Higher electron density → higher pressure
        if delta_v_per_e < 5.0:  # Very dense electron packing
            pressure += 50.0
        elif delta_v_per_e < 10.0:  # Moderate density
            pressure += 20.0
        
        return pressure
    
    def classify_material_type(self, composition: Dict[str, int]) -> str:
        """Classify material by superconductor category"""
        elements = set(composition.keys())
        h_fraction = composition.get('H', 0) / sum(composition.values())
        
        if 'Ca' in elements and h_fraction > 0.3:
            return "Electride"
        elif len(elements & {'Y', 'La'}) > 0 and 'Be' in elements:
            return "Clathrate"
        elif elements & {'Li', 'B', 'C', 'N'} and h_fraction < 0.4:
            return "Li-B-C-N"
        elif h_fraction > 0.2:
            return "Hydride"
        else:
            return "Other"
    
    def calculate_chemical_pressure(self, composition: Dict[str, int], 
                                  formula: str) -> ChemicalPressureData:
        """Calculate complete chemical pressure analysis"""
        
        # Basic quantities
        v_cell = self.estimate_unit_cell_volume(composition)
        v_ref = self.calculate_reference_volume(composition)
        n_electrons = self.calculate_total_electrons(composition)
        
        # Key metrics
        delta_v_per_e = v_cell / n_electrons if n_electrons > 0 else 100.0
        compression_ratio = v_cell / v_ref if v_ref > 0 else 1.0
        chemical_pressure = self.estimate_chemical_pressure(delta_v_per_e, compression_ratio)
        
        # Classification
        category = self.classify_material_type(composition)
        
        # Flag as self-compressed if meets criteria
        self_compressed = (delta_v_per_e < 5.0 and compression_ratio < 0.8 and 
                          chemical_pressure > 10.0)
        
        return ChemicalPressureData(
            formula=formula,
            v_cell=v_cell,
            n_electrons=n_electrons,
            delta_v_per_e=delta_v_per_e,
            chemical_pressure=chemical_pressure,
            compression_ratio=compression_ratio,
            self_compressed=self_compressed,
            category=category
        )
    
    def analyze_candidates(self, candidates_file: str) -> List[ChemicalPressureData]:
        """Analyze chemical pressure for enumerated candidates"""
        print(f"🔧 Analyzing chemical pressure from {candidates_file}")
        
        with open(candidates_file, 'r') as f:
            data = json.load(f)
        
        results = []
        candidates = data.get('candidates', [])
        
        for candidate in candidates:
            formula = candidate['formula']
            composition = candidate['elements']
            
            # Convert string keys to proper composition dict
            comp_dict = {}
            for element, count in composition.items():
                comp_dict[element] = int(count)
            
            # Calculate chemical pressure
            cp_data = self.calculate_chemical_pressure(comp_dict, formula)
            results.append(cp_data)
        
        print(f"✅ Analyzed {len(results)} candidates")
        return results
    
    def rank_by_chemical_pressure(self, results: List[ChemicalPressureData]) -> List[ChemicalPressureData]:
        """Rank candidates by chemical pressure potential"""
        print("📊 Ranking by chemical pressure...")
        
        # Score based on multiple factors
        for result in results:
            score = 0
            
            # High chemical pressure is good
            if result.chemical_pressure > 50:
                score += 3
            elif result.chemical_pressure > 20:
                score += 2
            elif result.chemical_pressure > 10:
                score += 1
            
            # Low ΔV per electron is good (dense packing)
            if result.delta_v_per_e < 3:
                score += 3
            elif result.delta_v_per_e < 5:
                score += 2
            elif result.delta_v_per_e < 8:
                score += 1
            
            # High compression ratio is good
            if result.compression_ratio < 0.5:
                score += 3
            elif result.compression_ratio < 0.7:
                score += 2
            elif result.compression_ratio < 0.9:
                score += 1
            
            # Self-compressed materials get bonus
            if result.self_compressed:
                score += 2
            
            # Category bonuses
            category_bonus = {
                'Electride': 2,
                'Clathrate': 2,
                'Li-B-C-N': 1,
                'Hydride': 1,
                'Other': 0
            }
            score += category_bonus.get(result.category, 0)
            
            result.pressure_score = score
        
        # Sort by score (descending)
        ranked = sorted(results, key=lambda x: x.pressure_score, reverse=True)
        
        print(f"🏆 Top candidate: {ranked[0].formula} (pressure score: {ranked[0].pressure_score})")
        return ranked

def main():
    parser = argparse.ArgumentParser(description='Calculate chemical pressure for superconductor candidates')
    parser.add_argument('--input', type=str, required=True,
                       help='Input JSON file with candidates')
    parser.add_argument('--output', type=str, required=True,
                       help='Output JSON file with chemical pressure data')
    parser.add_argument('--top-n', type=int, default=50,
                       help='Number of top candidates to analyze')
    
    args = parser.parse_args()
    
    print("🌡️ CHEMICAL PRESSURE ANALYSIS FOR AMBIENT SUPERCONDUCTORS")
    print("=" * 60)
    
    # Initialize calculator
    calculator = ChemicalPressureCalculator()
    
    # Analyze candidates
    results = calculator.analyze_candidates(args.input)
    
    # Rank by chemical pressure
    ranked_results = calculator.rank_by_chemical_pressure(results)
    
    # Take top N
    top_results = ranked_results[:args.top_n]
    
    # Prepare output data
    output_data = {
        'metadata': {
            'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'input_file': args.input,
            'total_analyzed': len(results),
            'top_n_saved': len(top_results)
        },
        'chemical_pressure_analysis': []
    }
    
    for i, result in enumerate(top_results):
        analysis_data = {
            'rank': i + 1,
            'formula': result.formula,
            'v_cell': round(result.v_cell, 2),
            'n_electrons': result.n_electrons,
            'delta_v_per_e': round(result.delta_v_per_e, 3),
            'chemical_pressure': round(result.chemical_pressure, 1),
            'compression_ratio': round(result.compression_ratio, 3),
            'self_compressed': result.self_compressed,
            'category': result.category,
            'pressure_score': result.pressure_score
        }
        output_data['chemical_pressure_analysis'].append(analysis_data)
    
    # Save results
    with open(args.output, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print("")
    print("🎯 CHEMICAL PRESSURE ANALYSIS COMPLETE!")
    print("=" * 40)
    print(f"📊 Candidates analyzed: {len(results)}")
    print(f"🏆 Top candidates saved: {len(top_results)}")
    print(f"💾 Output saved to: {args.output}")
    print("")
    
    # Show top 5 self-compressed candidates
    self_compressed = [r for r in top_results if r.self_compressed]
    print(f"🔥 TOP SELF-COMPRESSED CANDIDATES ({len(self_compressed)} found):")
    for i, result in enumerate(self_compressed[:5]):
        print(f"  {i+1}. {result.formula:15} | ΔV/e⁻={result.delta_v_per_e:.2f} | "
              f"P={result.chemical_pressure:.0f} GPa | {result.category}")
    
    if len(self_compressed) == 0:
        print("  No self-compressed candidates found with current criteria.")
        print("  Consider relaxing thresholds or expanding search space.")
    
    print("")
    print("🚀 Ready for xTB geometry optimization!")

if __name__ == "__main__":
    import time
    main() 