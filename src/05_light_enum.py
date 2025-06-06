#!/usr/bin/env python3
"""
AMBIENT-PRESSURE SUPERCONDUCTOR ENUMERATION PIPELINE

Generates candidate structures for room-temperature superconductors
that work at ambient pressure by focusing on:
- Hydrogen-rich covalent electrides 
- Metal-hydride clathrate alloys with chemical pre-compression
- Layered Li-B-C-N networks
- Interface-engineered materials

Author: RT-Superconductor Research Team
Date: 2025-06-07
"""

import itertools
import json
import argparse
import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Any, Tuple
import sys
import time

@dataclass
class Composition:
    """Store composition and computed properties"""
    formula: str
    elements: Dict[str, int]
    n_atoms: int
    tau: float  # Split-octet balance parameter
    kappa: float  # Electronic compressibility
    phase_lock: float  # Phase-lock metric
    v_per_h: float  # Volume per H atom (Å³)
    bulk_modulus: float  # GPa
    formation_energy: float  # eV/atom
    category: str  # Material category

class AmbientSuperconductorEnumerator:
    """Generate and score ambient-pressure superconductor candidates"""
    
    def __init__(self, elements: List[str]):
        self.elements = elements
        self.compositions = []
        
        # Element properties for scoring
        self.element_data = {
            'H': {'weight': 1.008, 'radius': 0.37, 'electronegativity': 2.20, 'valence': 1},
            'Li': {'weight': 6.94, 'radius': 1.52, 'electronegativity': 0.98, 'valence': 1},
            'Be': {'weight': 9.01, 'radius': 1.12, 'electronegativity': 1.57, 'valence': 2},
            'B': {'weight': 10.81, 'radius': 0.82, 'electronegativity': 2.04, 'valence': 3},
            'C': {'weight': 12.01, 'radius': 0.77, 'electronegativity': 2.55, 'valence': 4},
            'N': {'weight': 14.01, 'radius': 0.75, 'electronegativity': 3.04, 'valence': 3},
            'O': {'weight': 16.00, 'radius': 0.73, 'electronegativity': 3.44, 'valence': 2},
            'Ca': {'weight': 40.08, 'radius': 1.97, 'electronegativity': 1.00, 'valence': 2},
            'Y': {'weight': 88.91, 'radius': 1.80, 'electronegativity': 1.22, 'valence': 3},
            'La': {'weight': 138.91, 'radius': 1.87, 'electronegativity': 1.10, 'valence': 3}
        }
    
    def calculate_tau(self, composition: Dict[str, int]) -> float:
        """Calculate split-octet balance parameter τ"""
        total_valence = sum(self.element_data[el]['valence'] * count 
                          for el, count in composition.items())
        total_atoms = sum(composition.values())
        
        # Ideal electron count for superconductivity (empirical)
        ideal_electrons_per_atom = 2.5
        actual_electrons_per_atom = total_valence / total_atoms
        
        # τ = 0 is perfect balance, τ > 1 is poor
        tau = abs(actual_electrons_per_atom - ideal_electrons_per_atom) / ideal_electrons_per_atom
        return min(tau, 2.0)  # Cap at 2.0
    
    def calculate_kappa(self, composition: Dict[str, int]) -> float:
        """Calculate electronic compressibility κ"""
        # Based on electronegativity difference and bond covalency
        elements = list(composition.keys())
        if len(elements) < 2:
            return 1.0  # Single element = high compressibility
        
        electronegativities = [self.element_data[el]['electronegativity'] for el in elements]
        en_range = max(electronegativities) - min(electronegativities)
        
        # Lower electronegativity difference → lower κ (better for SC)
        kappa = en_range / 3.0  # Normalize by max possible difference
        return min(kappa, 1.0)
    
    def calculate_phase_lock(self, composition: Dict[str, int]) -> float:
        """Calculate phase-lock metric for superconducting stability"""
        # Based on structural rigidity and electronic homogeneity
        h_content = composition.get('H', 0) / sum(composition.values())
        light_content = sum(composition.get(el, 0) for el in ['H', 'Li', 'B', 'C', 'N']) / sum(composition.values())
        
        # High H content + light elements = good phase lock
        phase_lock = (h_content * 0.7 + light_content * 0.3)
        return min(phase_lock, 1.0)
    
    def estimate_bulk_modulus(self, composition: Dict[str, int]) -> float:
        """Estimate bulk modulus from composition"""
        # Weighted average based on element properties
        total_weight = 0
        weighted_modulus = 0
        
        # Approximate bulk moduli (GPa)
        moduli = {'H': 1.0, 'Li': 11, 'Be': 130, 'B': 320, 'C': 440, 
                 'N': 2.0, 'O': 1.0, 'Ca': 17, 'Y': 41, 'La': 28}
        
        for el, count in composition.items():
            weight = count * self.element_data[el]['weight']
            total_weight += weight
            weighted_modulus += weight * moduli.get(el, 50)
        
        return weighted_modulus / total_weight if total_weight > 0 else 50
    
    def estimate_volume_per_h(self, composition: Dict[str, int]) -> float:
        """Estimate volume per H atom (chemical pressure indicator)"""
        h_count = composition.get('H', 0)
        if h_count == 0:
            return 100.0  # No H = not relevant
        
        # Estimate unit cell volume from atomic radii
        total_volume = 0
        for el, count in composition.items():
            radius = self.element_data[el]['radius']
            atomic_volume = (4/3) * np.pi * radius**3
            total_volume += count * atomic_volume
        
        # Convert to Å³ per H atom
        volume_per_h = total_volume / h_count
        return volume_per_h
    
    def estimate_formation_energy(self, composition: Dict[str, int]) -> float:
        """Rough formation energy estimate (eV/atom)"""
        # Simple bond energy model
        total_atoms = sum(composition.values())
        
        # Favorable combinations get negative formation energy
        h_content = composition.get('H', 0) / total_atoms
        metal_content = sum(composition.get(el, 0) for el in ['Ca', 'Y', 'La']) / total_atoms
        
        # Metal hydrides are generally stable
        base_energy = -0.5 * h_content * metal_content
        
        # Add small random variation
        variation = np.random.normal(0, 0.1)
        return base_energy + variation
    
    def categorize_material(self, composition: Dict[str, int]) -> str:
        """Classify material into target categories"""
        elements = set(composition.keys())
        h_content = composition.get('H', 0) / sum(composition.values())
        
        if 'Ca' in elements and h_content > 0.3:
            return "Electride"
        elif len(elements & {'Y', 'La'}) > 0 and 'Be' in elements and h_content > 0.4:
            return "Clathrate"
        elif elements & {'Li', 'B', 'C', 'N'} and h_content < 0.3:
            return "Li-B-C-N"
        elif h_content > 0.2:
            return "Hydride"
        else:
            return "Other"
    
    def generate_compositions(self, max_atoms: int, mixed_cation: bool = True) -> List[Composition]:
        """Generate all possible compositions up to max_atoms"""
        print(f"🔧 Generating compositions with ≤{max_atoms} atoms...")
        
        compositions = []
        count = 0
        
        # Define composition ranges for each element
        max_counts = {
            'H': min(max_atoms // 2, 10),  # Up to 10 H atoms
            'Li': min(max_atoms // 3, 4),   # Up to 4 Li atoms
            'Be': min(max_atoms // 3, 3),   # Up to 3 Be atoms
            'B': min(max_atoms // 3, 3),    # Up to 3 B atoms  
            'C': min(max_atoms // 3, 4),    # Up to 4 C atoms
            'N': min(max_atoms // 3, 4),    # Up to 4 N atoms
            'O': min(max_atoms // 4, 3),    # Up to 3 O atoms
            'Ca': min(max_atoms // 4, 2),   # Up to 2 Ca atoms
            'Y': min(max_atoms // 5, 2),    # Up to 2 Y atoms
            'La': min(max_atoms // 6, 2)    # Up to 2 La atoms
        }
        
        # Generate compositions systematically
        for n_elements in range(2, min(5, len(self.elements) + 1)):  # 2-4 different elements
            for element_combo in itertools.combinations(self.elements, n_elements):
                # Must include H for most cases
                if 'H' not in element_combo and np.random.random() > 0.1:
                    continue
                
                # Generate count combinations
                count_ranges = [range(1, max_counts.get(el, 1) + 1) for el in element_combo]
                
                for counts in itertools.product(*count_ranges):
                    if sum(counts) > max_atoms or sum(counts) < 2:
                        continue
                    
                    # Create composition dictionary
                    comp_dict = dict(zip(element_combo, counts))
                    
                    # Skip if no mixed cations when required
                    if mixed_cation:
                        metals = sum(1 for el in ['Li', 'Be', 'Ca', 'Y', 'La'] if el in comp_dict)
                        if metals < 1:
                            continue
                    
                    # Calculate properties
                    formula = "".join(f"{el}{count}" if count > 1 else el 
                                    for el, count in sorted(comp_dict.items()))
                    
                    tau = self.calculate_tau(comp_dict)
                    kappa = self.calculate_kappa(comp_dict)  
                    phase_lock = self.calculate_phase_lock(comp_dict)
                    v_per_h = self.estimate_volume_per_h(comp_dict)
                    bulk_modulus = self.estimate_bulk_modulus(comp_dict)
                    formation_energy = self.estimate_formation_energy(comp_dict)
                    category = self.categorize_material(comp_dict)
                    
                    composition = Composition(
                        formula=formula,
                        elements=comp_dict,
                        n_atoms=sum(counts),
                        tau=tau,
                        kappa=kappa,
                        phase_lock=phase_lock,
                        v_per_h=v_per_h,
                        bulk_modulus=bulk_modulus,
                        formation_energy=formation_energy,
                        category=category
                    )
                    
                    compositions.append(composition)
                    count += 1
                    
                    if count % 1000 == 0:
                        print(f"  Generated {count} compositions...")
        
        print(f"✅ Generated {len(compositions)} total compositions")
        return compositions
    
    def filter_compositions(self, compositions: List[Composition], 
                          filter_str: str) -> List[Composition]:
        """Filter compositions based on criteria string"""
        print(f"🔍 Filtering with criteria: {filter_str}")
        
        filtered = []
        for comp in compositions:
            # Parse filter string (simplified)
            criteria = filter_str.split(' & ')
            passes = True
            
            for criterion in criteria:
                criterion = criterion.strip()
                if 'tau<' in criterion:
                    threshold = float(criterion.split('<')[1])
                    if comp.tau >= threshold:
                        passes = False
                        break
                elif 'kappa<' in criterion:
                    threshold = float(criterion.split('<')[1])
                    if comp.kappa >= threshold:
                        passes = False
                        break
                elif 'phase_lock>' in criterion:
                    threshold = float(criterion.split('>')[1])
                    if comp.phase_lock <= threshold:
                        passes = False
                        break
            
            if passes:
                filtered.append(comp)
        
        print(f"✅ {len(filtered)} compositions passed filters")
        return filtered
    
    def score_and_rank(self, compositions: List[Composition]) -> List[Composition]:
        """Score compositions for superconducting potential"""
        print("📊 Scoring and ranking compositions...")
        
        for comp in compositions:
            # Multi-objective scoring function
            tau_score = max(0, 1 - comp.tau)  # Lower τ is better
            kappa_score = max(0, 1 - comp.kappa)  # Lower κ is better  
            phase_lock_score = comp.phase_lock  # Higher phase_lock is better
            
            # Chemical pressure bonus (lower V/H = higher internal pressure)
            if comp.v_per_h < 5.0:  # Very compressed
                pressure_score = 1.0
            elif comp.v_per_h < 10.0:  # Moderately compressed
                pressure_score = 0.5
            else:
                pressure_score = 0.0
            
            # Formation energy penalty (too unstable = bad)
            stability_score = max(0, 1 + comp.formation_energy)  # More negative = better
            
            # Category bonuses
            category_bonus = {
                'Electride': 0.3,
                'Clathrate': 0.25, 
                'Li-B-C-N': 0.2,
                'Hydride': 0.1,
                'Other': 0.0
            }.get(comp.category, 0.0)
            
            # Overall score (0-1 scale)
            total_score = (tau_score * 0.25 + kappa_score * 0.25 + 
                          phase_lock_score * 0.20 + pressure_score * 0.15 +
                          stability_score * 0.10 + category_bonus * 0.05)
            
            comp.score = total_score
        
        # Sort by score (descending)
        ranked = sorted(compositions, key=lambda x: x.score, reverse=True)
        
        print(f"🏆 Top candidate: {ranked[0].formula} (score: {ranked[0].score:.3f})")
        return ranked

def main():
    parser = argparse.ArgumentParser(description='Enumerate ambient-pressure superconductor candidates')
    parser.add_argument('--elements', nargs='+', 
                       default=['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'Ca', 'Y', 'La'],
                       help='Elements to include in enumeration')
    parser.add_argument('--mixed-cation', type=bool, default=True,
                       help='Require mixed cation compositions')
    parser.add_argument('--max-atoms', type=int, default=12,
                       help='Maximum atoms per unit cell')
    parser.add_argument('--filter', type=str, 
                       default="tau<0.8 & kappa<0.2 & phase_lock>0.7",
                       help='Filter criteria string')
    parser.add_argument('--out', type=str, required=True,
                       help='Output JSON file')
    parser.add_argument('--top-n', type=int, default=100,
                       help='Number of top candidates to save')
    
    args = parser.parse_args()
    
    print("🌡️ AMBIENT-PRESSURE SUPERCONDUCTOR ENUMERATION")
    print("=" * 50)
    print(f"Elements: {args.elements}")
    print(f"Max atoms: {args.max_atoms}")
    print(f"Mixed cation: {args.mixed_cation}")
    print(f"Filter: {args.filter}")
    print("")
    
    start_time = time.time()
    
    # Initialize enumerator
    enumerator = AmbientSuperconductorEnumerator(args.elements)
    
    # Generate compositions
    compositions = enumerator.generate_compositions(args.max_atoms, args.mixed_cation)
    
    # Apply filters
    if args.filter:
        compositions = enumerator.filter_compositions(compositions, args.filter)
    
    # Score and rank
    ranked_compositions = enumerator.score_and_rank(compositions)
    
    # Take top N
    top_candidates = ranked_compositions[:args.top_n]
    
    # Convert to JSON-serializable format
    output_data = {
        'metadata': {
            'generated_at': time.strftime('%Y-%m-%d %H:%M:%S'),
            'elements': args.elements,
            'max_atoms': args.max_atoms,
            'mixed_cation': args.mixed_cation,
            'filter_criteria': args.filter,
            'total_generated': len(compositions),
            'top_n_saved': len(top_candidates)
        },
        'candidates': []
    }
    
    for i, comp in enumerate(top_candidates):
        candidate_data = {
            'rank': i + 1,
            'formula': comp.formula,
            'elements': comp.elements,
            'n_atoms': comp.n_atoms,
            'tau': round(comp.tau, 3),
            'kappa': round(comp.kappa, 3),
            'phase_lock': round(comp.phase_lock, 3),
            'v_per_h': round(comp.v_per_h, 2),
            'bulk_modulus': round(comp.bulk_modulus, 1),
            'formation_energy': round(comp.formation_energy, 3),
            'category': comp.category,
            'score': round(comp.score, 4)
        }
        output_data['candidates'].append(candidate_data)
    
    # Save results
    with open(args.out, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    elapsed = time.time() - start_time
    
    print("")
    print("🎯 ENUMERATION COMPLETE!")
    print("=" * 30)
    print(f"⏱️  Time elapsed: {elapsed:.1f} seconds")
    print(f"📊 Total candidates: {len(compositions)}")
    print(f"🏆 Top candidates saved: {len(top_candidates)}")
    print(f"💾 Output saved to: {args.out}")
    print("")
    
    # Show top 5 candidates
    print("🔥 TOP 5 CANDIDATES:")
    for i, comp in enumerate(top_candidates[:5]):
        print(f"  {i+1}. {comp.formula:15} | τ={comp.tau:.2f} κ={comp.kappa:.2f} | "
              f"Score={comp.score:.3f} | {comp.category}")
    
    print("")
    print("🚀 Ready for xTB filtering and QE calculations!")

if __name__ == "__main__":
    main() 