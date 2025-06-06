#!/usr/bin/env python
"""
Targeted enumeration for practical superconductor candidates.

Focuses on common, scalable elements while searching for RBT-favorable compositions.
Balances theoretical RBT principles with practical material constraints.
"""

import argparse
import json
import itertools
from pathlib import Path
from typing import List, Dict, Any, Set
import numpy as np
from pymatgen.core import Composition, Element
from tqdm import tqdm

import sys
sys.path.append('.')
from src.utils.ledger import tau_valence_mismatch, nearest_split

# Practical elements for scalable superconductors
PRACTICAL_ELEMENTS = [
    # Abundant light metals
    "Li", "Mg", "Al", "Ca", "Sr", "Ba",
    
    # Common transition metals
    "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Zr", "Nb", "Mo", "Ag", "Cd", "Sn",
    
    # p-block elements
    "B", "C", "N", "O", "F", "Si", "P", "S", "Cl",
    "Ga", "Ge", "As", "Se", "Br", "In", "Sb", "Te", "I",
    
    # Rare earths (some are accessible)
    "La", "Ce", "Pr", "Nd", "Y"
]

# Known superconductor element combinations to guide search
SC_ELEMENT_GROUPS = [
    ["Mg", "B"],           # MgB2 family
    ["La", "Fe", "As", "O"], # Pnictides  
    ["Y", "Ba", "Cu", "O"],  # Cuprates
    ["Ca", "C"],           # Intercalates
    ["Nb", "Ti"],          # A15 compounds
    ["Fe", "Se", "S"],     # Chalcogenides
    ["Li", "Ti", "O"],     # Titanates
    ["Sr", "Ti", "O"],     # Perovskites
    ["Ba", "Ti", "O"],     # Perovskites
]

def enumerate_guided_compositions(max_atoms: int = 8, max_candidates: int = 1000) -> List[Dict[str, Any]]:
    """Enumerate compositions using both RBT principles and superconductor families."""
    
    candidates = []
    seen_formulas = set()
    
    print("🎯 Generating RBT-guided practical superconductor candidates...")
    
    # Strategy 1: RBT split-octet favorable combinations
    print("Strategy 1: Split-octet favorable combinations")
    for n_species in range(2, 5):  # Binary to quaternary
        for elements in itertools.combinations(PRACTICAL_ELEMENTS, n_species):
            
            # Try different stoichiometries
            max_per_element = max_atoms // n_species + 2
            for stoich in itertools.product(range(1, max_per_element), repeat=n_species):
                if sum(stoich) > max_atoms or sum(stoich) < 2:
                    continue
                
                comp_dict = dict(zip(elements, stoich))
                composition = Composition(comp_dict)
                reduced = composition.reduced_formula
                
                if reduced in seen_formulas:
                    continue
                seen_formulas.add(reduced)
                
                # Calculate RBT metrics
                tau = tau_valence_mismatch(composition)
                
                # Add if reasonably good RBT score (relaxed threshold for practical materials)
                if tau < 2.0:  # Much more lenient than perfect split-octet
                    candidates.append({
                        'formula': composition.formula,
                        'reduced_formula': reduced,
                        'tau_valence': tau,
                        'n_elements': len(composition.elements),
                        'n_atoms': composition.num_atoms,
                        'strategy': 'split_octet',
                        'scalability': assess_scalability(composition)
                    })
    
    # Strategy 2: Guided by known superconductor families
    print("Strategy 2: Known superconductor family variations")
    for element_group in SC_ELEMENT_GROUPS:
        
        # Generate compositions within each family
        for n_select in range(2, min(5, len(element_group) + 1)):
            for elements in itertools.combinations(element_group, n_select):
                
                for stoich in itertools.product(range(1, 6), repeat=len(elements)):
                    if sum(stoich) > max_atoms:
                        continue
                    
                    comp_dict = dict(zip(elements, stoich))
                    composition = Composition(comp_dict)
                    reduced = composition.reduced_formula
                    
                    if reduced in seen_formulas:
                        continue
                    seen_formulas.add(reduced)
                    
                    tau = tau_valence_mismatch(composition)
                    
                    # Add if reasonable (prioritize known families even with higher τ)
                    if tau < 5.0:
                        candidates.append({
                            'formula': composition.formula,
                            'reduced_formula': reduced,
                            'tau_valence': tau,
                            'n_elements': len(composition.elements),
                            'n_atoms': composition.num_atoms,
                            'strategy': 'known_family',
                            'family': str(element_group),
                            'scalability': assess_scalability(composition)
                        })
    
    # Strategy 3: Hybrid compositions (mix RBT + known SC elements)
    print("Strategy 3: Hybrid RBT + superconductor compositions")
    sc_elements = set(['Cu', 'Fe', 'Nb', 'Ti', 'La', 'Y', 'Ba', 'Mg', 'B'])
    
    for sc_elem in sc_elements:
        for other_elem in PRACTICAL_ELEMENTS:
            if other_elem == sc_elem:
                continue
            
            # Binary combinations
            for stoich in itertools.product(range(1, 6), repeat=2):
                if sum(stoich) > max_atoms:
                    continue
                
                comp_dict = {sc_elem: stoich[0], other_elem: stoich[1]}
                composition = Composition(comp_dict)
                reduced = composition.reduced_formula
                
                if reduced in seen_formulas:
                    continue
                seen_formulas.add(reduced)
                
                tau = tau_valence_mismatch(composition)
                
                if tau < 3.0:
                    candidates.append({
                        'formula': composition.formula,
                        'reduced_formula': reduced,
                        'tau_valence': tau,
                        'n_elements': len(composition.elements),
                        'n_atoms': composition.num_atoms,
                        'strategy': 'hybrid',
                        'scalability': assess_scalability(composition)
                    })
    
    print(f"Generated {len(candidates)} guided candidates")
    
    # Sort by combined score (RBT + scalability)
    def combined_score(candidate):
        tau_score = np.exp(-candidate['tau_valence'] / 2.0)  # Gentle penalty for higher τ
        scalability_score = candidate['scalability']
        return tau_score * scalability_score
    
    candidates.sort(key=combined_score, reverse=True)
    
    return candidates[:max_candidates]

def assess_scalability(composition: Composition) -> float:
    """Assess how scalable/practical a composition is for mass production."""
    
    scalability = 1.0
    
    for element in composition.elements:
        symbol = element.symbol
        
        # Cost/abundance factors
        if symbol in ['H', 'Li', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K', 'Ca', 'Fe']:
            scalability *= 1.0  # Very abundant
        elif symbol in ['Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni', 'Cu', 'Zn', 'Sr', 'Zr', 'Nb', 'Mo', 'Ba', 'La']:
            scalability *= 0.8  # Common but more expensive
        elif symbol in ['Y', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te']:
            scalability *= 0.6  # Less common
        elif symbol in ['Ce', 'Pr', 'Nd']:
            scalability *= 0.4  # Rare earths
        else:
            scalability *= 0.2  # Expensive/rare elements
        
        # Toxicity factors
        if symbol in ['Hg', 'Pb', 'Cd', 'As']:
            scalability *= 0.1  # Toxic
        elif symbol in ['Cr', 'Ni', 'Co']:
            scalability *= 0.9  # Mildly concerning
    
    return scalability

def main():
    parser = argparse.ArgumentParser(description="Generate practical superconductor candidates")
    parser.add_argument("--max-atoms", type=int, default=8, help="Maximum atoms per formula")
    parser.add_argument("--max-candidates", type=int, default=500, help="Maximum candidates to generate")
    parser.add_argument("--out", required=True, help="Output JSON file")
    
    args = parser.parse_args()
    
    # Generate candidates
    candidates = enumerate_guided_compositions(args.max_atoms, args.max_candidates)
    
    # Save results
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(candidates, f, indent=2)
    
    print(f"\n✅ Saved {len(candidates)} practical candidates to {output_path}")
    
    # Show top candidates
    print("\n🏆 Top 15 practical superconductor candidates:")
    print("Rank | Formula | τ | Strategy | Scalability | Combined Score")
    print("-" * 65)
    
    for i, candidate in enumerate(candidates[:15]):
        tau_score = np.exp(-candidate['tau_valence'] / 2.0)
        combined = tau_score * candidate['scalability']
        
        print(f"{i+1:4d} | {candidate['reduced_formula']:12s} | "
              f"{candidate['tau_valence']:5.3f} | "
              f"{candidate['strategy']:12s} | "
              f"{candidate['scalability']:10.3f} | "
              f"{combined:12.4f}")

if __name__ == "__main__":
    main() 