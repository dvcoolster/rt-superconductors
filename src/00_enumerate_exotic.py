#!/usr/bin/env python
"""
Exotic superconductor composition enumeration.

Explores more unusual element combinations, rare earth elements, 
and complex multi-component alloys using RBT principles.
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

# Exotic element combinations for advanced superconductors
EXOTIC_ELEMENT_GROUPS = {
    # Rare Earth + Transition Metal combinations
    'rare_earth_tm': [
        ["La", "Ce", "Pr", "Nd", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y", "Sc"],
        ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag"]
    ],
    
    # Heavy element combinations  
    'heavy_metals': [
        ["Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi"],
        ["In", "Sn", "Sb", "Te", "I", "Cs", "Ba"]
    ],
    
    # Chalcogenide systems
    'chalcogenides': [
        ["Fe", "Co", "Ni", "Cu", "Ag", "Cd", "In", "Sn"],
        ["S", "Se", "Te"]
    ],
    
    # Pnictide systems
    'pnictides': [
        ["Ca", "Sr", "Ba", "Eu", "La", "Ce", "Pr", "Nd"],
        ["Fe", "Co", "Ni", "Cu", "Ru", "Rh", "Pd"],
        ["P", "As", "Sb", "Bi"]
    ],
    
    # Boride systems
    'borides': [
        ["Mg", "Al", "Ca", "Sr", "Ba", "Y", "La", "Ce", "Pr", "Nd"],
        ["B"]
    ],
    
    # Carbide systems
    'carbides': [
        ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Zr", "Nb", "Mo", "Hf", "Ta", "W"],
        ["C"]
    ],
    
    # Nitride systems
    'nitrides': [
        ["Ti", "V", "Cr", "Mn", "Zr", "Nb", "Mo", "Hf", "Ta", "W"],
        ["N"]
    ],
    
    # High-entropy alloys (5+ components)
    'high_entropy': [
        ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zr", "Nb", "Mo", "Hf", "Ta", "W", "Re"]
    ],
    
    # Quantum materials (exotic properties)
    'quantum_materials': [
        ["Li", "Na", "K", "Rb", "Cs"],  # Alkalis
        ["Be", "Mg", "Ca", "Sr", "Ba"],  # Alkaline earths  
        ["Sc", "Y", "La", "Ac"],  # Group 3
        ["C", "Si", "Ge", "Sn", "Pb"]  # Carbon group
    ],
    
    # Intermetallic systems
    'intermetallics': [
        ["Al", "Ga", "In", "Tl"],  # p-block metals
        ["Ti", "Zr", "Hf"],  # Group 4
        ["V", "Nb", "Ta"],   # Group 5
        ["Cr", "Mo", "W"]    # Group 6
    ]
}

# Scalability assessment for exotic elements
ELEMENT_RARITY = {
    # Common (relative abundance > 100 ppm in crust)
    'abundant': ['Al', 'Fe', 'Ca', 'Mg', 'Ti', 'Mn', 'P', 'Ba', 'Sr', 'S', 'C', 'Cr', 'V', 'Ni', 'Zn', 'Cu', 'Co'],
    
    # Uncommon (1-100 ppm)
    'uncommon': ['La', 'Ce', 'Nd', 'Y', 'Nb', 'Mo', 'Zr', 'Pb', 'Ga', 'As', 'Ge', 'Se', 'Br', 'Rb', 'Li', 'B'],
    
    # Rare (0.01-1 ppm)
    'rare': ['Pr', 'Gd', 'Dy', 'Er', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Cs'],
    
    # Very rare (< 0.01 ppm)
    'very_rare': ['Sc', 'Eu', 'Tb', 'Ho', 'Tm', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Bi', 'Tc'],
    
    # Special cases
    'radioactive': ['Tc', 'Pm', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U']
}

def assess_exotic_scalability(composition: Composition) -> Dict[str, float]:
    """Assess scalability of exotic compositions."""
    
    scalability = 1.0
    cost_factor = 1.0
    rarity_score = 0
    
    for element in composition.elements:
        symbol = element.symbol
        
        # Rarity assessment
        if symbol in ELEMENT_RARITY['abundant']:
            rarity_score += 4
            scalability *= 1.0
            cost_factor *= 1.0
        elif symbol in ELEMENT_RARITY['uncommon']:
            rarity_score += 3
            scalability *= 0.7
            cost_factor *= 2.0
        elif symbol in ELEMENT_RARITY['rare']:
            rarity_score += 2
            scalability *= 0.4
            cost_factor *= 10.0
        elif symbol in ELEMENT_RARITY['very_rare']:
            rarity_score += 1
            scalability *= 0.1
            cost_factor *= 100.0
        elif symbol in ELEMENT_RARITY['radioactive']:
            rarity_score += 0
            scalability *= 0.01
            cost_factor *= 1000.0
    
    # Normalize rarity score
    avg_rarity = rarity_score / len(composition.elements)
    
    return {
        'scalability': scalability,
        'cost_factor': cost_factor,
        'avg_rarity': avg_rarity,
        'total_elements': len(composition.elements)
    }

def enumerate_exotic_compositions(max_atoms: int = 10, max_candidates: int = 2000) -> List[Dict[str, Any]]:
    """Enumerate exotic compositions with advanced RBT scoring."""
    
    candidates = []
    seen_formulas = set()
    
    print("🔮 Generating exotic superconductor candidates...")
    
    # Strategy 1: Systematic rare earth combinations
    print("Strategy 1: Rare Earth + Transition Metal systems")
    rare_earths = EXOTIC_ELEMENT_GROUPS['rare_earth_tm'][0][:8]  # Limit for speed
    trans_metals = EXOTIC_ELEMENT_GROUPS['rare_earth_tm'][1][:8]
    
    for re_elem in tqdm(rare_earths, desc="Rare earths"):
        for tm_elem in trans_metals:
            for stoich in itertools.product(range(1, 4), repeat=2):
                if sum(stoich) > max_atoms:
                    continue
                
                comp_dict = {re_elem: stoich[0], tm_elem: stoich[1]}
                composition = Composition(comp_dict)
                reduced = composition.reduced_formula
                
                if reduced in seen_formulas:
                    continue
                seen_formulas.add(reduced)
                
                tau = tau_valence_mismatch(composition)
                exotic_metrics = assess_exotic_scalability(composition)
                
                if tau < 4.0:  # More lenient for exotic materials
                    candidates.append({
                        'formula': composition.formula,
                        'reduced_formula': reduced,
                        'tau_valence': tau,
                        'n_elements': len(composition.elements),
                        'n_atoms': composition.num_atoms,
                        'strategy': 'rare_earth_tm',
                        'family': f"{re_elem}-{tm_elem}",
                        **exotic_metrics
                    })
    
    # Strategy 2: Chalcogenide superconductors
    print("Strategy 2: Advanced chalcogenide systems")
    metals = EXOTIC_ELEMENT_GROUPS['chalcogenides'][0]
    chalcogens = EXOTIC_ELEMENT_GROUPS['chalcogenides'][1]
    
    for metal in metals:
        for chalcogen in chalcogens:
            # Try different stoichiometries
            for m_count in range(1, 4):
                for c_count in range(1, 4):
                    if m_count + c_count > max_atoms:
                        continue
                    
                    comp_dict = {metal: m_count, chalcogen: c_count}
                    composition = Composition(comp_dict)
                    reduced = composition.reduced_formula
                    
                    if reduced in seen_formulas:
                        continue
                    seen_formulas.add(reduced)
                    
                    tau = tau_valence_mismatch(composition)
                    exotic_metrics = assess_exotic_scalability(composition)
                    
                    if tau < 6.0:
                        candidates.append({
                            'formula': composition.formula,
                            'reduced_formula': reduced,
                            'tau_valence': tau,
                            'n_elements': len(composition.elements),
                            'n_atoms': composition.num_atoms,
                            'strategy': 'chalcogenide',
                            'family': f"{metal}-{chalcogen}",
                            **exotic_metrics
                        })
    
    # Strategy 3: Complex pnictide systems
    print("Strategy 3: Multi-component pnictide systems")
    alkaline_earths = EXOTIC_ELEMENT_GROUPS['pnictides'][0][:6]
    transition_metals = EXOTIC_ELEMENT_GROUPS['pnictides'][1][:6]
    pnictogens = EXOTIC_ELEMENT_GROUPS['pnictides'][2]
    
    for ae in alkaline_earths:
        for tm in transition_metals:
            for pn in pnictogens:
                # Try ternary combinations
                for stoich in itertools.product(range(1, 3), repeat=3):
                    if sum(stoich) > max_atoms:
                        continue
                    
                    comp_dict = {ae: stoich[0], tm: stoich[1], pn: stoich[2]}
                    composition = Composition(comp_dict)
                    reduced = composition.reduced_formula
                    
                    if reduced in seen_formulas:
                        continue
                    seen_formulas.add(reduced)
                    
                    tau = tau_valence_mismatch(composition)
                    exotic_metrics = assess_exotic_scalability(composition)
                    
                    if tau < 8.0:
                        candidates.append({
                            'formula': composition.formula,
                            'reduced_formula': reduced,
                            'tau_valence': tau,
                            'n_elements': len(composition.elements),
                            'n_atoms': composition.num_atoms,
                            'strategy': 'pnictide',
                            'family': f"{ae}-{tm}-{pn}",
                            **exotic_metrics
                        })
    
    # Strategy 4: High-entropy alloys (5 components)
    print("Strategy 4: High-entropy superconductor alloys")
    he_elements = EXOTIC_ELEMENT_GROUPS['high_entropy'][0][:10]
    
    for elements in itertools.combinations(he_elements, 5):
        # Equiatomic composition
        comp_dict = {elem: 1 for elem in elements}
        composition = Composition(comp_dict)
        reduced = composition.reduced_formula
        
        if reduced in seen_formulas:
            continue
        seen_formulas.add(reduced)
        
        tau = tau_valence_mismatch(composition)
        exotic_metrics = assess_exotic_scalability(composition)
        
        if tau < 2.0:  # Should be very good for high-entropy
            candidates.append({
                'formula': composition.formula,
                'reduced_formula': reduced,
                'tau_valence': tau,
                'n_elements': len(composition.elements),
                'n_atoms': composition.num_atoms,
                'strategy': 'high_entropy',
                'family': 'HEA-5',
                **exotic_metrics
            })
    
    print(f"Generated {len(candidates)} exotic candidates")
    
    # Sort by combined exotic score
    def exotic_score(candidate):
        tau_score = np.exp(-candidate['tau_valence'] / 3.0)  # Gentle penalty
        rarity_bonus = 1.0 + (candidate['avg_rarity'] - 2.0) * 0.1  # Bonus for interesting elements
        complexity_bonus = 1.0 + (candidate['n_elements'] - 2.0) * 0.05  # Bonus for complexity
        return tau_score * rarity_bonus * complexity_bonus
    
    candidates.sort(key=exotic_score, reverse=True)
    
    return candidates[:max_candidates]

def main():
    parser = argparse.ArgumentParser(description="Generate exotic superconductor candidates")
    parser.add_argument("--max-atoms", type=int, default=10, help="Maximum atoms per formula")
    parser.add_argument("--max-candidates", type=int, default=1000, help="Maximum candidates to generate")
    parser.add_argument("--out", required=True, help="Output JSON file")
    
    args = parser.parse_args()
    
    # Generate candidates
    candidates = enumerate_exotic_compositions(args.max_atoms, args.max_candidates)
    
    # Save results
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(candidates, f, indent=2)
    
    print(f"\n✅ Saved {len(candidates)} exotic candidates to {output_path}")
    
    # Analyze by strategy
    strategies = {}
    for candidate in candidates:
        strategy = candidate['strategy']
        if strategy not in strategies:
            strategies[strategy] = []
        strategies[strategy].append(candidate)
    
    print("\n📊 Candidates by strategy:")
    for strategy, cands in strategies.items():
        print(f"  {strategy:15s}: {len(cands):4d} candidates")
    
    # Show top candidates
    print("\n🏆 Top 20 exotic superconductor candidates:")
    print("Rank | Formula | τ | Strategy | Rarity | Cost | Exotic Score")
    print("-" * 75)
    
    for i, candidate in enumerate(candidates[:20]):
        tau_score = np.exp(-candidate['tau_valence'] / 3.0)
        rarity_bonus = 1.0 + (candidate['avg_rarity'] - 2.0) * 0.1
        complexity_bonus = 1.0 + (candidate['n_elements'] - 2.0) * 0.05
        exotic_score = tau_score * rarity_bonus * complexity_bonus
        
        print(f"{i+1:4d} | {candidate['reduced_formula']:15s} | "
              f"{candidate['tau_valence']:5.3f} | "
              f"{candidate['strategy']:12s} | "
              f"{candidate['avg_rarity']:6.1f} | "
              f"{candidate['cost_factor']:7.1f} | "
              f"{exotic_score:11.4f}")

if __name__ == "__main__":
    main() 