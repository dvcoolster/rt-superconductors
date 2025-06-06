#!/usr/bin/env python
"""
Generate valence-balanced binary/ternary/quaternary compositions for superconductor discovery.

This script implements the RBT principle that stable superconducting phases prefer
split-octet balanced configurations (τ ≈ 0). The enumeration focuses on compositions
where valence electrons can form Cooper-like pairs efficiently.

Usage:
    python 00_enumerate.py --max-atoms 10 --out data/generated/candidates.json
    
Reference: RBT Section 12.2 - "Tag permutations and the periodic table"
"""

import argparse
import json
import itertools
from pathlib import Path
from typing import List, Dict, Any, Set
import numpy as np
from pymatgen.core import Composition, Element
from tqdm import tqdm

from utils.ledger import (
    tau_valence_mismatch, 
    nearest_split, 
    enumerate_valence_balanced_compositions
)
from utils.io import save_candidates_json


# RBT-optimized element selection for superconductors
# Focus on elements that appear in known superconductors + RBT-favorable split-octet elements
RBT_ELEMENTS = [
    # Light metals (Group 1-2)
    "H", "Li", "Be", "Na", "Mg", "K", "Ca", "Sr", "Ba",
    
    # Transition metals (common in superconductors)
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    
    # p-block elements (important for electronic structure)
    "B", "C", "N", "O", "F", 
    "Al", "Si", "P", "S", "Cl",
    "Ga", "Ge", "As", "Se", "Br",
    "In", "Sn", "Sb", "Te", "I",
    "Tl", "Pb", "Bi",
    
    # Rare earths (for exotic superconductors)
    "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"
]


def calculate_composition_metrics(composition: Composition) -> Dict[str, float]:
    """
    Calculate RBT metrics for a composition without structure.
    
    Args:
        composition: Pymatgen Composition object
        
    Returns:
        Dictionary of calculated metrics
    """
    # Valence balance (τ metric)
    tau = tau_valence_mismatch(composition)
    
    # Electronic structure hints
    metal_count = sum(1 for el in composition.elements if el.is_metal)
    nonmetal_count = len(composition.elements) - metal_count
    
    # Average electronegativity and its spread
    electronegativities = [el.X for el in composition.elements if el.X is not None]
    if electronegativities:
        avg_electronegativity = np.mean(electronegativities)
        electronegativity_spread = np.std(electronegativities)
    else:
        avg_electronegativity = 2.0  # Default
        electronegativity_spread = 0.0
    
    # Average atomic mass
    avg_mass = composition.weight / composition.num_atoms
    
    # Density estimate (crude)
    avg_radius = np.mean([el.atomic_radius or 1.5 for el in composition.elements])
    density_estimate = avg_mass / (4/3 * np.pi * avg_radius**3) if avg_radius > 0 else 1.0
    
    return {
        'formula': composition.formula,
        'tau_valence': tau,
        'n_elements': len(composition.elements),
        'n_atoms': composition.num_atoms,
        'metal_fraction': metal_count / len(composition.elements),
        'avg_electronegativity': avg_electronegativity,
        'electronegativity_spread': electronegativity_spread,
        'avg_mass': avg_mass,
        'density_estimate': density_estimate,
        'reduced_formula': composition.reduced_formula
    }


def enumerate_rbt_compositions(
    elements: List[str],
    max_atoms: int = 10,
    max_species: int = 4,
    tau_threshold: float = 0.5
) -> List[Dict[str, Any]]:
    """
    Enumerate RBT-favorable compositions with valence balance.
    
    Args:
        elements: List of element symbols to consider
        max_atoms: Maximum atoms per formula unit
        max_species: Maximum number of different elements
        tau_threshold: Maximum allowed τ (valence mismatch)
        
    Returns:
        List of composition dictionaries with metrics
    """
    print(f"Enumerating compositions from {len(elements)} elements...")
    print(f"Constraints: max_atoms={max_atoms}, max_species={max_species}, τ≤{tau_threshold}")
    
    candidates = []
    compositions_seen: Set[str] = set()
    
    # Binary compositions (A-B)
    print("Generating binary compositions...")
    for elem_pair in tqdm(itertools.combinations(elements, 2)):
        for stoich in itertools.product(range(1, max_atoms), repeat=2):
            if sum(stoich) > max_atoms:
                continue
                
            comp_dict = dict(zip(elem_pair, stoich))
            composition = Composition(comp_dict)
            
            # Skip if we've seen this reduced formula
            reduced = composition.reduced_formula
            if reduced in compositions_seen:
                continue
            compositions_seen.add(reduced)
            
            # Calculate metrics
            metrics = calculate_composition_metrics(composition)
            
            # Apply RBT filter
            if metrics['tau_valence'] <= tau_threshold:
                candidates.append(metrics)
    
    # Ternary compositions (A-B-C)
    if max_species >= 3:
        print("Generating ternary compositions...")
        for elem_triple in tqdm(itertools.combinations(elements, 3)):
            for stoich in itertools.product(range(1, max_atoms), repeat=3):
                if sum(stoich) > max_atoms:
                    continue
                    
                comp_dict = dict(zip(elem_triple, stoich))
                composition = Composition(comp_dict)
                
                reduced = composition.reduced_formula
                if reduced in compositions_seen:
                    continue
                compositions_seen.add(reduced)
                
                metrics = calculate_composition_metrics(composition)
                
                if metrics['tau_valence'] <= tau_threshold:
                    candidates.append(metrics)
    
    # Quaternary compositions (A-B-C-D)
    if max_species >= 4:
        print("Generating quaternary compositions...")
        # Limit quaternary search to avoid combinatorial explosion
        selected_elements = elements[:30]  # Top 30 elements for quaternary
        
        for elem_quad in tqdm(itertools.combinations(selected_elements, 4)):
            for stoich in itertools.product(range(1, min(max_atoms//2, 4)), repeat=4):
                if sum(stoich) > max_atoms:
                    continue
                    
                comp_dict = dict(zip(elem_quad, stoich))
                composition = Composition(comp_dict)
                
                reduced = composition.reduced_formula
                if reduced in compositions_seen:
                    continue
                compositions_seen.add(reduced)
                
                metrics = calculate_composition_metrics(composition)
                
                if metrics['tau_valence'] <= tau_threshold:
                    candidates.append(metrics)
    
    print(f"Generated {len(candidates)} RBT-favorable compositions")
    return candidates


def add_known_superconductors(candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Add known high-Tc superconductors for comparison and training.
    
    Args:
        candidates: List of generated candidates
        
    Returns:
        Extended list including known superconductors
    """
    known_superconductors = [
        "YBa2Cu3O7",      # YBCO - 92K
        "Bi2Sr2CaCu2O8",  # BSCCO - 85K  
        "Tl2Ba2CaCu2O8",  # TBCCO - 110K
        "HgBa2Ca2Cu3O8",  # HBCCO - 133K
        "LaFeAsO",        # Iron pnictide - 26K
        "Ba0.6K0.4Fe2As2", # 122 Iron pnictide - 38K
        "MgB2",           # Magnesium diboride - 39K
        "H3S",            # Hydrogen sulfide - 203K (high P)
        "LaH10",          # Lanthanum hydride - 250K (high P)
        "NbTi",           # Niobium titanium - 10K
        "Nb3Sn",          # Niobium tin - 18K
        "V3Si",           # Vanadium silicon - 17K
        "PbMo6S8",        # Chevrel phase - 15K
        "CeCoIn5",        # Heavy fermion - 2.3K
        "Sr2RuO4",        # Ruthenate - 1.5K
        "UPt3",           # Actinide - 0.5K
        "CaC6",           # Graphite intercalate - 11.5K
        "Li0.68Mo6O17",   # Purple bronze - 1.9K
    ]
    
    print(f"Adding {len(known_superconductors)} known superconductors...")
    
    for formula in known_superconductors:
        try:
            composition = Composition(formula)
            metrics = calculate_composition_metrics(composition)
            metrics['is_known_superconductor'] = True
            metrics['source'] = 'literature'
            candidates.append(metrics)
        except Exception as e:
            print(f"Warning: Could not parse {formula}: {e}")
    
    return candidates


def filter_and_rank_candidates(
    candidates: List[Dict[str, Any]],
    top_n: int = 1000
) -> List[Dict[str, Any]]:
    """
    Apply additional filtering and rank by RBT favorability.
    
    Args:
        candidates: List of candidate compositions
        top_n: Number of top candidates to return
        
    Returns:
        Filtered and ranked candidates
    """
    print(f"Filtering and ranking {len(candidates)} candidates...")
    
    # Apply additional filters
    filtered = []
    for candidate in candidates:
        # Skip if too simple or too complex
        if candidate['n_elements'] < 2 or candidate['n_atoms'] > 20:
            continue
        
        # Require some metallic character for superconductivity
        if candidate['metal_fraction'] < 0.3:
            continue
        
        # Reasonable electronegativity spread (not too ionic)
        if candidate['electronegativity_spread'] > 2.0:
            continue
            
        # Reasonable density estimate
        if candidate['density_estimate'] > 50.0 or candidate['density_estimate'] < 0.1:
            continue
        
        filtered.append(candidate)
    
    print(f"After filtering: {len(filtered)} candidates remain")
    
    # Rank by combined RBT score
    def rbt_ranking_score(candidate):
        """Combined ranking score favoring RBT principles."""
        tau_score = np.exp(-candidate['tau_valence'])  # Lower τ is better
        
        # Favor moderate electronegativity (like cuprates)
        eneg_score = np.exp(-abs(candidate['avg_electronegativity'] - 2.5))
        
        # Favor moderate complexity
        complexity_score = np.exp(-abs(candidate['n_elements'] - 3))
        
        # Boost known superconductors for validation
        known_boost = 2.0 if candidate.get('is_known_superconductor', False) else 1.0
        
        return tau_score * eneg_score * complexity_score * known_boost
    
    # Sort by ranking score
    filtered.sort(key=rbt_ranking_score, reverse=True)
    
    # Add ranking scores to candidates
    for i, candidate in enumerate(filtered):
        candidate['rbt_ranking_score'] = rbt_ranking_score(candidate)
        candidate['rank'] = i + 1
    
    return filtered[:top_n]


def main():
    parser = argparse.ArgumentParser(
        description="Generate RBT-favorable compositions for superconductor discovery"
    )
    parser.add_argument(
        "--max-atoms", type=int, default=10,
        help="Maximum atoms per formula unit"
    )
    parser.add_argument(
        "--max-species", type=int, default=4,
        help="Maximum number of different elements"
    )
    parser.add_argument(
        "--tau-threshold", type=float, default=0.5,
        help="Maximum allowed τ (valence mismatch)"
    )
    parser.add_argument(
        "--elements", nargs="+", default=RBT_ELEMENTS,
        help="List of elements to consider"
    )
    parser.add_argument(
        "--top-n", type=int, default=1000,
        help="Number of top candidates to save"
    )
    parser.add_argument(
        "--out", required=True,
        help="Output JSON file path"
    )
    parser.add_argument(
        "--include-known", action="store_true",
        help="Include known superconductors for comparison"
    )
    
    args = parser.parse_args()
    
    # Generate compositions
    candidates = enumerate_rbt_compositions(
        elements=args.elements,
        max_atoms=args.max_atoms,
        max_species=args.max_species,
        tau_threshold=args.tau_threshold
    )
    
    # Add known superconductors if requested
    if args.include_known:
        candidates = add_known_superconductors(candidates)
    
    # Filter and rank
    top_candidates = filter_and_rank_candidates(candidates, args.top_n)
    
    # Save results
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(top_candidates, f, indent=2)
    
    print(f"\nSaved {len(top_candidates)} top candidates to {output_path}")
    
    # Print summary statistics
    print("\nTop 10 candidates by RBT score:")
    print("Rank | Formula | τ | Elements | RBT Score | Known?")
    print("-" * 60)
    
    for i, candidate in enumerate(top_candidates[:10]):
        known_marker = "✓" if candidate.get('is_known_superconductor', False) else ""
        print(f"{i+1:4d} | {candidate['reduced_formula']:15s} | "
              f"{candidate['tau_valence']:5.3f} | "
              f"{candidate['n_elements']:2d} | "
              f"{candidate['rbt_ranking_score']:8.4f} | {known_marker}")
    
    # Statistics
    tau_values = [c['tau_valence'] for c in top_candidates]
    print(f"\nτ statistics: mean={np.mean(tau_values):.3f}, "
          f"std={np.std(tau_values):.3f}, "
          f"min={np.min(tau_values):.3f}, "
          f"max={np.max(tau_values):.3f}")


if __name__ == "__main__":
    main() 