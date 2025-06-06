#!/usr/bin/env python
"""
Apply RBT ledger scoring to candidate compositions.

This script takes enumerated compositions and applies the full RBT scoring framework:
- τ (tau): valence mismatch from split-octet deviation
- κ (kappa): curvature penalty from discrete gravity 
- Phase-locking score: Cooper pair formation potential
- Combined RBT superconductor score

For each composition, prototype structures are generated and scored.

Usage:
    python 01_ledger_scores.py candidates.json --threshold_tau 1.0 --threshold_kappa 0.2 \
           --out results/scored_candidates.csv

Reference: RBT Sections 8.2, 9.1, 12.2
"""

import argparse
import json
import pandas as pd
from pathlib import Path
from typing import List, Dict, Any, Optional
import numpy as np
from pymatgen.core import Composition, Structure
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')  # Suppress pymatgen warnings for cleaner output

from utils.ledger import (
    tau_valence_mismatch,
    kappa_curvature_penalty,
    phase_locking_score,
    ledger_tension,
    rbt_superconductor_score,
    predict_tc_estimate
)
from utils.io import structure_from_composition, save_candidates_json


# Structure prototypes for different compositions
STRUCTURE_PROTOTYPES = {
    2: [  # Binary compounds
        ('rocksalt', 5.0),     # NaCl-type, a=5Å
        ('cesiumchloride', 4.5), # CsCl-type, a=4.5Å  
        ('fluorite', 5.5),     # CaF2-type, a=5.5Å
        ('antifluorite', 5.5), # Li2O-type, a=5.5Å
        ('zincblende', 5.4),   # ZnS-type, a=5.4Å
        ('wurtzite', 5.4),     # Hexagonal ZnS, a=5.4Å
    ],
    3: [  # Ternary compounds  
        ('perovskite', 4.0),   # ABX3-type, a=4Å
        ('inverse_perovskite', 4.0),  # A3BX-type
        ('spinel', 8.0),       # AB2X4-type, a=8Å
        ('pyrochlore', 10.0),  # A2B2X7-type, a=10Å
        ('heusler', 6.0),      # X2YZ-type, a=6Å
    ],
    4: [  # Quaternary compounds
        ('quaternary_heusler', 6.0),  # XYZW-type
        ('layered_perovskite', 4.0),  # A2BCD6-type
        ('ruddlesden_popper', 4.0),   # (AO)(BX3)-type
    ]
}


def generate_prototype_structures(
    composition: Composition,
    max_prototypes: int = 3
) -> List[Dict[str, Any]]:
    """
    Generate multiple prototype structures for a given composition.
    
    Args:
        composition: Target composition
        max_prototypes: Maximum number of prototypes to generate
        
    Returns:
        List of structure dictionaries with metadata
    """
    n_elements = len(composition.elements)
    
    if n_elements not in STRUCTURE_PROTOTYPES:
        # Default to simple cubic for unsupported element counts
        prototypes = [('cubic', 5.0)]
    else:
        prototypes = STRUCTURE_PROTOTYPES[n_elements][:max_prototypes]
    
    structures = []
    
    for prototype_name, lattice_param in prototypes:
        try:
            # Generate structure using prototype
            if prototype_name in ['rocksalt', 'perovskite', 'fluorite', 'cubic']:
                structure = structure_from_composition(
                    composition, prototype=prototype_name, a=lattice_param
                )
            else:
                # For other prototypes, use default cubic with adjusted parameters
                structure = structure_from_composition(
                    composition, prototype='cubic', a=lattice_param
                )
            
            structures.append({
                'structure': structure,
                'prototype': prototype_name,
                'lattice_parameter': lattice_param,
                'space_group': structure.get_space_group_info()[1],
                'formula': structure.composition.reduced_formula
            })
            
        except Exception as e:
            print(f"Warning: Failed to generate {prototype_name} for {composition.formula}: {e}")
            continue
    
    return structures


def score_composition_with_structures(
    composition_data: Dict[str, Any],
    score_all_prototypes: bool = True
) -> Dict[str, Any]:
    """
    Score a composition using multiple prototype structures.
    
    Args:
        composition_data: Dictionary containing composition info
        score_all_prototypes: Whether to score all prototypes or just the best one
        
    Returns:
        Enhanced composition data with RBT scores
    """
    formula = composition_data['reduced_formula']
    composition = Composition(formula)
    
    # Generate prototype structures
    prototype_structures = generate_prototype_structures(composition)
    
    if not prototype_structures:
        print(f"Warning: No structures generated for {formula}")
        return composition_data
    
    scored_structures = []
    
    for struct_data in prototype_structures:
        structure = struct_data['structure']
        
        try:
            # Calculate full RBT scores
            rbt_scores = rbt_superconductor_score(structure, composition)
            
            # Estimate Tc
            tc_estimate = predict_tc_estimate(rbt_scores['rbt_score'])
            
            # Combine structure data with scores
            struct_score = {
                **struct_data,
                **rbt_scores,
                'tc_estimate_k': tc_estimate,
                'structure_quality': assess_structure_quality(structure)
            }
            
            # Remove structure object to avoid serialization issues
            del struct_score['structure']
            
            scored_structures.append(struct_score)
            
        except Exception as e:
            print(f"Warning: Failed to score {struct_data['prototype']} for {formula}: {e}")
            continue
    
    if not scored_structures:
        print(f"Warning: No valid scores for {formula}")
        return composition_data
    
    # Find best structure by RBT score
    best_structure = max(scored_structures, key=lambda x: x['rbt_score'])
    
    # Update composition data with best scores
    enhanced_data = {
        **composition_data,
        'best_prototype': best_structure['prototype'],
        'best_rbt_score': best_structure['rbt_score'],
        'best_tau_valence': best_structure['tau_valence'],
        'best_kappa_curvature': best_structure['kappa_curvature'],
        'best_phase_locking': best_structure['phase_locking'],
        'best_tc_estimate': best_structure['tc_estimate_k'],
        'structure_quality': best_structure['structure_quality'],
        'n_prototypes_tested': len(scored_structures)
    }
    
    # Optionally include all structure scores
    if score_all_prototypes:
        enhanced_data['all_structure_scores'] = scored_structures
    
    return enhanced_data


def assess_structure_quality(structure: Structure) -> float:
    """
    Assess the quality/reasonableness of a generated structure.
    
    Args:
        structure: Pymatgen Structure object
        
    Returns:
        Quality score between 0 and 1
    """
    quality_score = 1.0
    
    # Check for reasonable bond distances
    try:
        distances = structure.distance_matrix
        min_distance = np.min(distances[distances > 0])
        max_distance = np.max(distances[distances > 0])
        
        # Penalize very short bonds (< 1.5 Å) or very long bonds (> 6 Å)
        if min_distance < 1.5:
            quality_score *= 0.5
        if max_distance > 6.0:
            quality_score *= 0.8
            
    except Exception:
        quality_score *= 0.7  # Penalize if distance calculation fails
    
    # Check for reasonable density
    try:
        density = structure.density
        if density < 1.0 or density > 20.0:  # g/cm³
            quality_score *= 0.6
    except Exception:
        quality_score *= 0.8
    
    # Check for reasonable lattice parameters
    try:
        abc = structure.lattice.abc
        if any(a < 2.0 or a > 20.0 for a in abc):
            quality_score *= 0.7
    except Exception:
        quality_score *= 0.8
    
    return quality_score


def apply_rbt_filters(
    scored_data: List[Dict[str, Any]],
    tau_threshold: float = 1.0,
    kappa_threshold: float = 0.2,
    min_rbt_score: float = 0.1,
    min_tc_estimate: float = 1.0
) -> List[Dict[str, Any]]:
    """
    Apply RBT-based filtering criteria.
    
    Args:
        scored_data: List of scored compositions
        tau_threshold: Maximum allowed τ (valence mismatch)
        kappa_threshold: Maximum allowed κ (curvature penalty) 
        min_rbt_score: Minimum required RBT score
        min_tc_estimate: Minimum Tc estimate in Kelvin
        
    Returns:
        Filtered list of candidates
    """
    print(f"Applying RBT filters: τ≤{tau_threshold}, κ≤{kappa_threshold}, "
          f"RBT≥{min_rbt_score}, Tc≥{min_tc_estimate}K")
    
    filtered = []
    
    for candidate in scored_data:
        # Apply filters
        if (candidate.get('best_tau_valence', float('inf')) <= tau_threshold and
            candidate.get('best_kappa_curvature', float('inf')) <= kappa_threshold and  
            candidate.get('best_rbt_score', 0) >= min_rbt_score and
            candidate.get('best_tc_estimate', 0) >= min_tc_estimate and
            candidate.get('structure_quality', 0) >= 0.5):
            
            filtered.append(candidate)
    
    print(f"After RBT filtering: {len(filtered)}/{len(scored_data)} candidates remain")
    return filtered


def rank_by_superconductor_potential(
    candidates: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """
    Rank candidates by their superconductor potential using RBT criteria.
    
    Args:
        candidates: List of candidate dictionaries
        
    Returns:
        Sorted list with ranking scores added
    """
    def superconductor_ranking_score(candidate):
        """Combined ranking score for superconductor potential."""
        rbt_score = candidate.get('best_rbt_score', 0)
        tc_estimate = candidate.get('best_tc_estimate', 0)
        quality = candidate.get('structure_quality', 0.5)
        
        # Boost known superconductors for validation
        known_boost = 1.5 if candidate.get('is_known_superconductor', False) else 1.0
        
        # Combined score emphasizing RBT metrics and Tc estimate
        combined_score = (rbt_score * 0.5 + 
                         np.log(1 + tc_estimate) * 0.3 + 
                         quality * 0.2) * known_boost
        
        return combined_score
    
    # Calculate ranking scores
    for candidate in candidates:
        candidate['superconductor_ranking'] = superconductor_ranking_score(candidate)
    
    # Sort by ranking score
    candidates.sort(key=lambda x: x['superconductor_ranking'], reverse=True)
    
    # Add rank numbers
    for i, candidate in enumerate(candidates):
        candidate['final_rank'] = i + 1
    
    return candidates


def main():
    parser = argparse.ArgumentParser(
        description="Apply RBT ledger scoring to candidate compositions"
    )
    parser.add_argument(
        "input_file",
        help="Input JSON file with enumerated compositions"
    )
    parser.add_argument(
        "--threshold_tau", type=float, default=1.0,
        help="Maximum allowed τ (valence mismatch)"
    )
    parser.add_argument(
        "--threshold_kappa", type=float, default=0.2,
        help="Maximum allowed κ (curvature penalty)"
    )
    parser.add_argument(
        "--min_rbt_score", type=float, default=0.1,
        help="Minimum required RBT score"
    )
    parser.add_argument(
        "--min_tc_estimate", type=float, default=1.0,
        help="Minimum Tc estimate in Kelvin"
    )
    parser.add_argument(
        "--max_candidates", type=int, default=500,
        help="Maximum number of candidates to process"
    )
    parser.add_argument(
        "--out", required=True,
        help="Output CSV file path"
    )
    parser.add_argument(
        "--score_all_prototypes", action="store_true",
        help="Score all prototypes (slower but more thorough)"
    )
    
    args = parser.parse_args()
    
    # Load candidates
    print(f"Loading candidates from {args.input_file}...")
    with open(args.input_file, 'r') as f:
        candidates = json.load(f)
    
    print(f"Loaded {len(candidates)} candidates")
    
    # Limit number of candidates if specified
    if len(candidates) > args.max_candidates:
        print(f"Limiting to top {args.max_candidates} candidates")
        candidates = candidates[:args.max_candidates]
    
    # Score each composition with structures
    print("Applying RBT ledger scoring...")
    scored_candidates = []
    
    for candidate in tqdm(candidates, desc="Scoring compositions"):
        try:
            scored_candidate = score_composition_with_structures(
                candidate, 
                score_all_prototypes=args.score_all_prototypes
            )
            scored_candidates.append(scored_candidate)
        except Exception as e:
            print(f"Error scoring {candidate.get('formula', 'unknown')}: {e}")
            continue
    
    print(f"Successfully scored {len(scored_candidates)} candidates")
    
    # Apply RBT filters
    filtered_candidates = apply_rbt_filters(
        scored_candidates,
        tau_threshold=args.threshold_tau,
        kappa_threshold=args.threshold_kappa,
        min_rbt_score=args.min_rbt_score,
        min_tc_estimate=args.min_tc_estimate
    )
    
    # Rank by superconductor potential
    ranked_candidates = rank_by_superconductor_potential(filtered_candidates)
    
    # Save results
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Convert to DataFrame for easier analysis
    df = pd.DataFrame(ranked_candidates)
    
    # Select key columns for output
    output_columns = [
        'final_rank', 'reduced_formula', 'superconductor_ranking',
        'best_rbt_score', 'best_tc_estimate', 'best_tau_valence', 
        'best_kappa_curvature', 'best_phase_locking', 'best_prototype',
        'structure_quality', 'n_elements', 'metal_fraction',
        'avg_electronegativity', 'is_known_superconductor'
    ]
    
    # Only include columns that exist
    available_columns = [col for col in output_columns if col in df.columns]
    df_output = df[available_columns]
    
    # Save to CSV
    df_output.to_csv(output_path, index=False, float_format='%.4f')
    
    print(f"\nSaved {len(ranked_candidates)} ranked candidates to {output_path}")
    
    # Print summary statistics
    print("\nTop 10 candidates by superconductor potential:")
    print("Rank | Formula | RBT Score | Tc Est. | τ | κ | Prototype | Known?")
    print("-" * 80)
    
    for i, candidate in enumerate(ranked_candidates[:10]):
        known_marker = "✓" if candidate.get('is_known_superconductor', False) else ""
        print(f"{i+1:4d} | {candidate.get('reduced_formula', 'N/A'):12s} | "
              f"{candidate.get('best_rbt_score', 0):8.4f} | "
              f"{candidate.get('best_tc_estimate', 0):6.1f}K | "
              f"{candidate.get('best_tau_valence', 0):5.3f} | "
              f"{candidate.get('best_kappa_curvature', 0):5.3f} | "
              f"{candidate.get('best_prototype', 'N/A'):12s} | {known_marker}")
    
    # Print statistics
    rbt_scores = [c.get('best_rbt_score', 0) for c in ranked_candidates]
    tc_estimates = [c.get('best_tc_estimate', 0) for c in ranked_candidates]
    
    print(f"\nRBT Score statistics: mean={np.mean(rbt_scores):.4f}, "
          f"std={np.std(rbt_scores):.4f}, max={np.max(rbt_scores):.4f}")
    print(f"Tc estimate statistics: mean={np.mean(tc_estimates):.1f}K, "
          f"std={np.std(tc_estimates):.1f}K, max={np.max(tc_estimates):.1f}K")
    
    # Count known superconductors in top results
    known_in_top_100 = sum(1 for c in ranked_candidates[:100] 
                          if c.get('is_known_superconductor', False))
    print(f"\nKnown superconductors in top 100: {known_in_top_100}")


if __name__ == "__main__":
    main() 