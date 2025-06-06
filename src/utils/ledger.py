"""
Recursive Becoming Theory (RBT) ledger mathematics for superconductor discovery.

Based on the axiom that reality emerges from recursive bit-flips,
this module implements the key RBT constructs:
- Split-octet valence balance (τ metric)
- Curvature penalty (κ metric) 
- Phase-locking assessment for superconducting gaps
- Ledger tension calculations for coupling unification

Reference: Chauhan & Chouhan (2025) "Recursive Becoming: From Nothingness to Everything"
DOI: 10.5281/zenodo.15391360
"""

import itertools
import math
import numpy as np
from typing import Dict, List, Tuple, Union
from pymatgen.core import Element, Structure, Composition
from pymatgen.analysis.local_env import CrystalNN
import warnings


# RBT Constants from the theory
_SPLIT_OCTET = [2, 4, 8, 12]  # Recursive branching points in number tower
_PRIMITIVE_ENERGY = 5.34e7  # GeV, ε₀ from RBT Section 5.3
_MASS_QUANTUM = 5.0  # GeV, m* after gauge sharing (Section 9.1)
_HORIZON_KAPPA = 1.2e-3  # Entropy leakage coefficient (Section 10.4)


def nearest_split(v: int) -> int:
    """
    Return smallest split-octet value ≥ v.
    
    In RBT, the number tower construction creates natural
    branching points at 2, 4, 8, 12 corresponding to 
    ℕ→ℤ→ℚ→ℝ→ℂ→ℍ→𝕆 transitions.
    
    Args:
        v: Valence or group number
        
    Returns:
        Nearest split-octet value
    """
    for s in _SPLIT_OCTET:
        if v <= s:
            return s
    return _SPLIT_OCTET[-1]


def tau_valence_mismatch(composition: Union[Composition, Dict[str, float]]) -> float:
    """
    Calculate ledger-valence mismatch τ (average over atoms).
    
    In RBT Section 12.2, the periodic table emerges from tag-axis
    permutations. Superconducting phases prefer valence-balanced
    configurations where electrons can form Cooper-like pairs.
    
    τ = (1/N) Σᵢ nᵢ(vᵢ - nearest_split(vᵢ))²
    
    Args:
        composition: Pymatgen Composition or dict {element: amount}
        
    Returns:
        Valence mismatch parameter τ
    """
    if isinstance(composition, Composition):
        comp_dict = composition.as_dict()
    else:
        comp_dict = composition
    
    total_weighted_mismatch = 0.0
    total_atoms = 0.0
    
    for element_str, amount in comp_dict.items():
        element = Element(element_str)
        valence = element.group if element.group is not None else element.Z
        split_val = nearest_split(valence)
        mismatch = (valence - split_val) ** 2
        
        total_weighted_mismatch += amount * mismatch
        total_atoms += amount
    
    return total_weighted_mismatch / total_atoms if total_atoms > 0 else 0.0


def kappa_curvature_penalty(structure: Structure, rcut: float = 3.0) -> float:
    """
    Calculate curvature penalty κ from RBT discrete gravity.
    
    From Section 8.2, the ledger metric creates curvature that
    affects superconducting coherence length. High curvature
    disrupts Cooper pair formation.
    
    κ = Σᵢⱼ (Zᵢ² + Zⱼ²) / dᵢⱼ² for pairs within rcut
    
    Args:
        structure: Pymatgen Structure object
        rcut: Cutoff radius in Angstroms
        
    Returns:
        Curvature penalty κ
    """
    atomic_numbers = {site.specie.symbol: site.specie.Z for site in structure}
    penalty = 0.0
    
    for i, j in itertools.combinations(range(len(structure)), 2):
        try:
            distance = structure.get_distance(i, j, mic=True)  # minimum image convention
            if distance < rcut and distance > 0.1:  # Avoid singularities
                zi = atomic_numbers[structure[i].specie.symbol]
                zj = atomic_numbers[structure[j].specie.symbol]
                penalty += (zi**2 + zj**2) / distance**2
        except Exception:
            continue  # Skip problematic pairs
    
    return penalty


def phase_locking_score(structure: Structure, composition: Composition) -> float:
    """
    Estimate phase-locking strength for superconducting gap formation.
    
    From RBT Section 9.1, mass generation occurs through phase-locking
    mechanisms. For superconductors, this manifests as the ability to
    form coherent Cooper pairs across the crystal.
    
    Args:
        structure: Crystal structure
        composition: Chemical composition
        
    Returns:
        Phase-locking score (higher = better for superconductivity)
    """
    # Base score from valence balance
    tau = tau_valence_mismatch(composition)
    valence_score = np.exp(-tau)  # Exponential penalty for mismatch
    
    # Coordination environment contribution
    try:
        cnn = CrystalNN()
        coord_nums = []
        
        for i, site in enumerate(structure):
            try:
                cn_dict = cnn.get_cn_dict(structure, i)
                total_cn = sum(cn_dict.values())
                coord_nums.append(total_cn)
            except Exception:
                coord_nums.append(6)  # Default coordination
        
        # Prefer moderate, uniform coordination (like cuprate superconductors)
        avg_coord = np.mean(coord_nums)
        coord_uniformity = 1.0 / (1.0 + np.std(coord_nums))
        coord_score = np.exp(-abs(avg_coord - 6.0) / 2.0) * coord_uniformity
        
    except Exception:
        coord_score = 0.5  # Default if coordination analysis fails
    
    # Electronic structure hint from composition
    metal_fraction = sum(
        composition[el] for el in composition.elements 
        if el.is_metal
    ) / composition.num_atoms
    
    electronic_score = metal_fraction * (1 - metal_fraction) * 4  # Optimal ~0.5
    
    return valence_score * coord_score * electronic_score


def ledger_tension(composition: Composition, temperature: float = 300.0) -> float:
    """
    Calculate ledger tension that drives coupling unification.
    
    From RBT Section 10.1, one-loop running emerges from ledger tension.
    At the superconducting transition, this tension reorganizes to
    minimize free energy.
    
    Args:
        composition: Chemical composition
        temperature: Temperature in Kelvin
        
    Returns:
        Ledger tension parameter
    """
    # Average atomic mass contribution
    avg_mass = composition.weight / composition.num_atoms
    mass_factor = np.log(avg_mass / _MASS_QUANTUM) if avg_mass > 0 else 0
    
    # Thermal contribution
    thermal_factor = _HORIZON_KAPPA * temperature / 300.0
    
    # Valence contribution
    tau = tau_valence_mismatch(composition)
    valence_factor = np.exp(-tau / 2.0)
    
    return mass_factor * thermal_factor * valence_factor


def rbt_superconductor_score(
    structure: Structure, 
    composition: Composition,
    weights: Tuple[float, float, float] = (0.4, 0.3, 0.3)
) -> Dict[str, float]:
    """
    Comprehensive RBT-based superconductor scoring function.
    
    Combines all RBT metrics to predict superconducting potential:
    - τ: valence balance (lower is better)
    - κ: curvature penalty (lower is better)  
    - Phase-locking strength (higher is better)
    
    Args:
        structure: Crystal structure
        composition: Chemical composition
        weights: (τ_weight, κ_weight, phase_weight) for final score
        
    Returns:
        Dictionary with individual scores and combined RBT score
    """
    # Individual metrics
    tau = tau_valence_mismatch(composition)
    kappa = kappa_curvature_penalty(structure)
    phase_lock = phase_locking_score(structure, composition)
    tension = ledger_tension(composition)
    
    # Normalize scores to [0,1] range
    tau_norm = np.exp(-tau)  # Lower τ is better
    kappa_norm = np.exp(-kappa / 1000.0)  # Lower κ is better
    phase_norm = min(phase_lock, 1.0)  # Higher phase-locking is better
    
    # Combined RBT score
    w_tau, w_kappa, w_phase = weights
    rbt_score = w_tau * tau_norm + w_kappa * kappa_norm + w_phase * phase_norm
    
    return {
        'tau_valence': tau,
        'kappa_curvature': kappa,
        'phase_locking': phase_lock,
        'ledger_tension': tension,
        'tau_normalized': tau_norm,
        'kappa_normalized': kappa_norm,
        'phase_normalized': phase_norm,
        'rbt_score': rbt_score
    }


def predict_tc_estimate(rbt_score: float, base_temp: float = 300.0) -> float:
    """
    Rough Tc estimate from RBT score using empirical scaling.
    
    This is a heuristic mapping that should be calibrated against
    known superconductors. The RBT theory suggests Tc scales with
    the ledger reorganization energy.
    
    Args:
        rbt_score: Combined RBT superconductor score
        base_temp: Reference temperature scale
        
    Returns:
        Estimated critical temperature in Kelvin
    """
    # Empirical scaling - to be refined with training data
    tc_estimate = base_temp * rbt_score**2 * 0.3
    return max(tc_estimate, 0.0)


# Utility functions for enumeration and filtering
def enumerate_valence_balanced_compositions(
    elements: List[str], 
    max_atoms: int = 10,
    max_species: int = 4
) -> List[Composition]:
    """
    Generate compositions with zero net valence mismatch.
    
    This implements the RBT principle that stable phases
    prefer split-octet balanced configurations.
    
    Args:
        elements: List of element symbols to consider
        max_atoms: Maximum total atoms per formula unit
        max_species: Maximum number of different elements
        
    Returns:
        List of valence-balanced compositions
    """
    balanced_compositions = []
    
    for n_species in range(2, min(max_species + 1, len(elements) + 1)):
        for elem_combo in itertools.combinations(elements, n_species):
            # Try different stoichiometries
            for stoich in itertools.product(range(1, max_atoms), repeat=n_species):
                if sum(stoich) > max_atoms:
                    continue
                
                # Create composition
                comp_dict = dict(zip(elem_combo, stoich))
                comp = Composition(comp_dict)
                
                # Check valence balance
                if tau_valence_mismatch(comp) < 0.1:  # Nearly balanced
                    balanced_compositions.append(comp)
    
    return balanced_compositions


# Export key functions
__all__ = [
    'tau_valence_mismatch',
    'kappa_curvature_penalty', 
    'phase_locking_score',
    'ledger_tension',
    'rbt_superconductor_score',
    'predict_tc_estimate',
    'enumerate_valence_balanced_compositions'
] 