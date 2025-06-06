#!/usr/bin/env python
"""Test practical superconductor compositions with common elements."""

import sys
sys.path.append('src')

from pymatgen.core import Composition
from utils.ledger import tau_valence_mismatch

# Test practical compositions with common elements
practical_formulas = [
    'MgB2',      # Known superconductor
    'LiTi2', 
    'CaFe2As2',  # Iron pnictide family
    'LaFeO3', 
    'SrTiO3', 
    'BaTiO3',
    'LiFePO4', 
    'NaCoO2',
    'LiCoO2',
    'FeSe',      # Known superconductor
    'CaC6',      # Known intercalate
    'Li2CuO2',
    'Ba2CuO3',
    'Sr2CuO3',
    'Ca2CuO3',
    'La2CuO4',   # Cuprate family
    'YBa2Cu3O7', # Known high-Tc
    'AlB2',      # Known structure type
    'TiN',       # Simple binary
    'NbN',       # Known superconductor
    'MoS2',      # Dichalcogenide
    'TaS2',      # Known superconductor
    'CrSi2',
    'FeSi2',
    'CoSi2',
    'NiSi2'
]

print("Testing practical superconductor candidates (common, scalable elements):")
print("Formula       τ (valence)  RBT Assessment")
print("-" * 50)

good_candidates = []

for formula in practical_formulas:
    try:
        comp = Composition(formula)
        tau = tau_valence_mismatch(comp)
        
        # RBT assessment
        if tau < 0.1:
            assessment = "EXCELLENT"
        elif tau < 0.25:
            assessment = "GOOD"
        elif tau < 0.5:
            assessment = "MODERATE"
        else:
            assessment = "POOR"
        
        print(f"{formula:12s}  {tau:8.4f}     {assessment}")
        
        if tau < 0.3:  # Collect good candidates
            good_candidates.append((formula, tau))
            
    except Exception as e:
        print(f"{formula:12s}  ERROR: {str(e)[:20]}")

print(f"\n🎯 Top candidates with τ < 0.3:")
good_candidates.sort(key=lambda x: x[1])
for formula, tau in good_candidates[:10]:
    print(f"  {formula}: τ = {tau:.4f}") 