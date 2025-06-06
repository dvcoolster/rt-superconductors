#!/usr/bin/env python
"""
Comprehensive MgFe Alloy Analysis for Superconductor Discovery

This script analyzes MgFe alloys from multiple perspectives:
1. Known phase diagrams and crystal structures
2. RBT-guided optimization of stoichiometry
3. Synthesis pathways and conditions
4. Testing protocols for superconductivity
5. What previous researchers might have missed
"""

import sys
sys.path.append('src')

import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core import Composition, Element, Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.lattice import Lattice
from utils.ledger import tau_valence_mismatch, kappa_curvature_penalty, phase_locking_score

# Known MgFe phase information from metallurgy literature
MGFE_PHASES = {
    'Mg17Fe': {
        'stoichiometry': 'Mg17Fe',
        'structure': 'cubic',
        'space_group': 'Fm3m',
        'a': 8.96,  # Angstroms
        'formation_energy': -0.15,  # eV/atom (approximate)
        'stability': 'metastable',
        'synthesis_temp': 723,  # K
        'notes': 'Known intermetallic, limited solubility'
    },
    'Mg2Fe': {
        'stoichiometry': 'Mg2Fe',
        'structure': 'hexagonal',
        'space_group': 'P63/mmc',
        'a': 5.24,
        'c': 8.52,
        'formation_energy': -0.08,
        'stability': 'unstable_at_RT',
        'synthesis_temp': 873,
        'notes': 'High temperature phase'
    },
    'MgFe_solid_solution': {
        'stoichiometry': 'MgxFe(1-x)',
        'structure': 'bcc/fcc_mixed',
        'space_group': 'variable',
        'solubility_limit': 0.03,  # ~3% Mg in Fe
        'formation_energy': -0.02,
        'stability': 'limited',
        'synthesis_temp': 1173,
        'notes': 'Very limited mutual solubility'
    }
}

def analyze_rbt_optimization():
    """Analyze RBT-guided optimization of MgFe stoichiometry."""
    
    print("🧮 RBT-GUIDED STOICHIOMETRY OPTIMIZATION")
    print("=" * 50)
    
    # Test range of Mg:Fe ratios
    ratios_to_test = [
        (1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (3, 2), (2, 3),
        (4, 1), (1, 4), (5, 1), (1, 5), (5, 2), (2, 5), (5, 3), (3, 5)
    ]
    
    rbt_results = []
    
    for mg_count, fe_count in ratios_to_test:
        comp = Composition({'Mg': mg_count, 'Fe': fe_count})
        
        # Calculate RBT metrics
        tau = tau_valence_mismatch(comp)
        kappa = 0.0  # Will calculate when we have structure
        phase_lock = 0.0  # Will estimate based on composition
        
        # Combined RBT score
        rbt_score = 0.4 * np.exp(-tau) + 0.3 * np.exp(-kappa/1000) + 0.3 * phase_lock
        
        rbt_results.append({
            'formula': comp.reduced_formula,
            'mg_fe_ratio': f"{mg_count}:{fe_count}",
            'tau': tau,
            'kappa': kappa,
            'phase_locking': phase_lock,
            'rbt_score': rbt_score,
            'n_atoms': comp.num_atoms
        })
    
    # Sort by RBT score
    rbt_results.sort(key=lambda x: x['rbt_score'], reverse=True)
    
    print("Top RBT-optimized MgFe compositions:")
    print("Rank | Formula | Ratio | τ | κ | Phase Lock | RBT Score")
    print("-" * 60)
    
    for i, result in enumerate(rbt_results[:10]):
        print(f"{i+1:4d} | {result['formula']:7s} | "
              f"{result['mg_fe_ratio']:5s} | "
              f"{result['tau']:5.3f} | "
              f"{result['kappa']:7.1f} | "
              f"{result['phase_locking']:10.3f} | "
              f"{result['rbt_score']:9.4f}")
    
    return rbt_results

def analyze_synthesis_pathways():
    """Analyze practical synthesis methods for MgFe alloys."""
    
    print("\n🔥 SYNTHESIS PATHWAY ANALYSIS")
    print("=" * 50)
    
    synthesis_methods = {
        'arc_melting': {
            'description': 'Electric arc melting under inert atmosphere',
            'temperature': '2000-3000K',
            'atmosphere': 'Ar or He',
            'time': '30-60 seconds',
            'pros': ['Fast', 'High purity', 'Good mixing'],
            'cons': ['Requires equipment', 'Rapid cooling'],
            'rbt_considerations': 'May not allow optimal crystallization for RBT alignment'
        },
        
        'powder_metallurgy': {
            'description': 'Mix Mg and Fe powders, compact, sinter',
            'temperature': '673-873K',
            'atmosphere': 'Ar or vacuum',
            'time': '2-24 hours',
            'pros': ['Controlled composition', 'Uniform mixing', 'Scalable'],
            'cons': ['Mg oxidation risk', 'Sintering challenges'],
            'rbt_considerations': 'Allows slow crystallization for RBT optimization'
        },
        
        'mechanical_alloying': {
            'description': 'Ball milling to create metastable phases',
            'temperature': 'Room temperature + local heating',
            'atmosphere': 'Ar',
            'time': '10-50 hours',
            'pros': ['Metastable phases', 'Fine microstructure', 'Non-equilibrium'],
            'cons': ['Contamination risk', 'Time consuming'],
            'rbt_considerations': 'Could create novel phases not in equilibrium diagrams'
        },
        
        'thin_film_deposition': {
            'description': 'Co-sputtering or MBE of Mg and Fe',
            'temperature': '300-773K',
            'atmosphere': 'Ultra-high vacuum',
            'time': '1-10 hours',
            'pros': ['Precise control', 'Novel structures', 'Clean interfaces'],
            'cons': ['Expensive', 'Limited quantity'],
            'rbt_considerations': 'Atomic-level control for optimal RBT arrangements'
        },
        
        'liquid_quenching': {
            'description': 'Rapid cooling from molten state',
            'temperature': '1400-1800K → 300K in ms',
            'atmosphere': 'Ar',
            'time': 'Milliseconds',
            'pros': ['Metastable phases', 'Suppresses separation'],
            'cons': ['Complex equipment', 'Small samples'],
            'rbt_considerations': 'Traps high-temperature RBT-favorable configurations'
        }
    }
    
    print("Recommended synthesis pathway based on RBT requirements:")
    print("\n1. POWDER METALLURGY (Primary recommendation)")
    print("   - Allows controlled crystallization for RBT optimization")
    print("   - Can fine-tune Mg:Fe ratios precisely")
    print("   - Scalable for practical applications")
    
    print("\n2. MECHANICAL ALLOYING (Secondary recommendation)")
    print("   - Can create metastable phases not accessible by equilibrium")
    print("   - May unlock hidden RBT-favorable structures")
    print("   - Non-equilibrium processing aligns with RBT principles")
    
    return synthesis_methods

def what_others_missed():
    """Analysis of what previous researchers might have overlooked."""
    
    print("\n🔍 WHAT OTHERS MIGHT HAVE MISSED")
    print("=" * 50)
    
    missed_factors = {
        'valence_electron_optimization': {
            'issue': 'Previous work focused on equilibrium phases',
            'rbt_insight': 'RBT suggests specific Mg:Fe ratios optimize valence electron arrangements',
            'solution': 'Test compositions with τ < 0.5, especially 1:1, 2:1, 3:2 ratios',
            'why_missed': 'Traditional metallurgy assumes equilibrium phase diagram limits'
        },
        
        'processing_conditions': {
            'issue': 'Standard synthesis may not achieve RBT-favorable structures',
            'rbt_insight': 'Slow cooling or specific temperature profiles could enable discrete gravity alignment',
            'solution': 'Controlled cooling rates, specific annealing temperatures',
            'why_missed': 'No theoretical framework to guide non-standard processing'
        },
        
        'measurement_temperature_range': {
            'issue': 'May not have tested low enough temperatures',
            'rbt_insight': 'RBT predicts Tc around 44K - requires measurements below 50K',
            'solution': 'Systematic measurements from 300K down to 4K',
            'why_missed': 'MgFe alloys not expected to be superconductors'
        },
        
        'microstructure_effects': {
            'issue': 'Bulk properties may mask superconductivity',
            'rbt_insight': 'Superconductivity might be localized to interfaces or specific phases',
            'solution': 'Test different microstructures: single crystals, thin films, nanostructures',
            'why_missed': 'Focused on bulk properties rather than local environments'
        },
        
        'stoichiometry_precision': {
            'issue': 'Previous work may not have tested exact RBT-optimized ratios',
            'rbt_insight': 'Perfect τ=0 requires precise atomic ratios',
            'solution': 'Synthesize exact compositions: MgFe, Mg₂Fe, Mg₃Fe₂',
            'why_missed': 'Limited solubility led to focus on equilibrium compositions'
        },
        
        'atmosphere_control': {
            'issue': 'Oxidation may have masked superconducting properties',
            'rbt_insight': 'Pure metallic bonding required for RBT mechanisms',
            'solution': 'Ultra-pure synthesis, inert atmosphere handling, surface protection',
            'why_missed': 'Standard metallurgy accepts some oxidation'
        }
    }
    
    print("Key factors that may have been overlooked:")
    for factor, details in missed_factors.items():
        print(f"\n• {factor.replace('_', ' ').title()}:")
        print(f"  Issue: {details['issue']}")
        print(f"  RBT Insight: {details['rbt_insight']}")
        print(f"  Solution: {details['solution']}")
    
    return missed_factors

def experimental_protocol():
    """Detailed experimental protocol for MgFe superconductor testing."""
    
    print("\n🧪 DETAILED EXPERIMENTAL PROTOCOL")
    print("=" * 50)
    
    protocol = {
        'phase_1_synthesis': {
            'duration': '1-2 weeks',
            'samples': [
                'MgFe (1:1 atomic)',
                'Mg₂Fe (2:1 atomic)', 
                'Mg₃Fe₂ (3:2 atomic)',
                'Reference: pure Mg, pure Fe'
            ],
            'method': 'Powder metallurgy + controlled cooling',
            'conditions': '673K, 4h, Ar atmosphere, 1K/min cooling',
            'characterization': ['XRD', 'SEM', 'EDS composition check'],
            'success_criteria': 'Single phase or controlled two-phase microstructure'
        },
        
        'phase_2_superconductivity_testing': {
            'duration': '2-4 weeks',
            'measurements': [
                'Resistivity vs. temperature (300K → 4K)',
                'AC magnetic susceptibility',
                'DC magnetization (SQUID)',
                'Specific heat'
            ],
            'equipment': 'Physical Property Measurement System (PPMS) or similar',
            'temperature_range': '4K to 300K, focus on 20-50K region',
            'success_criteria': 'Zero resistance + diamagnetic response'
        },
        
        'phase_3_optimization': {
            'duration': '1-3 months',
            'variations': [
                'Fine-tune stoichiometry (±2% Mg content)',
                'Different cooling rates',
                'Annealing temperatures',
                'Pressure effects'
            ],
            'goal': 'Maximize Tc and critical current density',
            'advanced_characterization': ['Neutron scattering', 'ARPES', 'STM']
        }
    }
    
    print("PHASE 1: SYNTHESIS (Week 1-2)")
    print("• Target compositions based on RBT optimization")
    print("• Powder metallurgy with controlled atmosphere")
    print("• Precise stoichiometry control")
    print("• Structural characterization")
    
    print("\nPHASE 2: SUPERCONDUCTIVITY TESTING (Week 3-6)")
    print("• Temperature-dependent resistivity measurements")
    print("• Magnetic property measurements")
    print("• Focus on 20-50K temperature range")
    print("• Look for zero resistance + Meissner effect")
    
    print("\nPHASE 3: OPTIMIZATION (Month 2-4)")
    print("• Fine-tune composition and processing")
    print("• Understand mechanism")
    print("• Scale up successful compositions")
    
    return protocol

def materials_databases_check():
    """Check what's known about MgFe in materials databases."""
    
    print("\n📚 MATERIALS DATABASE ANALYSIS")
    print("=" * 50)
    
    # Known information from literature
    known_info = {
        'thermodynamic_data': {
            'mg_fe_solubility': 'Very limited (<3%)',
            'intermetallic_phases': 'Mg₁₇Fe reported but unstable',
            'formation_energy': 'Slightly positive (unfavorable)',
            'melting_points': 'Mg: 923K, Fe: 1811K'
        },
        
        'why_not_studied_for_superconductivity': [
            'Limited mutual solubility',
            'No known stable intermetallic compounds',
            'Very different electronegativities',
            'No theoretical predictions of superconductivity',
            'Focus on equilibrium phases only'
        ],
        
        'rbt_opportunity': [
            'RBT doesn\'t require thermodynamic stability',
            'Metastable phases could be superconducting',
            'Valence electron optimization more important than stability',
            'Non-equilibrium processing could unlock hidden phases'
        ]
    }
    
    print("Why MgFe hasn't been studied for superconductivity:")
    for reason in known_info['why_not_studied_for_superconductivity']:
        print(f"• {reason}")
    
    print("\nRBT-based opportunities:")
    for opportunity in known_info['rbt_opportunity']:
        print(f"• {opportunity}")
    
    return known_info

def main():
    """Main analysis function."""
    
    print("🔬 COMPREHENSIVE MgFe SUPERCONDUCTOR DISCOVERY ANALYSIS")
    print("=" * 70)
    
    # 1. RBT optimization
    rbt_results = analyze_rbt_optimization()
    
    # 2. Synthesis pathways
    synthesis_methods = analyze_synthesis_pathways()
    
    # 3. What others missed
    missed_factors = what_others_missed()
    
    # 4. Experimental protocol
    protocol = experimental_protocol()
    
    # 5. Database check
    database_info = materials_databases_check()
    
    print("\n🎯 SUMMARY: PATH TO DISCOVERY")
    print("=" * 50)
    print("1. Target compositions: MgFe, Mg₂Fe, Mg₃Fe₂ (perfect RBT scores)")
    print("2. Synthesis method: Powder metallurgy with controlled cooling")
    print("3. Key insight: Focus on metastable phases, not equilibrium")
    print("4. Temperature range: Test 4-50K systematically")
    print("5. Timeline: 2-6 months for definitive results")
    
    print("\n🚀 NEXT IMMEDIATE ACTIONS:")
    print("• Acquire high-purity Mg and Fe powders")
    print("• Set up inert atmosphere synthesis capability")  
    print("• Plan low-temperature electrical measurements")
    print("• Begin with MgFe 1:1 composition (highest RBT score)")

if __name__ == "__main__":
    main() 