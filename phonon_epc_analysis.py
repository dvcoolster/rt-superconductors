#!/usr/bin/env python
"""
Phonon and Electron-Phonon Coupling (EPC) Analysis for MgFe Superconductors

This script calculates:
- λ (electron-phonon coupling constant)
- ωlog (logarithmic average phonon frequency)
- Conventional Tc using McMillan equation
- Allen-Dynes correction
- Phonon DOS and coupling functions

Authors: Dharamveer Chouhan et al.
Purpose: Provide EPC parameters for experimental validation
"""

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Physical constants
KB = 8.617e-5  # eV/K (Boltzmann constant)
HBAR = 6.582e-16  # eV⋅s (reduced Planck constant)

class PhononEPCAnalyzer:
    """Comprehensive phonon and EPC analysis for MgFe alloys."""
    
    def __init__(self):
        self.structures = {
            'MgFe_cubic': {'mg_frac': 0.5, 'fe_frac': 0.5},
            'Mg2Fe_layered': {'mg_frac': 0.67, 'fe_frac': 0.33},
            'MgFe2_ordered': {'mg_frac': 0.33, 'fe_frac': 0.67},
            'Mg3Fe_intermetallic': {'mg_frac': 0.75, 'fe_frac': 0.25}
        }
        self.results = {}
        
    def calculate_phonon_frequencies(self, structure_name, composition):
        """Calculate phonon frequency spectrum for each structure."""
        
        print(f"\n🌊 PHONON FREQUENCY ANALYSIS: {structure_name}")
        print("=" * 50)
        
        # Element-specific Debye frequencies (experimental values)
        debye_freqs = {
            'Mg': 380,  # cm⁻¹
            'Fe': 470   # cm⁻¹
        }
        
        # Average Debye frequency (composition weighted)
        mg_frac = composition['mg_frac']
        fe_frac = composition['fe_frac']
        
        avg_debye = mg_frac * debye_freqs['Mg'] + fe_frac * debye_freqs['Fe']
        
        # Generate phonon spectrum (simplified Debye model)
        freq_max = avg_debye * 1.2  # Extend beyond Debye
        frequencies = np.linspace(1, freq_max, 1000)  # cm⁻¹
        
        # Debye density of states: g(ω) ∝ ω² for ω < ωD
        phonon_dos = np.zeros_like(frequencies)
        debye_cutoff = avg_debye
        
        for i, freq in enumerate(frequencies):
            if freq <= debye_cutoff:
                phonon_dos[i] = 3 * freq**2 / debye_cutoff**3
            else:
                # Exponential decay beyond Debye cutoff
                phonon_dos[i] = 3 * np.exp(-(freq - debye_cutoff) / (0.1 * debye_cutoff))
        
        # Normalize DOS
        phonon_dos = phonon_dos / np.trapz(phonon_dos, frequencies)
        
        # Calculate key phonon parameters
        
        # 1. Average phonon frequency <ω>
        omega_avg = np.trapz(frequencies * phonon_dos, frequencies)
        
        # 2. Logarithmic average frequency ωlog
        # ωlog = exp(∫ ln(ω) g(ω) dω / ∫ g(ω) dω)
        log_omega_weighted = np.trapz(np.log(frequencies) * phonon_dos, frequencies)
        omega_log = np.exp(log_omega_weighted)
        
        # 3. Second moment <ω²>
        omega_sq_avg = np.trapz(frequencies**2 * phonon_dos, frequencies)
        
        print(f"  Average Debye frequency: {avg_debye:.1f} cm⁻¹")
        print(f"  <ω> (first moment): {omega_avg:.1f} cm⁻¹")
        print(f"  ωlog (logarithmic avg): {omega_log:.1f} cm⁻¹")
        print(f"  <ω²>^(1/2) (RMS): {np.sqrt(omega_sq_avg):.1f} cm⁻¹")
        
        return {
            'frequencies': frequencies,
            'phonon_dos': phonon_dos,
            'avg_debye': avg_debye,
            'omega_avg': omega_avg,
            'omega_log': omega_log,
            'omega_sq_avg': omega_sq_avg
        }
    
    def calculate_electron_phonon_coupling(self, structure_name, composition, phonon_data):
        """Calculate electron-phonon coupling parameters."""
        
        print(f"\n⚡ ELECTRON-PHONON COUPLING: {structure_name}")
        print("=" * 50)
        
        # Electronic density of states at Fermi level (from band structure)
        # Estimated from composition (states/eV/unit cell)
        mg_dos = 2.5  # states/eV (s,p electrons)
        fe_dos = 8.5  # states/eV (s,p,d electrons, larger d-DOS)
        
        mg_frac = composition['mg_frac']
        fe_frac = composition['fe_frac']
        
        n_ef = mg_frac * mg_dos + fe_frac * fe_dos  # Total DOS at E_F
        
        # Electron-phonon matrix elements (composition dependent)
        # |M|² depends on orbital overlap and atomic mass
        mg_mass = 24.305  # amu
        fe_mass = 55.845  # amu
        avg_mass = mg_frac * mg_mass + fe_frac * fe_mass
        
        # Matrix element squared (empirical scaling)
        # Stronger for d-electrons (Fe) due to localization
        mg_matrix_sq = 0.15  # eV²/amu (s,p coupling)
        fe_matrix_sq = 0.35  # eV²/amu (d-electron coupling stronger)
        
        avg_matrix_sq = mg_frac * mg_matrix_sq + fe_frac * fe_matrix_sq
        
        # Electron-phonon coupling function α²F(ω)
        frequencies = phonon_data['frequencies']
        phonon_dos = phonon_data['phonon_dos']
        
        # α²F(ω) = N(E_F) |M|² g(ω) / M̄ω
        # Convert frequencies to eV (1 cm⁻¹ = 1.24e-4 eV)
        freq_ev = frequencies * 1.24e-4
        
        alpha2f = np.zeros_like(frequencies)
        for i, (freq_cm, freq_eV, gw) in enumerate(zip(frequencies, freq_ev, phonon_dos)):
            if freq_eV > 0:
                alpha2f[i] = n_ef * avg_matrix_sq * gw / (avg_mass * freq_eV)
        
        # Calculate λ (total electron-phonon coupling constant)
        # λ = 2 ∫ α²F(ω)/ω dω
        lambda_ep = 2 * np.trapz(alpha2f / freq_ev, freq_ev)
        
        # Calculate logarithmic average frequency for superconductivity
        # ωlog = exp((2/λ) ∫ (α²F(ω)/ω) ln(ω) dω)
        if lambda_ep > 0:
            log_omega_weighted = np.trapz((alpha2f / freq_ev) * np.log(freq_ev), freq_ev)
            omega_log_sc = np.exp((2 / lambda_ep) * log_omega_weighted)
            omega_log_sc_cm = omega_log_sc / 1.24e-4  # Convert back to cm⁻¹
        else:
            omega_log_sc_cm = 0
        
        # Effective Coulomb repulsion μ* (typical values)
        mu_star = 0.13  # Typical for metals
        
        print(f"  N(E_F): {n_ef:.2f} states/eV/unit cell")
        print(f"  <|M|²>: {avg_matrix_sq:.3f} eV²/amu")
        print(f"  λ (EPC constant): {lambda_ep:.3f}")
        print(f"  ωlog (SC): {omega_log_sc_cm:.1f} cm⁻¹")
        print(f"  μ*: {mu_star:.3f}")
        
        return {
            'alpha2f': alpha2f,
            'lambda_ep': lambda_ep,
            'omega_log_sc': omega_log_sc_cm,
            'n_ef': n_ef,
            'mu_star': mu_star,
            'avg_matrix_sq': avg_matrix_sq
        }
    
    def calculate_superconducting_tc(self, structure_name, phonon_data, epc_data):
        """Calculate superconducting Tc using various formulas."""
        
        print(f"\n🌡️ SUPERCONDUCTING TC CALCULATIONS: {structure_name}")
        print("=" * 50)
        
        lambda_ep = epc_data['lambda_ep']
        omega_log = epc_data['omega_log_sc']  # cm⁻¹
        mu_star = epc_data['mu_star']
        
        # Convert ωlog to K
        omega_log_k = omega_log * 1.44  # cm⁻¹ to K conversion
        
        results = {}
        
        # 1. McMillan equation (1968)
        if lambda_ep > 0 and lambda_ep < 1.5:
            tc_mcmillan = (omega_log_k / 1.2) * np.exp(-1.04 * (1 + lambda_ep) / (lambda_ep - mu_star * (1 + 0.62 * lambda_ep)))
            tc_mcmillan = max(tc_mcmillan, 0)
        else:
            tc_mcmillan = 0
        
        # 2. Allen-Dynes equation (1975) - more accurate for strong coupling
        if lambda_ep > 0:
            # Allen-Dynes correction factor
            f1 = (1 + (lambda_ep / 2.46) * (1 + 3.8 * mu_star))**(1/3)
            f2 = 1 + ((lambda_ep - mu_star * (1 + lambda_ep)) / (lambda_ep + 0.1))**2 / ((lambda_ep - mu_star * (1 + lambda_ep))**2 + 0.25)
            
            tc_allen_dynes = (f1 * f2 * omega_log_k / 1.2) * np.exp(-1.04 * (1 + lambda_ep) / (lambda_ep - mu_star * (1 + 0.62 * lambda_ep)))
            tc_allen_dynes = max(tc_allen_dynes, 0)
        else:
            tc_allen_dynes = 0
        
        # 3. Strong coupling limit (Carbotte, 1990)
        if lambda_ep > 1.5:
            tc_strong = omega_log_k * np.exp(-1 / lambda_ep) / 3
        else:
            tc_strong = 0
        
        # 4. Simplified BCS estimate
        tc_bcs_simple = 1.14 * omega_log_k * np.exp(-1 / lambda_ep) if lambda_ep > 0.3 else 0
        
        results = {
            'tc_mcmillan': tc_mcmillan,
            'tc_allen_dynes': tc_allen_dynes,
            'tc_strong_coupling': tc_strong,
            'tc_bcs_simple': tc_bcs_simple,
            'lambda_ep': lambda_ep,
            'omega_log_k': omega_log_k,
            'coupling_regime': self._classify_coupling(lambda_ep)
        }
        
        print(f"  McMillan Tc: {tc_mcmillan:.2f} K")
        print(f"  Allen-Dynes Tc: {tc_allen_dynes:.2f} K")
        print(f"  Strong coupling Tc: {tc_strong:.2f} K")
        print(f"  BCS simple Tc: {tc_bcs_simple:.2f} K")
        print(f"  Coupling regime: {results['coupling_regime']}")
        
        return results
    
    def _classify_coupling(self, lambda_ep):
        """Classify electron-phonon coupling strength."""
        if lambda_ep < 0.3:
            return "Weak coupling"
        elif lambda_ep < 0.7:
            return "Intermediate coupling"
        elif lambda_ep < 1.5:
            return "Strong coupling"
        else:
            return "Very strong coupling"
    
    def run_complete_analysis(self):
        """Run complete phonon/EPC analysis for all structures."""
        
        print("🔬 COMPREHENSIVE PHONON/EPC ANALYSIS FOR MGFE ALLOYS")
        print("=" * 65)
        print("Calculating λ, ωlog, and conventional Tc predictions")
        print("=" * 65)
        
        for structure_name, composition in self.structures.items():
            print(f"\n📊 ANALYZING: {structure_name}")
            print("=" * 60)
            
            # Calculate phonon properties
            phonon_data = self.calculate_phonon_frequencies(structure_name, composition)
            
            # Calculate electron-phonon coupling
            epc_data = self.calculate_electron_phonon_coupling(structure_name, composition, phonon_data)
            
            # Calculate superconducting Tc
            tc_data = self.calculate_superconducting_tc(structure_name, phonon_data, epc_data)
            
            # Store complete results
            self.results[structure_name] = {
                'composition': composition,
                'phonon': phonon_data,
                'epc': epc_data,
                'tc': tc_data
            }
        
        return self.results
    
    def save_results(self, output_dir="results/phonon_epc"):
        """Save all phonon/EPC results to files."""
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save complete results as JSON
        results_clean = {}
        for structure, data in self.results.items():
            results_clean[structure] = {
                'composition': data['composition'],
                'phonon': {
                    'avg_debye': data['phonon']['avg_debye'],
                    'omega_avg': data['phonon']['omega_avg'],
                    'omega_log': data['phonon']['omega_log'],
                    'omega_sq_avg': data['phonon']['omega_sq_avg']
                },
                'epc': {
                    'lambda_ep': data['epc']['lambda_ep'],
                    'omega_log_sc': data['epc']['omega_log_sc'],
                    'n_ef': data['epc']['n_ef'],
                    'mu_star': data['epc']['mu_star']
                },
                'tc': data['tc']
            }
        
        with open(output_path / "phonon_epc_results.json", 'w') as f:
            json.dump(results_clean, f, indent=2, default=str)
        
        # Create summary CSV for experimental reference
        import csv
        
        with open(output_path / "epc_parameters_summary.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Structure', 'Mg_fraction', 'Fe_fraction', 'lambda_ep', 
                'omega_log_cm', 'omega_log_K', 'Tc_McMillan_K', 
                'Tc_AllenDynes_K', 'Coupling_regime', 'N_EF'
            ])
            
            for structure, data in self.results.items():
                comp = data['composition']
                epc = data['epc']
                tc = data['tc']
                
                writer.writerow([
                    structure,
                    f"{comp['mg_frac']:.3f}",
                    f"{comp['fe_frac']:.3f}",
                    f"{epc['lambda_ep']:.3f}",
                    f"{epc['omega_log_sc']:.1f}",
                    f"{tc['omega_log_k']:.1f}",
                    f"{tc['tc_mcmillan']:.2f}",
                    f"{tc['tc_allen_dynes']:.2f}",
                    tc['coupling_regime'],
                    f"{epc['n_ef']:.2f}"
                ])
        
        print(f"\n💾 RESULTS SAVED TO: {output_path}")
        print("📄 Files created:")
        print("  • phonon_epc_results.json - Complete analysis data")
        print("  • epc_parameters_summary.csv - Key parameters for experiments")
        
        return output_path
    
    def create_plots(self, output_dir="results/phonon_epc"):
        """Create visualization plots for phonon/EPC analysis."""
        
        output_path = Path(output_dir)
        
        # Plot 1: Lambda vs Tc comparison
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        structures = list(self.results.keys())
        lambdas = [self.results[s]['epc']['lambda_ep'] for s in structures]
        tc_mcmillan = [self.results[s]['tc']['tc_mcmillan'] for s in structures]
        tc_allen_dynes = [self.results[s]['tc']['tc_allen_dynes'] for s in structures]
        
        # Lambda comparison
        x_pos = np.arange(len(structures))
        ax1.bar(x_pos, lambdas, color='skyblue', alpha=0.7)
        ax1.set_xlabel('Structure')
        ax1.set_ylabel('λ (EPC constant)')
        ax1.set_title('Electron-Phonon Coupling Strength')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels([s.replace('_', '\n') for s in structures], rotation=45)
        ax1.grid(True, alpha=0.3)
        
        # Tc comparison
        width = 0.35
        ax2.bar(x_pos - width/2, tc_mcmillan, width, label='McMillan', alpha=0.7)
        ax2.bar(x_pos + width/2, tc_allen_dynes, width, label='Allen-Dynes', alpha=0.7)
        ax2.set_xlabel('Structure')
        ax2.set_ylabel('Tc (K)')
        ax2.set_title('Superconducting Tc Predictions')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels([s.replace('_', '\n') for s in structures], rotation=45)
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_path / "phonon_epc_analysis.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"📊 Plot saved: {output_path}/phonon_epc_analysis.png")

def main():
    """Main phonon/EPC analysis workflow."""
    
    # Initialize analyzer
    analyzer = PhononEPCAnalyzer()
    
    # Run complete analysis
    results = analyzer.run_complete_analysis()
    
    # Save results
    output_path = analyzer.save_results()
    
    # Create plots
    analyzer.create_plots()
    
    # Print summary
    print("\n🎯 PHONON/EPC ANALYSIS SUMMARY")
    print("=" * 50)
    
    for structure, data in results.items():
        tc = data['tc']
        epc = data['epc']
        print(f"\n{structure}:")
        print(f"  λ = {epc['lambda_ep']:.3f}")
        print(f"  ωlog = {epc['omega_log_sc']:.1f} cm⁻¹")
        print(f"  Tc (McMillan) = {tc['tc_mcmillan']:.2f} K")
        print(f"  Tc (Allen-Dynes) = {tc['tc_allen_dynes']:.2f} K")
        print(f"  Regime: {tc['coupling_regime']}")
    
    print(f"\n📁 All outputs saved in: {output_path}")
    print("🚀 Ready for experimental validation!")

if __name__ == "__main__":
    main() 