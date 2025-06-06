#!/usr/bin/env python3
"""
Room-Temperature Superconductor DFT Master Script

Coordinates DFPT+EPW calculations and analyzes results for:
- YH9, ScH9, CaH10 (hydride superconductors)  
- Li2BH12, Be4CH10 (complex hydrides)
- MgFe (RBT prediction control)

Real DFT calculations with actual λ, ω_log, and Tc values.
"""

import os
import sys
import json
import subprocess
import numpy as np
import re
import math
from pathlib import Path

def check_qe_installation():
    """Check if Quantum ESPRESSO is properly installed."""
    
    print("🔬 CHECKING QUANTUM ESPRESSO INSTALLATION")
    print("=" * 50)
    
    required_programs = ['pw.x', 'ph.x', 'epw.x', 'q2r.x', 'matdyn.x']
    missing = []
    
    for program in required_programs:
        if subprocess.run(['which', program], capture_output=True).returncode != 0:
            missing.append(program)
    
    if missing:
        print(f"❌ Missing: {', '.join(missing)}")
        print("\n📦 INSTALLATION INSTRUCTIONS:")
        print("1. Download QE 7.2+ from quantum-espresso.org")
        print("2. Configure with: ./configure --enable-epw")
        print("3. Compile with: make all epw")
        print("4. Add to PATH")
        return False
    
    print("✅ All QE programs found")
    return True

def download_pseudopotentials():
    """Download required pseudopotentials."""
    
    pseudo_dir = Path('pseudos')
    pseudo_dir.mkdir(exist_ok=True)
    
    pseudos = [
        'H.pbe-rrkjus_psl.1.0.0.UPF',
        'Y.pbe-spn-rrkjus_psl.1.0.0.UPF', 
        'Sc.pbe-spn-rrkjus_psl.1.0.0.UPF',
        'Ca.pbe-spn-rrkjus_psl.1.0.0.UPF',
        'Li.pbe-s-rrkjus_psl.1.0.0.UPF',
        'B.pbe-n-rrkjus_psl.1.0.0.UPF',
        'Be.pbe-n-rrkjus_psl.1.0.0.UPF',
        'C.pbe-n-rrkjus_psl.1.0.0.UPF',
        'Mg.pbe-nspn-rrkjus_psl.1.0.0.UPF',
        'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'
    ]
    
    print("📁 DOWNLOADING PSEUDOPOTENTIALS")
    print("-" * 40)
    
    base_url = "https://www.quantum-espresso.org/upf_files/"
    
    for pseudo in pseudos:
        pseudo_path = pseudo_dir / pseudo
        if not pseudo_path.exists():
            print(f"Downloading {pseudo}...")
            try:
                subprocess.run(['curl', '-o', str(pseudo_path), f"{base_url}{pseudo}"], 
                              check=True, capture_output=True)
                print(f"✅ {pseudo}")
            except:
                print(f"❌ Failed to download {pseudo}")
                print(f"Manual download from: {base_url}{pseudo}")
        else:
            print(f"✅ {pseudo} (already exists)")

def submit_slurm_job(system, script_type='qe'):
    """Submit SLURM job for a system."""
    
    system_dir = Path(system)
    if not system_dir.exists():
        print(f"❌ Directory not found: {system}")
        return None
    
    os.chdir(system_dir)
    
    script = f"../tools/{script_type}_run.sh"
    result = subprocess.run(['sbatch', script], capture_output=True, text=True)
    
    if result.returncode == 0:
        job_id = result.stdout.strip().split()[-1]
        print(f"✅ {system} {script_type} job: {job_id}")
        return job_id
    else:
        print(f"❌ Failed to submit {system} {script_type}")
        return None

def analyze_epw_results(system):
    """Analyze EPW output and calculate Tc."""
    
    epw_file = Path(f"{system}/epw.out")
    if not epw_file.exists():
        print(f"❌ {system}: EPW output not found")
        return None
    
    with open(epw_file, 'r') as f:
        content = f.read()
    
    # Extract lambda and omega_log
    lambda_match = re.search(r'lambda\s*=\s*([\d.]+)', content)
    omega_match = re.search(r'omega_log\s*=\s*([\d.]+)', content)
    
    if not (lambda_match and omega_match):
        print(f"❌ {system}: Could not extract λ and ω_log")
        return None
    
    λ = float(lambda_match.group(1))
    ω = float(omega_match.group(1))  # in K
    μ_star = 0.10
    
    # Calculate Tc using Allen-Dynes formula
    if λ > μ_star:
        f1 = (1 + (λ/2.46) * (1 + 3.8*μ_star))**(1/3)
        f2 = 1 + ((λ - μ_star*(1+λ))/(λ + 0.15))**2 / 15
        tc = (f1 * f2 * ω / 1.2) * math.exp(-1.04 * (1+λ) / (λ - μ_star*(1+0.62*λ)))
        tc = max(tc, 0)
    else:
        tc = 0
    
    # Coupling regime
    if λ < 0.3:
        regime = "Weak"
    elif λ < 0.7:
        regime = "Intermediate" 
    elif λ < 1.5:
        regime = "Strong"
    else:
        regime = "Very strong"
    
    results = {
        'system': system,
        'lambda': λ,
        'omega_log_K': ω,
        'tc_allen_dynes_K': tc,
        'coupling_regime': regime,
        'rt_potential': tc > 273  # Room temperature
    }
    
    print(f"\n📊 {system.upper()} RESULTS:")
    print(f"  λ = {λ:.3f}")
    print(f"  ω_log = {ω:.1f} K") 
    print(f"  Tc = {tc:.1f} K")
    print(f"  Regime: {regime} coupling")
    print(f"  RT potential: {'🔥 YES' if tc > 273 else '❄️ NO'}")
    
    return results

def main_workflow():
    """Main calculation workflow."""
    
    print("🚀 RT-SUPERCONDUCTOR DFT CALCULATION SUITE")
    print("=" * 60)
    print("Real DFPT+EPW calculations for room-temperature superconductors")
    print("=" * 60)
    
    # Systems to calculate
    systems = ['YH9', 'ScH9', 'CaH10', 'Li2BH12', 'Be4CH10', 'MgFe']
    
    # Check prerequisites
    if not check_qe_installation():
        print("\n❌ Quantum ESPRESSO not found. Please install first.")
        return
    
    # Download pseudopotentials
    download_pseudopotentials()
    
    # Submit calculations
    print(f"\n🔄 SUBMITTING CALCULATIONS FOR {len(systems)} SYSTEMS")
    print("-" * 50)
    
    job_ids = {}
    for system in systems:
        qe_job = submit_slurm_job(system, 'qe')
        if qe_job:
            epw_job = submit_slurm_job(system, 'epw')
            job_ids[system] = {'qe': qe_job, 'epw': epw_job}
    
    print(f"\n✅ Submitted {len(job_ids)} calculation pairs")
    print("💡 Monitor with: squeue -u $USER")
    print("🔍 Analyze when complete: python3 master_calculations.py --analyze")

def analyze_workflow():
    """Analyze completed calculations."""
    
    print("🔍 ANALYZING COMPLETED DFT CALCULATIONS")
    print("=" * 50)
    
    systems = ['YH9', 'ScH9', 'CaH10', 'Li2BH12', 'Be4CH10', 'MgFe']
    all_results = []
    
    for system in systems:
        result = analyze_epw_results(system)
        if result:
            all_results.append(result)
    
    if not all_results:
        print("❌ No results found. Check if calculations completed.")
        return
    
    # Sort by Tc
    all_results.sort(key=lambda x: x['tc_allen_dynes_K'], reverse=True)
    
    print("\n🏆 ROOM-TEMPERATURE SUPERCONDUCTOR RANKINGS")
    print("=" * 60)
    print(f"{'System':<12} {'λ':<8} {'ω_log(K)':<10} {'Tc(K)':<8} {'Regime':<12} {'RT?'}")
    print("-" * 60)
    
    for r in all_results:
        rt_status = "🔥 YES" if r['rt_potential'] else "❄️ NO"
        print(f"{r['system']:<12} {r['lambda']:<8.3f} {r['omega_log_K']:<10.1f} "
              f"{r['tc_allen_dynes_K']:<8.1f} {r['coupling_regime']:<12} {rt_status}")
    
    # Save results
    output_dir = Path('../results/dft_real')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    with open(output_dir / 'rt_superconductor_dft_results.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n💾 Results saved to: {output_dir}")
    
    # Identify room-temperature candidates
    rt_candidates = [r for r in all_results if r['rt_potential']]
    
    if rt_candidates:
        print(f"\n🔥 ROOM-TEMPERATURE SUPERCONDUCTOR CANDIDATES: {len(rt_candidates)}")
        for candidate in rt_candidates:
            print(f"  🌟 {candidate['system']}: Tc = {candidate['tc_allen_dynes_K']:.1f} K")
    else:
        print("\n❄️ No room-temperature candidates found")
        print("🔍 Highest Tc candidate:", all_results[0]['system'], 
              f"({all_results[0]['tc_allen_dynes_K']:.1f} K)")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == '--analyze':
        analyze_workflow()
    else:
        main_workflow() 