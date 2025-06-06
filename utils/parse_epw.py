#!/usr/bin/env python3
"""
EPW RESULTS PARSER FOR AMBIENT-PRESSURE SUPERCONDUCTORS

Collects and analyzes results from QE+EPW calculations to identify
the most promising room-temperature superconductor candidates.

Author: RT-Superconductor Research Team
Date: 2025-06-07
"""

import os
import json
import glob
import pandas as pd
import numpy as np
from pathlib import Path
import re

def parse_epw_output(epw_file):
    """Parse EPW output file to extract key parameters"""
    
    results = {
        'lambda': None,
        'omega_log': None,
        'tc_allen_dynes': None,
        'status': 'UNKNOWN'
    }
    
    if not os.path.exists(epw_file):
        results['status'] = 'NOT_FOUND'
        return results
    
    try:
        with open(epw_file, 'r') as f:
            content = f.read()
        
        # Extract λ (electron-phonon coupling)
        lambda_match = re.search(r'lambda\s*=\s*([\d\.-]+)', content)
        if lambda_match:
            results['lambda'] = float(lambda_match.group(1))
        
        # Extract ω_log (logarithmic average phonon frequency)
        omega_match = re.search(r'omega_log\s*=\s*([\d\.-]+)', content)
        if omega_match:
            results['omega_log'] = float(omega_match.group(1))
        
        # Calculate Allen-Dynes Tc
        if results['lambda'] is not None and results['omega_log'] is not None:
            lambda_val = results['lambda']
            omega_log = results['omega_log']
            mu_star = 0.10  # Coulomb pseudopotential
            
            if lambda_val > mu_star:
                # Allen-Dynes formula with strong coupling corrections
                f1 = (1 + (lambda_val/2.46) * (1 + 3.8*mu_star)**0.25)**0.25
                f2 = 1 + (lambda_val**2 / (lambda_val**2 + 2.5**2))**1.5 * \
                     (omega_log/250 - 1)
                
                tc = (f1 * f2 * omega_log / 1.2) * \
                     np.exp(-1.04 * (1 + lambda_val) / 
                           (lambda_val - mu_star * (1 + 0.62 * lambda_val)))
                
                results['tc_allen_dynes'] = tc
            else:
                results['tc_allen_dynes'] = 0.0
        
        # Check completion status
        if 'Ended the EPW calculation' in content:
            results['status'] = 'COMPLETED'
        elif 'ERROR' in content.upper():
            results['status'] = 'ERROR'
        else:
            results['status'] = 'RUNNING'
            
    except Exception as e:
        results['status'] = f'PARSE_ERROR: {str(e)}'
    
    return results

def parse_calculation_status(calc_dir):
    """Parse calculation status and metadata"""
    
    status_file = Path(calc_dir) / 'calculation.status'
    log_file = Path(calc_dir) / 'calculation.log'
    summary_file = Path(calc_dir) / 'results_summary.json'
    
    metadata = {
        'directory': str(calc_dir),
        'structure': Path(calc_dir).name,
        'overall_status': 'UNKNOWN',
        'gate1_status': 'UNKNOWN',
        'gate2_status': 'UNKNOWN',
        'calculation_time': None,
        'last_update': None
    }
    
    # Read status file
    if status_file.exists():
        with open(status_file, 'r') as f:
            metadata['overall_status'] = f.read().strip()
    
    # Read summary JSON if available
    if summary_file.exists():
        try:
            with open(summary_file, 'r') as f:
                summary = json.load(f)
                metadata.update(summary)
        except:
            pass
    
    # Read log file for timing
    if log_file.exists():
        try:
            with open(log_file, 'r') as f:
                lines = f.readlines()
                if lines:
                    metadata['last_update'] = lines[-1].split(':')[0]
        except:
            pass
    
    return metadata

def collect_results(qe_runs_dir):
    """Collect results from all QE calculation directories"""
    
    print("📊 COLLECTING EPW RESULTS FROM QE CALCULATIONS")
    print("=" * 60)
    
    qe_runs_path = Path(qe_runs_dir)
    if not qe_runs_path.exists():
        print(f"❌ QE runs directory not found: {qe_runs_dir}")
        return None
    
    # Find all calculation directories
    calc_dirs = [d for d in qe_runs_path.iterdir() if d.is_dir()]
    print(f"🔍 Found {len(calc_dirs)} calculation directories")
    
    results = []
    
    for calc_dir in sorted(calc_dirs):
        print(f"  📁 Processing {calc_dir.name}...")
        
        # Parse calculation metadata
        metadata = parse_calculation_status(calc_dir)
        
        # Parse EPW results
        epw_file = calc_dir / 'epw.out'
        epw_results = parse_epw_output(epw_file)
        
        # Combine all data
        combined_results = {**metadata, **epw_results}
        
        # Determine Gate-1 status (check for gamma phonon completion)
        gamma_file = calc_dir / 'gamma.out'
        if gamma_file.exists():
            try:
                with open(gamma_file, 'r') as f:
                    gamma_content = f.read()
                if 'negative' in gamma_content.lower():
                    combined_results['gate1_status'] = 'FAIL'
                elif 'End of phonon calculation' in gamma_content:
                    combined_results['gate1_status'] = 'PASS'
            except:
                pass
        
        # Determine Gate-2 status
        if (combined_results['lambda'] is not None and 
            combined_results['omega_log'] is not None):
            if (combined_results['lambda'] >= 2.0 and 
                combined_results['omega_log'] >= 800):
                combined_results['gate2_status'] = 'PASS'
            else:
                combined_results['gate2_status'] = 'FAIL'
        
        results.append(combined_results)
        
        # Print status
        status = combined_results['overall_status']
        lambda_val = combined_results['lambda']
        tc = combined_results['tc_allen_dynes']
        
        if status == 'COMPLETED' and lambda_val is not None:
            if combined_results['gate2_status'] == 'PASS':
                print(f"    🔥 λ={lambda_val:.2f}, Tc={tc:.1f}K - GATE-2 PASS!")
            else:
                print(f"    ✅ λ={lambda_val:.2f}, Tc={tc:.1f}K")
        else:
            print(f"    ⏳ Status: {status}")
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Save results
    output_file = Path(qe_runs_dir).parent / 'results' / 'ambient_leaderboard.csv'
    output_file.parent.mkdir(exist_ok=True)
    
    df.to_csv(output_file, index=False)
    
    print("")
    print("📊 RESULTS SUMMARY")
    print("=" * 30)
    print(f"📁 Total calculations: {len(results)}")
    print(f"✅ Completed: {len(df[df['overall_status'] == 'COMPLETED'])}")
    print(f"⏳ Running: {len(df[df['overall_status'].isin(['STARTED', 'RUNNING'])])}")
    print(f"❌ Failed: {len(df[df['overall_status'].str.contains('FAILED')])}")
    
    # Count Gate passes
    gate1_pass = len(df[df['gate1_status'] == 'PASS'])
    gate2_pass = len(df[df['gate2_status'] == 'PASS'])
    
    print(f"🚪 Gate-1 pass: {gate1_pass}")
    print(f"🚪 Gate-2 pass: {gate2_pass}")
    
    if gate2_pass > 0:
        print("")
        print("🔥 GATE-2 WINNERS (Potential RT-Superconductors):")
        winners = df[df['gate2_status'] == 'PASS'].sort_values('tc_allen_dynes', ascending=False)
        for _, row in winners.head(10).iterrows():
            print(f"  🏆 {row['structure']:15} | λ={row['lambda']:.2f} | Tc={row['tc_allen_dynes']:.1f}K")
    
    print(f"\n💾 Results saved to: {output_file}")
    
    return df

def generate_leaderboard_report(df, output_dir):
    """Generate comprehensive leaderboard report"""
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Filter completed calculations with valid results
    completed = df[
        (df['overall_status'] == 'COMPLETED') & 
        (df['lambda'].notna()) & 
        (df['tc_allen_dynes'].notna())
    ].copy()
    
    if len(completed) == 0:
        print("⚠️  No completed calculations with valid results found")
        return
    
    # Sort by Tc
    completed = completed.sort_values('tc_allen_dynes', ascending=False)
    
    # Generate detailed report
    report_file = output_path / 'leaderboard_report.md'
    
    with open(report_file, 'w') as f:
        f.write("# 🏆 Ambient-Pressure Superconductor Leaderboard\n\n")
        f.write("**Generated**: {}\n\n".format(pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')))
        
        # Summary statistics
        f.write("## 📊 Summary Statistics\n\n")
        f.write(f"- **Total candidates tested**: {len(df)}\n")
        f.write(f"- **Completed calculations**: {len(completed)}\n")
        f.write(f"- **Gate-1 passes**: {len(df[df['gate1_status'] == 'PASS'])}\n")
        f.write(f"- **Gate-2 passes**: {len(df[df['gate2_status'] == 'PASS'])}\n")
        f.write(f"- **Highest Tc**: {completed['tc_allen_dynes'].max():.1f} K\n")
        f.write(f"- **Average λ**: {completed['lambda'].mean():.2f}\n\n")
        
        # Top candidates
        f.write("## 🔥 Top Candidates\n\n")
        f.write("| Rank | Structure | λ | ω_log (K) | Tc (K) | Gate-1 | Gate-2 |\n")
        f.write("|------|-----------|---|-----------|--------|--------|--------|\n")
        
        for i, (_, row) in enumerate(completed.head(20).iterrows()):
            g1_icon = "✅" if row['gate1_status'] == 'PASS' else "❌"
            g2_icon = "🔥" if row['gate2_status'] == 'PASS' else "❌"
            
            f.write(f"| {i+1:2d} | {row['structure']:12} | "
                   f"{row['lambda']:.2f} | {row['omega_log']:.0f} | "
                   f"{row['tc_allen_dynes']:.1f} | {g1_icon} | {g2_icon} |\n")
        
        # Gate-2 winners detailed analysis
        winners = completed[completed['gate2_status'] == 'PASS']
        if len(winners) > 0:
            f.write("\n## 🎯 Room-Temperature Candidates (Gate-2 Winners)\n\n")
            for _, row in winners.iterrows():
                f.write(f"### {row['structure']}\n")
                f.write(f"- **λ**: {row['lambda']:.3f}\n")
                f.write(f"- **ω_log**: {row['omega_log']:.1f} K\n")
                f.write(f"- **Tc (Allen-Dynes)**: {row['tc_allen_dynes']:.1f} K\n")
                f.write(f"- **Status**: {row['overall_status']}\n\n")
    
    print(f"📄 Leaderboard report saved to: {report_file}")

def main():
    """Main function for standalone execution"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Parse EPW results and generate leaderboard')
    parser.add_argument('--qe-runs-dir', type=str, default='calc/qe_runs',
                       help='Directory containing QE calculation results')
    parser.add_argument('--output-dir', type=str, default='results',
                       help='Output directory for reports')
    
    args = parser.parse_args()
    
    # Collect results
    df = collect_results(args.qe_runs_dir)
    
    if df is not None:
        # Generate detailed report
        generate_leaderboard_report(df, args.output_dir)
        
        print("\n🚀 Results parsing complete!")
        print("   Check results/ambient_leaderboard.csv for full data")
        print("   Check results/leaderboard_report.md for summary")

if __name__ == "__main__":
    main() 