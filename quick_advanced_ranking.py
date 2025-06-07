#!/usr/bin/env python3
"""
Quick Advanced Ranking: Chemical Pressure + Electronegativity Mismatch
"""

import json

# Pauling electronegativity values
ELECTRONEGATIVITY = {'H': 2.20, 'Be': 1.57, 'B': 2.04, 'C': 2.55}

def calculate_electronegativity_mismatch(formula):
    """Calculate max electronegativity difference in formula"""
    elements = []
    for element in ['B', 'Be', 'C', 'H']:
        if element in formula:
            elements.append(ELECTRONEGATIVITY[element])
    
    if len(elements) >= 2:
        return max(elements) - min(elements)
    return 0.0

def main():
    print("🧠 ADVANCED RANKING: Chemical Pressure + Electronegativity Mismatch")
    print("=" * 65)
    
    # Load survivors
    with open('survivors.txt') as f:
        survivors = [line.strip() for line in f if line.strip()]
    
    # Load chemical pressure data
    with open('data/generated/chemical_pressure_analysis.json') as f:
        data = json.load(f)
    pressure_data = {c['formula']: c['chemical_pressure'] for c in data['chemical_pressure_analysis']}
    
    print(f"📊 Analyzing {len(survivors)} survivors with enhanced methodology...")
    print("")
    print("Formula      | P(GPa) | Δχ    | Score | Category")
    print("-------------|--------|-------|-------|----------")
    
    results = []
    
    for formula in survivors:
        P_internal = pressure_data.get(formula, 0.0)
        dchi = calculate_electronegativity_mismatch(formula)
        
        # Enhanced scoring: P/200 + 0.5*(Δχ/1.0)
        score = (P_internal / 200) + 0.5 * (dchi / 1.0)
        
        results.append({
            'formula': formula,
            'P_internal': P_internal,
            'dchi': dchi,
            'score': score
        })
    
    # Sort by score
    results.sort(key=lambda x: x['score'], reverse=True)
    
    # Display top 15
    breakthrough_candidates = []
    
    for i, r in enumerate(results[:15], 1):
        # Enhanced categorization
        if r['score'] >= 1.6 and r['P_internal'] >= 215:
            category = "🔥 BREAKTHROUGH"
            breakthrough_candidates.append(r['formula'])
        elif r['score'] >= 1.3 and r['P_internal'] >= 200:
            category = "⭐ HIGH"
        elif r['score'] >= 1.0:
            category = "💡 MEDIUM"
        else:
            category = "🔍 LOW"
        
        print(f"{r['formula']:12s} | {r['P_internal']:6.1f} | {r['dchi']:5.2f} | {r['score']:5.2f} | {category}")
    
    print("")
    print(f"🔥 BREAKTHROUGH CANDIDATES ({len(breakthrough_candidates)}):")
    print("=" * 35)
    for formula in breakthrough_candidates:
        result = next(r for r in results if r['formula'] == formula)
        print(f"   🔥 {formula:12s} (Score: {result['score']:.2f}, P: {result['P_internal']:.1f} GPa, Δχ: {result['dchi']:.2f})")
    
    print("")
    print("📊 COMPARISON WITH PREVIOUS RANKING:")
    print("=" * 35)
    
    # Previous top 3 (chemical pressure only)
    previous_top3 = ['BBeC3H5', 'BBeC4H6', 'BeC3H4']
    current_top3 = [r['formula'] for r in results[:3]]
    
    print("Previous Top 3 (Chemical Pressure Only):")
    for i, formula in enumerate(previous_top3, 1):
        result = next((r for r in results if r['formula'] == formula), None)
        if result:
            new_rank = [r['formula'] for r in results].index(formula) + 1
            print(f"   {i}. {formula:12s} → Rank {new_rank:2d} (Score: {result['score']:.2f})")
    
    print("")
    print("NEW Top 3 (Enhanced Ranking):")
    for i, formula in enumerate(current_top3, 1):
        result = next(r for r in results if r['formula'] == formula)
        print(f"   {i}. {formula:12s} (Score: {result['score']:.2f}, P: {result['P_internal']:.1f} GPa, Δχ: {result['dchi']:.2f})")
    
    # Save new priority file
    with open('priority_advanced.txt', 'w') as f:
        for r in results[:10]:  # Top 10
            f.write(f"{r['formula']}\n")
    
    print(f"\n✅ Enhanced top 10 saved to priority_advanced.txt")
    print("🚀 Deploy these optimized candidates for maximum breakthrough probability!")

if __name__ == "__main__":
    main() 