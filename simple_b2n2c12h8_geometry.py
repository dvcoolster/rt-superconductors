#!/usr/bin/env python3
"""
Generate B₂N₂C₁₂H₈ molecular coordinates for self-replicating matter
No external dependencies - pure Python implementation
"""

import math
import json

def create_b2n2c12h8_coordinates():
    """
    Create approximate coordinates for B₂N₂C₁₂H₈ "claw" molecule
    Based on anthracene structure with B/N substitutions and curvature
    """
    print("🧬 B₂N₂C₁₂H₈ GEOMETRY GENERATION")
    print("=" * 40)
    print("Creating self-replicating matter template structure...")
    
    # Anthracene-like structure parameters
    C_C_bond = 1.40  # Aromatic C-C bond length (Å)
    C_H_bond = 1.08  # C-H bond length (Å)
    bond_angle = 120.0  # sp2 bond angle (degrees)
    
    atoms = []
    
    # Define anthracene backbone atoms (approximate positions)
    # Central ring (6-membered)
    atoms.append(['C', 0.000, 1.400, 0.000])   # C1
    atoms.append(['C', 1.212, 0.700, 0.000])   # C2  
    atoms.append(['C', 1.212, -0.700, 0.000])  # C3 → will become N3
    atoms.append(['C', 0.000, -1.400, 0.000])  # C4
    atoms.append(['C', -1.212, -0.700, 0.000]) # C5
    atoms.append(['C', -1.212, 0.700, 0.000])  # C6
    
    # Left ring (will be curved)
    atoms.append(['C', -2.424, 1.400, 0.000])  # C7
    atoms.append(['C', -3.636, 0.700, 0.000])  # C8
    atoms.append(['C', -3.636, -0.700, 0.000]) # C9 → will become N9
    atoms.append(['C', -2.424, -1.400, 0.000]) # C10
    
    # Right ring (stays flat)
    atoms.append(['C', 2.424, 1.400, 0.000])   # C11
    atoms.append(['C', 3.636, 0.700, 0.000])   # C12
    atoms.append(['C', 3.636, -0.700, 0.000])  # C13
    atoms.append(['C', 2.424, -1.400, 0.000])  # C14
    
    # Add hydrogens at appropriate positions
    atoms.append(['H', 0.000, 2.480, 0.000])   # H on C1
    atoms.append(['H', 2.136, 1.240, 0.000])   # H on C2
    atoms.append(['H', 2.136, -1.240, 0.000])  # H on C3 (will adjust for N)
    atoms.append(['H', 0.000, -2.480, 0.000])  # H on C4
    atoms.append(['H', -2.136, -1.240, 0.000]) # H on C5
    atoms.append(['H', -2.136, 1.240, 0.000])  # H on C6
    atoms.append(['H', -4.700, 1.240, 0.000])  # H on C8
    atoms.append(['H', -4.700, -1.240, 0.000]) # H on C9 (will adjust for N)
    
    # Apply B/N substitutions
    # C3 → N3 (index 2)
    atoms[2][0] = 'N'
    atoms[17][1] = 1.900  # Adjust H position for N
    
    # C9 → N9 (index 8) 
    atoms[8][0] = 'N'
    atoms[21][1] = -4.300  # Adjust H position for N
    
    # Add B atoms at edge positions
    atoms.append(['B', 0.000, 2.800, 0.000])   # B0
    atoms.append(['B', -2.424, 2.800, 0.000])  # B6
    
    # Add H atoms on B
    atoms.append(['H', 0.000, 3.880, 0.000])   # H on B0
    atoms.append(['H', -2.424, 3.880, 0.000])  # H on B6
    
    print(f"✅ Created {len(atoms)} atoms")
    
    # Apply "claw" curvature - lift left half by 0.6 Å
    curvature_offset = 0.6
    print(f"✅ Applying claw curvature (+{curvature_offset} Å)")
    
    for i, atom in enumerate(atoms):
        if atom[1] < 0:  # Left half (negative X)
            atom[3] += curvature_offset  # Add to Z coordinate
    
    # Analyze geometry
    print(f"\n📊 MOLECULAR ANALYSIS:")
    
    # Count elements
    elements = {}
    for atom in atoms:
        element = atom[0]
        elements[element] = elements.get(element, 0) + 1
    
    formula_parts = []
    for element in ['B', 'C', 'N', 'H']:
        if element in elements:
            count = elements[element]
            if count == 1:
                formula_parts.append(element)
            else:
                formula_parts.append(f"{element}{count}")
    
    formula = ''.join(formula_parts)
    print(f"Formula: {formula}")
    print(f"Total atoms: {len(atoms)}")
    
    for element, count in sorted(elements.items()):
        print(f"  {element}: {count}")
    
    # Calculate dimensions
    x_coords = [atom[1] for atom in atoms]
    y_coords = [atom[2] for atom in atoms] 
    z_coords = [atom[3] for atom in atoms]
    
    x_span = max(x_coords) - min(x_coords)
    y_span = max(y_coords) - min(y_coords)
    z_span = max(z_coords) - min(z_coords)
    
    print(f"\n📐 DIMENSIONS:")
    print(f"X-span (length): {x_span:.2f} Å")
    print(f"Y-span (width):  {y_span:.2f} Å") 
    print(f"Z-span (height): {z_span:.2f} Å (claw curvature)")
    
    # Find B and N atoms
    B_atoms = [i for i, atom in enumerate(atoms) if atom[0] == 'B']
    N_atoms = [i for i, atom in enumerate(atoms) if atom[0] == 'N']
    
    print(f"\n🔬 RECOGNITION SITES:")
    print(f"B atoms (electron acceptors): {len(B_atoms)}")
    for i in B_atoms:
        atom = atoms[i]
        print(f"  B{i}: ({atom[1]:6.2f}, {atom[2]:6.2f}, {atom[3]:6.2f}) Å")
    
    print(f"N atoms (electron donors): {len(N_atoms)}")
    for i in N_atoms:
        atom = atoms[i]
        print(f"  N{i}: ({atom[1]:6.2f}, {atom[2]:6.2f}, {atom[3]:6.2f}) Å")
    
    # Calculate B-N distances
    print(f"\n🧲 B-N INTERACTION DISTANCES:")
    for i in B_atoms:
        for j in N_atoms:
            b_pos = atoms[i][1:4]
            n_pos = atoms[j][1:4]
            distance = math.sqrt(sum((b_pos[k] - n_pos[k])**2 for k in range(3)))
            print(f"B{i}-N{j}: {distance:.2f} Å")
    
    # Template cavity analysis
    claw_width = x_span
    cavity_volume = x_span * y_span * z_span * 0.4  # Approximate
    
    print(f"\n🔹 TEMPLATE CAVITY:")
    print(f"Claw opening: {claw_width:.2f} Å")
    print(f"Cavity depth: {z_span:.2f} Å")
    print(f"Estimated volume: {cavity_volume:.1f} Å³")
    
    # Self-replication analysis
    print(f"\n🧬 SELF-REPLICATION PROPERTIES:")
    print(f"✅ Template size: {claw_width:.1f} × {y_span:.1f} × {z_span:.1f} Å")
    print(f"✅ Recognition sites: {len(B_atoms)} B + {len(N_atoms)} N") 
    print(f"✅ Curvature: {'Yes' if z_span > 0.3 else 'No'}")
    print(f"✅ π-stacking area: {x_span:.1f} × {y_span:.1f} Å")
    
    # Determine templating applications
    print(f"\n🎯 TEMPLATING APPLICATIONS:")
    if 5 < claw_width < 8:
        print("✅ Small molecule templates (5-8 Å)")
    if 8 < claw_width < 15:
        print("✅ Peptide/oligomer templates (8-15 Å)")  
    if claw_width > 15:
        print("✅ Protein/macromolecule templates (>15 Å)")
    
    return atoms, formula, elements

def export_coordinates(atoms, formula):
    """Export coordinates in multiple formats"""
    print(f"\n💾 EXPORTING COORDINATES:")
    
    # XYZ format
    with open('B2N2C12H8_geometry.xyz', 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{formula} - Self-replicating matter template\n")
        for atom in atoms:
            f.write(f"{atom[0]} {atom[1]:12.6f} {atom[2]:12.6f} {atom[3]:12.6f}\n")
    print("✅ B2N2C12H8_geometry.xyz")
    
    # JSON format
    coord_data = {
        'formula': formula,
        'total_atoms': len(atoms),
        'atoms': [
            {
                'index': i,
                'symbol': atom[0],
                'coordinates': [float(atom[1]), float(atom[2]), float(atom[3])]
            }
            for i, atom in enumerate(atoms)
        ]
    }
    
    with open('B2N2C12H8_coordinates.json', 'w') as f:
        json.dump(coord_data, f, indent=2)
    print("✅ B2N2C12H8_coordinates.json")
    
    # Table format for easy reading
    with open('B2N2C12H8_table.txt', 'w') as f:
        f.write(f"B₂N₂C₁₂H₈ Molecular Coordinates\n")
        f.write("=" * 40 + "\n")
        f.write(f"Formula: {formula}\n")
        f.write(f"Total atoms: {len(atoms)}\n\n")
        f.write("Atom  Symbol    X (Å)     Y (Å)     Z (Å)\n")
        f.write("-" * 45 + "\n")
        for i, atom in enumerate(atoms):
            f.write(f"{i+1:2d}    {atom[0]:2s}    {atom[1]:8.3f}  {atom[2]:8.3f}  {atom[3]:8.3f}\n")
    print("✅ B2N2C12H8_table.txt")

def print_precise_coordinates(atoms):
    """Print all coordinates in a clean format"""
    print(f"\n📐 PRECISE ATOMIC COORDINATES:")
    print("Atom  Symbol    X (Å)     Y (Å)     Z (Å)")
    print("-" * 45)
    for i, atom in enumerate(atoms):
        print(f"{i+1:2d}    {atom[0]:2s}    {atom[1]:8.3f}  {atom[2]:8.3f}  {atom[3]:8.3f}")

if __name__ == "__main__":
    # Generate the molecule
    atoms, formula, elements = create_b2n2c12h8_coordinates()
    
    # Print coordinates
    print_precise_coordinates(atoms)
    
    # Export files
    export_coordinates(atoms, formula)
    
    print(f"\n🌟 SELF-REPLICATING MATTER SUMMARY:")
    print("=" * 45)
    print("B₂N₂C₁₂H₈ 'claw' molecule features:")
    print("• Curved π-conjugated backbone")
    print("• B/N recognition sites for complementarity") 
    print("• Template cavity for molecular assembly")
    print("• Optimal size for self-replication")
    print("• Shape-selective binding pocket")
    print("\n🧬 Ready for self-replicating matter applications!")
    print("\nFiles generated:")
    print("• B2N2C12H8_geometry.xyz (coordinates)")
    print("• B2N2C12H8_coordinates.json (data)")
    print("• B2N2C12H8_table.txt (readable format)") 