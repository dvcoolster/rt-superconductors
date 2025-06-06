#!/usr/bin/env python
"""
Create Excel version of measurement checklist template for partner labs.
"""

import pandas as pd
from pathlib import Path
import numpy as np

def create_excel_template():
    """Create comprehensive Excel measurement template."""
    
    print("📊 Creating Excel measurement template for partner labs...")
    
    # Sample data with examples
    sample_data = [
        {
            'Sample_ID': 'MgFe-001',
            'Composition': 'Mg50Fe50',
            'Synthesis_Method': 'Arc_Melting',
            'Synthesis_Date': '2024-01-15',
            'Operator': 'J.Smith',
            'Purity_Mg_Percent': 99.9,
            'Purity_Fe_Percent': 99.95,
            'Crucible_Material': 'Tungsten',
            'Atmosphere': 'Argon',
            'Max_Temperature_C': 1800,
            'Cooling_Rate_K_min': 1.0,
            'Annealing_Temperature_C': 800,
            'Annealing_Time_hours': 24,
            'Final_Mass_g': 2.45,
            'Color_Appearance': 'Dark_Gray',
            'Crystal_Structure_XRD': 'Cubic',
            'Lattice_Parameters_a_b_c': '4.02_4.02_4.02',
            'Density_g_cm3': 4.2,
            'Resistivity_300K_ohm_cm': 0.025,
            'Resistivity_77K_ohm_cm': 0.018,
            'Resistivity_4K_ohm_cm': 0.0,
            'Tc_onset_K': 42.3,
            'Tc_zero_K': 41.8,
            'Tc_midpoint_K': 42.0,
            'Critical_Field_Hc1_T': 0.05,
            'Critical_Field_Hc2_T': 15.2,
            'Magnetic_Susceptibility_emu_g': 0.95,
            'Heat_Capacity_Jump_J_molK': 125,
            'Isotope_Effect_alpha': 0.45,
            'Phonon_Frequencies_cm': '420-450',
            'EPC_Lambda': 0.72,
            'Notes': 'Clean_transition_sharp_onset',
            'Data_Quality': 'Excellent',
            'Review_Status': 'Approved',
            'Publication_Ready': 'Yes'
        },
        {
            'Sample_ID': 'MgFe-002',
            'Composition': 'Mg67Fe33',
            'Synthesis_Method': 'Arc_Melting',
            'Synthesis_Date': '2024-01-16',
            'Operator': 'J.Smith',
            'Purity_Mg_Percent': 99.9,
            'Purity_Fe_Percent': 99.95,
            'Crucible_Material': 'Tungsten',
            'Atmosphere': 'Argon',
            'Max_Temperature_C': 1750,
            'Cooling_Rate_K_min': 1.0,
            'Annealing_Temperature_C': 750,
            'Annealing_Time_hours': 24,
            'Final_Mass_g': 1.85,
            'Color_Appearance': 'Light_Gray',
            'Crystal_Structure_XRD': 'Layered',
            'Lattice_Parameters_a_b_c': '5.25_5.25_8.60',
            'Density_g_cm3': 3.8,
            'Resistivity_300K_ohm_cm': 0.032,
            'Resistivity_77K_ohm_cm': 0.024,
            'Resistivity_4K_ohm_cm': 0.0,
            'Tc_onset_K': 38.7,
            'Tc_zero_K': 38.1,
            'Tc_midpoint_K': 38.4,
            'Critical_Field_Hc1_T': 0.04,
            'Critical_Field_Hc2_T': 12.8,
            'Magnetic_Susceptibility_emu_g': 0.88,
            'Heat_Capacity_Jump_J_molK': 108,
            'Isotope_Effect_alpha': 0.48,
            'Phonon_Frequencies_cm': '400-430',
            'EPC_Lambda': 0.65,
            'Notes': 'Broad_transition_multiple_phases',
            'Data_Quality': 'Good',
            'Review_Status': 'Under_Review',
            'Publication_Ready': 'No'
        }
    ]
    
    # Create blank template rows
    blank_template = {key: '' for key in sample_data[0].keys()}
    blank_template['Sample_ID'] = 'TEMPLATE_ENTRY'
    blank_template['Composition'] = 'MgXFeY'
    
    # Add 10 blank rows for data entry
    for i in range(10):
        blank_row = {key: '' for key in sample_data[0].keys()}
        sample_data.append(blank_row)
    
    # Create DataFrame
    df = pd.DataFrame(sample_data)
    
    # Create instructions DataFrame
    instructions = {
        'Parameter': [
            'Sample_ID', 'Composition', 'Synthesis_Method', 'Synthesis_Date',
            'Operator', 'Purity_Mg_Percent', 'Purity_Fe_Percent', 'Crucible_Material',
            'Atmosphere', 'Max_Temperature_C', 'Cooling_Rate_K_min', 'Annealing_Temperature_C',
            'Annealing_Time_hours', 'Final_Mass_g', 'Color_Appearance', 'Crystal_Structure_XRD',
            'Lattice_Parameters_a_b_c', 'Density_g_cm3', 'Resistivity_300K_ohm_cm',
            'Resistivity_77K_ohm_cm', 'Resistivity_4K_ohm_cm', 'Tc_onset_K',
            'Tc_zero_K', 'Tc_midpoint_K', 'Critical_Field_Hc1_T', 'Critical_Field_Hc2_T',
            'Magnetic_Susceptibility_emu_g', 'Heat_Capacity_Jump_J_molK', 'Isotope_Effect_alpha',
            'Phonon_Frequencies_cm', 'EPC_Lambda', 'Notes', 'Data_Quality',
            'Review_Status', 'Publication_Ready'
        ],
        'Description': [
            'Unique identifier for each sample',
            'MgXFeY format with atomic percentages',
            'Arc_Melting / Induction_Melting / Powder_Metallurgy / Other',
            'YYYY-MM-DD format',
            'Name of person conducting synthesis',
            'Starting material purity for Mg (e.g. 99.9)',
            'Starting material purity for Fe (e.g. 99.95)',
            'Tungsten / Alumina / Graphite / Other',
            'Argon / Vacuum / Air / Other',
            'Maximum synthesis temperature in Celsius',
            'Cooling rate in K/min (CRITICAL: Use 1.0 K/min for RBT phase)',
            'Post-synthesis annealing temperature',
            'Duration of annealing in hours',
            'Final sample mass in grams',
            'Visual description of final sample',
            'X-ray diffraction determined structure',
            'Unit cell parameters in Angstrom (a_b_c format)',
            'Measured density in g/cm³',
            'Electrical resistivity at room temperature',
            'Electrical resistivity at liquid nitrogen temperature',
            'Electrical resistivity at liquid helium temperature',
            'Superconducting transition onset temperature',
            'Zero resistance temperature',
            'Midpoint of superconducting transition',
            'Lower critical field in Tesla',
            'Upper critical field in Tesla',
            'Magnetic susceptibility measurement',
            'Specific heat jump at Tc',
            'Isotope effect exponent',
            'Measured phonon frequency range',
            'Electron-phonon coupling constant',
            'Detailed observations and comments',
            'Excellent / Good / Fair / Poor',
            'Approved / Under_Review / Needs_Revision / Rejected',
            'Yes / No - Ready for publication'
        ],
        'Expected_Range': [
            'MgFe-XXX format',
            'Mg25Fe75 to Mg75Fe25',
            'Prefer Arc_Melting',
            'Current date',
            'Lab member name',
            '>99.9% required',
            '>99.9% required',
            'Tungsten preferred',
            'Argon preferred',
            '1700-1900°C',
            '1.0 K/min CRITICAL',
            '700-900°C',
            '12-48 hours',
            '1-5 grams typical',
            'Gray to dark metal',
            'Cubic/Layered expected',
            '4-6 Angstrom typical',
            '3-6 g/cm³',
            '0.01-0.1 Ω⋅cm',
            '0.005-0.05 Ω⋅cm',
            '0 (superconductor)',
            '35-45 K expected',
            '34-44 K expected',
            '35-45 K expected',
            '0.01-0.1 Tesla',
            '5-20 Tesla',
            'Negative (diamagnetic)',
            '50-200 J/mol⋅K',
            '0.3-0.6 typical',
            '320-470 cm⁻¹',
            '0.3-0.9 range',
            'Detailed notes',
            'Quality assessment',
            'Review status',
            'Publication ready?'
        ]
    }
    
    instructions_df = pd.DataFrame(instructions)
    
    # Critical parameters DataFrame
    critical_params = {
        'Parameter': [
            'Expected_Tc_Range',
            'Required_Cooling_Rate',
            'Purity_Requirements',
            'Testing_Temperature_Range',
            'Expected_Structure',
            'Interface_Effects',
            'Pressure_Sensitivity',
            'RBT_Enhancement',
            'Contact_Information'
        ],
        'Value': [
            '35-45 K based on theoretical predictions',
            '1.0 K/min for optimal RBT phase formation',
            '>99.9% to avoid oxide contamination',
            'Measure resistivity from 300K down to 4K',
            'Cubic or layered phases most promising',
            'Superconductivity may be enhanced at grain boundaries',
            'Test under slight compression if ambient fails',
            '1.5x enhancement expected from RBT theory',
            'Contact: rt-superconductors-project@research.org'
        ]
    }
    
    critical_df = pd.DataFrame(critical_params)
    
    # Save to Excel with multiple sheets
    output_file = Path("measurement_checklist_template.xlsx")
    
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # Main data sheet
        df.to_excel(writer, sheet_name='Measurement_Data', index=False)
        
        # Instructions sheet
        instructions_df.to_excel(writer, sheet_name='Instructions', index=False)
        
        # Critical parameters sheet
        critical_df.to_excel(writer, sheet_name='Critical_Parameters', index=False)
    
    print(f"✅ Excel template created: {output_file}")
    print("📋 Sheets included:")
    print("  • Measurement_Data - Main data entry")
    print("  • Instructions - Detailed parameter descriptions")
    print("  • Critical_Parameters - Key experimental requirements")
    
    return output_file

if __name__ == "__main__":
    create_excel_template() 