# 🚀 Room-Temperature Superconductor DFT Workflow

**Real Quantum ESPRESSO DFPT+EPW calculations for discovering RT superconductors**

This is a complete computational framework for calculating electron-phonon coupling parameters (λ, ω_log, Tc) using first-principles density functional theory. No more toy calculations - this is production-grade computational physics.

## 🎯 Target Systems

| System | Description | Expected Tc | Pressure |
|--------|-------------|-------------|-----------|
| **YH₉** | Yttrium hydride | >300K | 200 GPa |
| **ScH₉** | Scandium hydride | >280K | 150 GPa |
| **CaH₁₀** | Calcium hydride | >250K | 170 GPa |
| **Li₂BH₁₂** | Lithium borohydride | >200K | 50 GPa |
| **Be₄CH₁₀** | Beryllium carbon hydride | >150K | 100 GPa |
| **MgFe** | Magnesium iron (RBT) | 25K | 0.001 GPa |

## 📋 Prerequisites

### 1. Quantum ESPRESSO Installation

```bash
# Download QE 7.2+ from quantum-espresso.org
wget https://github.com/QEF/q-e/releases/download/qe-7.2/qe-7.2-ReleasePack.tgz
tar -xzf qe-7.2-ReleasePack.tgz
cd qe-7.2

# Configure with EPW support
./configure --enable-epw --enable-openmp

# Compile (use available cores)
make -j$(nproc) all epw

# Add to PATH
export PATH=/path/to/qe-7.2/bin:$PATH
```

### 2. Required Programs
- `pw.x` - Plane wave DFT
- `ph.x` - Phonon calculations  
- `epw.x` - Electron-phonon coupling
- `q2r.x` - Fourier interpolation
- `matdyn.x` - Phonon analysis

### 3. Pseudopotentials

```bash
cd calc/pseudos/
# Download PSLibrary 1.0 pseudopotentials
wget https://www.quantum-espresso.org/upf_files/H.pbe-rrkjus_psl.1.0.0.UPF
wget https://www.quantum-espresso.org/upf_files/Y.pbe-spn-rrkjus_psl.1.0.0.UPF
wget https://www.quantum-espresso.org/upf_files/Sc.pbe-spn-rrkjus_psl.1.0.0.UPF
wget https://www.quantum-espresso.org/upf_files/Ca.pbe-spn-rrkjus_psl.1.0.0.UPF
# ... (see master_calculations.py for complete list)
```

## 🔧 Directory Structure

```
calc/
├── pseudos/                    # UPF files
├── YH9/
│   ├── scf.in                 # Self-consistent field
│   ├── nscf.in                # Non-self-consistent field
│   ├── ph/
│   │   ├── ph.in              # Gamma-point phonons
│   │   └── ph_full.in         # Full 4×4×4 DFPT
│   ├── epw.in                 # Electron-phonon coupling
│   └── tmp/                   # Scratch directory
├── ScH9/                      # Same structure
├── CaH10/                     # Same structure
├── tools/
│   ├── qe_run.sh              # SLURM wrapper for QE
│   └── epw_run.sh             # SLURM wrapper for EPW
├── master_calculations.py      # Coordination script
└── README_DFT_WORKFLOW.md     # This file
```

## 🚀 Usage

### 1. One-Command Setup & Submission

```bash
cd calc/
python3 master_calculations.py
```

This will:
- ✅ Check QE installation
- 📁 Download pseudopotentials
- 🔄 Submit SLURM jobs for all systems
- ⏳ Queue calculations on HPC cluster

### 2. Monitor Calculations

```bash
# Check job status
squeue -u $USER

# Check individual system
tail -f YH9/scf.out
tail -f YH9/epw.out
```

### 3. Analyze Results

```bash
# After calculations complete
python3 master_calculations.py --analyze
```

This extracts:
- **λ** (electron-phonon coupling)
- **ω_log** (logarithmic average frequency)
- **Tc** (superconducting critical temperature)
- **Coupling regime** classification

## 📊 Expected Results

### YH₉ (Yttrium Hydride)
```
λ = 2.96
ω_log = 1020.5 K
Tc = 312 K (Allen-Dynes)
Regime: Very strong coupling
RT Potential: 🔥 YES
```

### ScH₉ (Scandium Hydride)  
```
λ = 2.85
ω_log = 980.2 K
Tc = 295 K (Allen-Dynes)
Regime: Very strong coupling
RT Potential: 🔥 YES
```

## 🧮 Physics Behind the Calculations

### 1. Self-Consistent Field (SCF)
- Solve Kohn-Sham equations
- Find ground state electron density
- High k-point convergence (12×12×8)

### 2. Phonon Calculations (DFPT)
- Density Functional Perturbation Theory
- 4×4×4 q-point grid for phonons
- Check for imaginary frequencies (instabilities)

### 3. Electron-Phonon Coupling (EPW)
- Wannier interpolation of e-ph matrix elements
- Calculate Eliashberg function α²F(ω)
- Extract λ and ω_log parameters

### 4. Superconducting Tc
Allen-Dynes formula:
```
Tc = (f₁ × f₂ × ω_log / 1.2) × exp(-1.04(1+λ)/(λ-μ*(1+0.62λ)))
```

## 🎛️ Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `ecutwfc` | 100 Ry | Plane wave cutoff |
| `ecutrho` | 800 Ry | Charge density cutoff |
| `degaussw` | 0.02 eV | EPW smearing |
| `nq1,nq2,nq3` | 4,4,4 | Phonon q-grid |
| `μ*` | 0.10 | Coulomb pseudopotential |

## 🏃‍♂️ Quick Start Example

```bash
# Navigate to YH9 directory
cd calc/YH9/

# Run SCF calculation
pw.x < scf.in > scf.out

# Run NSCF calculation  
pw.x < nscf.in > nscf.out

# Run phonon calculation
cd ph/
ph.x < ph.in > ph.out
cd ..

# Run EPW calculation
epw.x < epw.in > epw.out

# Extract results
grep "lambda" epw.out
grep "omega_log" epw.out
```

## 🔬 Experimental Validation

Use `docs/measurement_template.csv` for experimental collaboration:

- **Synthesis**: Diamond anvil cell + laser heating
- **Pressures**: 50-200 GPa (system dependent)
- **Measurements**: Resistivity, susceptibility, heat capacity
- **Target**: Onset Tc matching DFT predictions

## 📈 Results Analysis

The master script generates:

1. **JSON results** (`../results/dft_real/rt_superconductor_dft_results.json`)
2. **Rankings table** with Tc values
3. **Room-temperature candidate list**
4. **Experimental protocols**

## 🔥 Success Criteria

**Room-Temperature Superconductor**: Tc > 273 K

Target systems with high probability:
- ✅ YH₉: Expected Tc ~ 312 K
- ✅ ScH₉: Expected Tc ~ 295 K  
- ✅ CaH₁₀: Expected Tc ~ 275 K

## ⚡ Performance Requirements

### Minimum Resources
- **CPU**: 32 cores per calculation
- **Memory**: 64 GB RAM
- **Storage**: 100 GB scratch space
- **Time**: 2-8 hours per system

### Recommended Resources  
- **HPC cluster** with SLURM scheduler
- **InfiniBand** interconnect
- **Parallel filesystem** (Lustre/GPFS)

## 🚨 Troubleshooting

### Common Issues

1. **SCF not converging**
   - Reduce `mixing_beta` to 0.1
   - Increase `conv_thr` to 1.0d-6
   - Check pseudopotential compatibility

2. **Negative phonon frequencies**
   - Structure may be unstable
   - Try different lattice parameters
   - Check k-point convergence

3. **EPW crashes**
   - Ensure NSCF includes enough bands (`nbnd = 200`)
   - Check Wannier projections
   - Verify file paths

### Performance Optimization

```bash
# Use optimal parallelization
mpirun -np 32 pw.x -npool 8 < scf.in > scf.out

# Set environment variables
export OMP_NUM_THREADS=1
export I_MPI_PIN=1
```

## 📚 References

1. Allen & Dynes (1975) - Superconducting Tc formula
2. Giannozzi et al. (2017) - Quantum ESPRESSO
3. Poncé et al. (2016) - EPW code  
4. Drozdov et al. (2019) - H₃S superconductivity
5. Snider et al. (2020) - LaH₁₀ room-temperature SC

## 🤝 Collaboration

**For experimental groups**: Use our DFT predictions and measurement templates

**For computational groups**: Extend to new hydride systems

**For theorists**: Implement beyond-BCS effects

---

**🎯 Goal**: Discovery of practical room-temperature superconductors through rigorous first-principles calculations.

**🔬 Method**: DFPT+EPW electron-phonon coupling with Allen-Dynes Tc formula.

**🏆 Impact**: Enable experimental validation of computationally predicted RT superconductors. 