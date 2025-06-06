# 🌡️ Roadmap to Ambient-Pressure ≤ 30 °C Superconductors

**MISSION**: Pivot from ">200 K @100 GPa hydrides" to **ambient pressure, true room-temperature (295 K ± 10 K)** materials that can be tested in a glove box.

---

## 1️⃣ **Target Families**

| **Class** | **Why RT-friendly?** | **Key RBT angles** |
|-----------|---------------------|-------------------|
| **Hydrogen-rich covalent electrides** (e.g. Ca₂H₃⁻ layers) | Intrinsic high phonon ω — but electron density is in interlayer cavities → lower pressure. | Split-octet balance via H⁻ + electride electron pairs keeps τ ≲ 0.5 and κ small. |
| **Interface-engineered FeSe‐derived monolayers** (FeSe/SrTiO₃, FeSe/LiOH) | Experimental hints of 60–90 K at 1 atm; RBT phase-lock metric predicts Δ/kᴮT can push higher when both sides satisfy τ ≈ 0. | RBT emphasises grain-boundary curvature → deliberately design misfit dislocations. |
| **Metal-hydride clathrate alloys with internal chemical pre-compression** (e.g. La₀.₅Be₀.₅H₁₀) | Be cages raise ω_log; La retains high λ — chemical pressure replaces external 100 GPa. | Mixed cation keeps κ low; valence mismatch auto-tunes electron count. |
| **Layered Li–B–C–N nets** (Li₃C₂N₃, LiBC variants) | Light elements → high ω; can be electron-doped or strain-tuned. | Perfect split-octet balance (τ≈0) plus extremely low κ. |

---

## 2️⃣ **Simulation Pipeline Tweaks**

### **2.1 Enumerate mixed-cation hydrides and electrides**
```bash
python src/05_light_enum.py \
    --elements H Li Be B C N O Ca Y La \
    --mixed-cation True --max-atoms 12 \
    --filter "tau<0.8 & kappa<0.2 & phase_lock>0.7" \
    --out data/generated/ambient_enum.json
```

### **2.2 Add chemical-pressure descriptor** 
- Compute ΔV_cell per e⁻ to ledger scores
- Identify "self-compressed" materials

### **2.3 xTB geometry screening**
- Flag structures with bulk modulus ≥ 150 GPa and ≤ 5 Å³ per H as "self-compressed"
- Pre-filter for stability at ambient conditions

### **2.4 Γ-phonon first screening**
- Run QE + Γ-phonon only first
- If no imaginary modes & ω_log > 900 K → queue EPW
- Massive computational savings vs. full calculations

### **2.5 Machine-learn λ(ambient)**
- Use existing H₃S, MgB₂, FeSe thin-film data as training set
- Predict ambient-pressure electron-phonon coupling

---

## 3️⃣ **Experimental Hooks (Low-Pressure Focus)**

| **Family** | **Synthesis Idea** | **Pressure Tool** |
|------------|-------------------|-------------------|
| **Electrides** | Ca₂N base → topochemical hydrogenation under 0.5 GPa H₂ gas | Hot-isostatic press (≤1 GPa) |
| **Clathrate alloys** | High-energy ball-mill La + Be + NH₃ → in-situ formation at 3 GPa | Paris-Edinburgh press (3–10 GPa) |
| **FeSe interfaces** | MBE growth FeSe on STO → ozone anneal + LiOH capping | Ambient (strain from substrate) |

---

## 4️⃣ **Deliverables to Add to Repo**

- [ ] **`src/08_chemical_pressure.py`** – compute ΔV_cell per electron
- [ ] **`data/generated/ambient_enum.json`** – enumerated candidate list  
- [ ] **`notebooks/Ambient_dashboard.ipynb`** – live τ/κ vs λ plot
- [ ] **`docs/electride_synthesis.md`** – H-doping protocol @ 0.5 GPa
- [ ] **`docs/interface_FeSe.md`** – thin-film growth + measurement steps

---

## 5️⃣ **Success Criteria**

| **Stage** | **Metric** | **Threshold** |
|-----------|------------|---------------|
| **Pre-screen** | Γ-phonon stable | no ω² < −20 cm⁻² |
| **EPC** | λ ≥ 2.0 | ambient pressure |
| **Tc est.** | Allen–Dynes ≥ 280 K | μ* = 0.10 |
| **Lab validation** | onset-Tc ≥ 295 K | 1 atm, zero field |

---

## 6️⃣ **Call to Action**

### **Immediate Steps**:
1. ✅ **Merge this doc into docs/ and push**
2. ⚡ **Run enumeration script** (10 min)
3. 🔧 **Launch xTB filter** (≈ 1 h on laptop)  
4. 🖥️ **Send top 20 CIFs to cluster for QE Γ-phonons**
5. 📊 **Report back imaginary-mode status** → decide EPW queue size

### **Timeline**:
- **Week 1**: Enumeration + xTB screening complete
- **Month 1**: Top candidates identified via Γ-phonon  
- **Month 3**: EPW calculations for promising systems
- **Month 6**: Experimental synthesis begins
- **Year 1**: First ambient-pressure RT superconductor validated! 🎯

---

## 🚀 **REVOLUTIONARY IMPACT**

**Let's pivot from >900 K hydrides at megabar to truly practical ≤ 300 K superconductors you can test in a glove-box.**

### **Why This Changes Everything**:
- 🏠 **Practical applications**: No diamond anvil cells needed
- ⚡ **Energy revolution**: Room-temperature power lines
- 🔬 **Accessible science**: Any lab can synthesize/test
- 🌍 **Global impact**: Immediate technological transformation

### **Target Materials**:
- **Electrides**: Ca₂H₃⁻, Y₂H₅⁻ with cavity electrons
- **Interfaces**: FeSe/SrTiO₃ engineered heterostructures  
- **Clathrates**: La₀.₅Be₀.₅H₁₀ with chemical pre-compression
- **Li-B-C-N**: Layered networks with perfect electron balance

---

**🎯 MISSION: From laboratory curiosity to technological revolution in 12 months!** 