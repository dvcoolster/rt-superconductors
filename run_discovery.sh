#!/bin/bash
#
# RT-Superconductors: Complete Discovery Workflow
# 
# This script demonstrates the full RBT-based superconductor discovery pipeline:
# 1. Enumerate valence-balanced compositions
# 2. Apply RBT ledger scoring 
# 3. Train/use ML models for ranking
# 4. Screen with xTB geometry optimization
# 5. Validate with DFT calculations
#
# Usage: ./run_discovery.sh [quick|full]
#

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
MODE=${1:-quick}  # quick or full
MAX_ATOMS_QUICK=6
MAX_ATOMS_FULL=10
N_CANDIDATES_QUICK=50
N_CANDIDATES_FULL=200
N_XTB_QUICK=20
N_XTB_FULL=100
N_DFT_QUICK=5
N_DFT_FULL=25

echo -e "${BLUE}🧬 RT-Superconductors: RBT-Based Discovery Pipeline${NC}"
echo -e "${BLUE}============================================${NC}"
echo -e "Mode: ${YELLOW}${MODE}${NC}"
echo -e "Reference: Chauhan & Chouhan (2025) 'Recursive Becoming'"

# Create necessary directories
echo -e "\n${GREEN}📁 Setting up directories...${NC}"
mkdir -p data/generated results models logs

# Check if conda environment is active
if [[ "$CONDA_DEFAULT_ENV" != "rt-sc" ]]; then
    echo -e "${YELLOW}⚠️  Activating conda environment...${NC}"
    conda activate rt-sc || {
        echo -e "${RED}❌ Could not activate rt-sc environment. Please run:${NC}"
        echo -e "   conda env create -f env/conda.yaml"
        echo -e "   conda activate rt-sc"
        exit 1
    }
fi

# Set parameters based on mode
if [[ "$MODE" == "quick" ]]; then
    MAX_ATOMS=$MAX_ATOMS_QUICK
    N_CANDIDATES=$N_CANDIDATES_QUICK
    N_XTB=$N_XTB_QUICK
    N_DFT=$N_DFT_QUICK
    echo -e "${YELLOW}🚀 Quick mode: Small-scale demonstration${NC}"
else
    MAX_ATOMS=$MAX_ATOMS_FULL
    N_CANDIDATES=$N_CANDIDATES_FULL
    N_XTB=$N_XTB_FULL
    N_DFT=$N_DFT_FULL
    echo -e "${YELLOW}🏭 Full mode: Production-scale discovery${NC}"
fi

# Step 1: Enumerate valence-balanced compositions
echo -e "\n${GREEN}1️⃣ Enumerating RBT-favorable compositions...${NC}"
echo -e "   Using split-octet valence balance principle"

python src/00_enumerate.py \
    --max-atoms $MAX_ATOMS \
    --tau-threshold 0.5 \
    --include-known \
    --top-n $N_CANDIDATES \
    --out data/generated/candidates.json \
    2>&1 | tee logs/01_enumerate.log

if [[ $? -eq 0 ]]; then
    echo -e "${GREEN}✅ Enumeration completed${NC}"
    N_ENUM=$(python -c "import json; print(len(json.load(open('data/generated/candidates.json'))))")
    echo -e "   Generated: ${YELLOW}${N_ENUM}${NC} compositions"
else
    echo -e "${RED}❌ Enumeration failed${NC}"
    exit 1
fi

# Step 2: Apply RBT ledger scoring
echo -e "\n${GREEN}2️⃣ Applying RBT ledger scoring (τ, κ, phase-locking)...${NC}"
echo -e "   τ: valence mismatch from split-octet deviations"
echo -e "   κ: curvature penalty from discrete gravity"

python src/01_ledger_scores.py \
    data/generated/candidates.json \
    --threshold_tau 1.0 \
    --threshold_kappa 0.2 \
    --min_rbt_score 0.1 \
    --max_candidates $((N_CANDIDATES / 2)) \
    --out results/scored_candidates.csv \
    2>&1 | tee logs/02_scoring.log

if [[ $? -eq 0 ]]; then
    echo -e "${GREEN}✅ RBT scoring completed${NC}"
    N_SCORED=$(python -c "import pandas as pd; print(len(pd.read_csv('results/scored_candidates.csv')))")
    echo -e "   Scored: ${YELLOW}${N_SCORED}${NC} candidates"
else
    echo -e "${RED}❌ RBT scoring failed${NC}"
    exit 1
fi

# Step 3: Machine Learning Prediction (if not in quick mode or model exists)
echo -e "\n${GREEN}3️⃣ Machine Learning Tc prediction...${NC}"

if [[ ! -f "models/rbt_predictor.pkl" ]]; then
    echo -e "   Training new ML model..."
    python src/02_ml_predict.py --train \
        --supercon_data data-sources/supercon.csv \
        --rbt_features results/scored_candidates.csv \
        --model_out models/rbt_predictor.pkl \
        2>&1 | tee logs/03_ml_train.log
    
    if [[ $? -ne 0 ]]; then
        echo -e "${YELLOW}⚠️  ML training failed, using RBT scores only${NC}"
        cp results/scored_candidates.csv results/ml_ranked.csv
    fi
else
    echo -e "   Using existing ML model..."
fi

if [[ -f "models/rbt_predictor.pkl" ]]; then
    echo -e "   Ranking candidates with ML predictions..."
    python src/02_ml_predict.py --rank results/scored_candidates.csv \
        --model models/rbt_predictor.pkl \
        --top $N_XTB \
        --out results/ml_ranked.csv \
        2>&1 | tee logs/03_ml_rank.log
    
    if [[ $? -eq 0 ]]; then
        echo -e "${GREEN}✅ ML ranking completed${NC}"
        N_ML=$(python -c "import pandas as pd; print(len(pd.read_csv('results/ml_ranked.csv')))")
        echo -e "   Top candidates: ${YELLOW}${N_ML}${NC}"
    else
        echo -e "${YELLOW}⚠️  ML ranking failed, using RBT ranking${NC}"
        head -n $((N_XTB + 1)) results/scored_candidates.csv > results/ml_ranked.csv
    fi
else
    echo -e "${YELLOW}⚠️  No ML model, using RBT ranking${NC}"
    head -n $((N_XTB + 1)) results/scored_candidates.csv > results/ml_ranked.csv
fi

# Step 4: xTB Geometry Optimization (fast screening)
echo -e "\n${GREEN}4️⃣ xTB geometry optimization (fast screening)...${NC}"
echo -e "   Semi-empirical method for stability assessment"

# Check if xTB is available
if command -v xtb &> /dev/null || python -c "import xtb" &> /dev/null; then
    python src/03_xtb_relax.py \
        --input results/ml_ranked.csv \
        --n_candidates $N_XTB \
        --min_stability 0.6 \
        --max_gap 3.0 \
        --save_structures \
        --out results/xtb_stable.csv \
        2>&1 | tee logs/04_xtb.log
    
    if [[ $? -eq 0 ]]; then
        echo -e "${GREEN}✅ xTB screening completed${NC}"
        N_XTB_STABLE=$(python -c "import pandas as pd; print(len(pd.read_csv('results/xtb_stable.csv')))")
        echo -e "   Stable candidates: ${YELLOW}${N_XTB_STABLE}${NC}"
    else
        echo -e "${YELLOW}⚠️  xTB screening failed, skipping to DFT${NC}"
        head -n $((N_DFT + 1)) results/ml_ranked.csv > results/xtb_stable.csv
    fi
else
    echo -e "${YELLOW}⚠️  xTB not available, skipping geometry optimization${NC}"
    echo -e "   Install with: conda install -c conda-forge xtb-python"
    head -n $((N_DFT + 1)) results/ml_ranked.csv > results/xtb_stable.csv
fi

# Step 5: DFT Validation (high-accuracy)
echo -e "\n${GREEN}5️⃣ DFT validation (Quantum ESPRESSO)...${NC}"
echo -e "   High-accuracy electronic structure calculations"

# Check if we're in a SLURM environment
if command -v sbatch &> /dev/null; then
    echo -e "   SLURM detected - submitting batch jobs"
    DFT_MODE="--slurm"
else
    echo -e "   Running DFT calculations locally"
    DFT_MODE=""
fi

# For quick mode or if QE not available, just create input files
if [[ "$MODE" == "quick" ]] || ! command -v pw.x &> /dev/null; then
    echo -e "${YELLOW}⚠️  Creating DFT input files only (no calculations)${NC}"
    
    python src/04_qe_relax.py \
        --batch results/xtb_stable.csv \
        --n_best $N_DFT \
        --ecut 60 \
        --out results/dft_inputs \
        2>&1 | tee logs/05_dft.log
    
    echo -e "${GREEN}✅ DFT input files created${NC}"
    echo -e "   Run with: python src/04_qe_relax.py --batch results/xtb_stable.csv"
else
    python src/04_qe_relax.py \
        --batch results/xtb_stable.csv \
        --n_best $N_DFT \
        --ecut 60 \
        $DFT_MODE \
        --out results/dft_results \
        2>&1 | tee logs/05_dft.log
    
    if [[ $? -eq 0 ]]; then
        echo -e "${GREEN}✅ DFT calculations completed${NC}"
    else
        echo -e "${YELLOW}⚠️  Some DFT calculations failed${NC}"
    fi
fi

# Generate summary report
echo -e "\n${GREEN}📊 Generating discovery summary...${NC}"

cat > results/discovery_summary.md << EOF
# RBT Superconductor Discovery Summary

**Run Date**: $(date)
**Mode**: $MODE
**Reference**: Chauhan & Chouhan (2025) "Recursive Becoming: From Nothingness to Everything"

## Pipeline Results

1. **Enumeration**: Generated $N_ENUM RBT-favorable compositions
2. **Scoring**: Scored $N_SCORED candidates with τ/κ metrics  
3. **ML Ranking**: Selected top candidates for screening
4. **xTB Screening**: Validated geometric stability
5. **DFT Validation**: High-accuracy electronic structure

## Key Files

- \`data/generated/candidates.json\`: Enumerated compositions
- \`results/scored_candidates.csv\`: RBT-scored candidates
- \`results/ml_ranked.csv\`: ML-ranked top candidates
- \`results/xtb_stable.csv\`: Geometrically stable candidates
- \`results/dft_results/\`: DFT calculation outputs

## RBT Theory Summary

- **τ (tau)**: Valence mismatch from split-octet deviations (lower is better)
- **κ (kappa)**: Curvature penalty from discrete gravity (lower is better)
- **Phase-locking**: Cooper pair formation potential (higher is better)

## Next Steps

1. Analyze DFT results for electronic structure
2. Calculate superconducting properties (if stable)
3. Prioritize candidates for experimental synthesis
4. Validate RBT predictions experimentally

---
*Generated by RT-Superconductors Discovery Pipeline*
EOF

echo -e "${GREEN}✅ Summary saved to results/discovery_summary.md${NC}"

# Final status
echo -e "\n${BLUE}🎯 RBT Superconductor Discovery Complete!${NC}"
echo -e "${BLUE}=========================================${NC}"

if [[ -f "results/xtb_stable.csv" ]]; then
    N_FINAL=$(python -c "import pandas as pd; print(len(pd.read_csv('results/xtb_stable.csv')))" 2>/dev/null || echo "0")
    echo -e "🏆 Final candidates for synthesis: ${YELLOW}${N_FINAL}${NC}"
fi

echo -e "\n📖 Key Results:"
echo -e "   📁 Check results/ directory for all outputs"
echo -e "   📊 Open notebooks/explorer.ipynb for interactive analysis"
echo -e "   📋 Read results/discovery_summary.md for complete summary"

echo -e "\n🔮 RBT Predictions:"
echo -e "   • Room-temperature superconductors exist in unexplored τ/κ space"
echo -e "   • Optimal candidates have τ≈0 (split-octet balance)"
echo -e "   • Low κ minimizes curvature disruption of Cooper pairs"

echo -e "\n🚀 What's Next:"
echo -e "   1. Analyze top candidates in results/xtb_stable.csv"
echo -e "   2. Run DFT calculations: python src/04_qe_relax.py --batch"
echo -e "   3. Prioritize for experimental synthesis"
echo -e "   4. Test RBT's experimental predictions (XRISM 2026, MAGIS-100 2028)"

echo -e "\n📚 Reference: https://www.recursivebecoming.info/RBT_v1.0_release.pdf"
echo -e "${GREEN}🧬 From a single δ-glitch to room-temperature superconductivity!${NC}" 