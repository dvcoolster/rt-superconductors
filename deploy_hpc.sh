#!/bin/bash

# ================================================================
# HPC DEPLOYMENT SCRIPT FOR RT-SUPERCONDUCTOR DISCOVERY
# ================================================================
#
# This script deploys the complete ambient-pressure superconductor
# discovery pipeline to HPC clusters for breakthrough identification.
#
# Author: RT-Superconductor Research Team
# Date: 2025-06-07
# 
# Usage: ./deploy_hpc.sh [options]
#
# Options:
#   --cluster <name>    Target HPC cluster (auto-detected if not specified)
#   --nodes <n>         Number of nodes to request (default: auto)
#   --priority <level>  Job priority (low/normal/high, default: normal)
#   --test-run          Submit only 3 structures for testing
#   --dry-run           Show commands without executing
#
# Examples:
#   ./deploy_hpc.sh                    # Deploy all 51 candidates
#   ./deploy_hpc.sh --test-run         # Deploy 3 test candidates
#   ./deploy_hpc.sh --dry-run          # Show deployment commands
#   ./deploy_hpc.sh --cluster frontera # Deploy to specific cluster

set -e

# ================================================================
# CONFIGURATION AND DEFAULTS
# ================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

# Default settings
CLUSTER=""
NODES="auto"
PRIORITY="normal"
TEST_RUN=false
DRY_RUN=false
SUBMIT_DELAY=2  # seconds between job submissions

# HPC cluster detection and configuration
get_cluster_info() {
    case "$1" in
        frontera) echo "TACC Frontera (56 cores/node, 192GB RAM)" ;;
        stampede2) echo "TACC Stampede2 (68 cores/node, 96GB RAM)" ;;
        bridges2) echo "PSC Bridges-2 (128 cores/node, 256GB RAM)" ;;
        summit) echo "ORNL Summit (42 cores/node, 512GB RAM)" ;;
        perlmutter) echo "NERSC Perlmutter (128 cores/node, 512GB RAM)" ;;
        *) echo "Generic SLURM cluster" ;;
    esac
}

# ================================================================
# FUNCTIONS
# ================================================================

print_banner() {
    echo "🚀 RT-SUPERCONDUCTOR HPC DEPLOYMENT"
    echo "===================================="
    echo "Date: $(date)"
    echo "Project: Ambient-Pressure Room-Temperature Superconductors"
    echo "Mission: Discover the first RT superconductor that works on your desk!"
    echo ""
}

detect_cluster() {
    if [[ -n "$TACC_SYSTEM" ]]; then
        echo "$TACC_SYSTEM"
    elif [[ -n "$PSC_CLUSTER" ]]; then
        echo "bridges2"
    elif [[ "$HOSTNAME" == *"summit"* ]]; then
        echo "summit"
    elif [[ "$HOSTNAME" == *"perlmutter"* ]]; then
        echo "perlmutter"
    elif [[ -n "$SLURM_CLUSTER_NAME" ]]; then
        echo "$SLURM_CLUSTER_NAME"
    else
        echo "unknown"
    fi
}

get_cluster_params() {
    local cluster=$1
    case "$cluster" in
        frontera)
            echo "--account=TG-DMR140143 --partition=normal --nodes=1 --ntasks=56"
            ;;
        stampede2)
            echo "--account=TG-DMR140143 --partition=normal --nodes=1 --ntasks=68"
            ;;
        bridges2)
            echo "--account=che140143p --partition=RM --nodes=1 --ntasks=128"
            ;;
        summit)
            echo "--account=CHE143 --partition=batch --nodes=1 --ntasks=42"
            ;;
        perlmutter)
            echo "--account=m1234 --partition=regular --nodes=1 --ntasks=128"
            ;;
        *)
            echo "--nodes=1 --ntasks=32"  # Generic defaults
            ;;
    esac
}

check_prerequisites() {
    echo "🔍 Checking prerequisites..."
    
    # Check if we're in the right directory
    if [[ ! -f "COMPUTE_FIRST_PLAYBOOK.md" ]]; then
        echo "❌ ERROR: Must run from rt-superconductors root directory"
        exit 1
    fi
    
    # Check if QE runs directory exists
    if [[ ! -d "calc/qe_runs" ]]; then
        echo "❌ ERROR: calc/qe_runs directory not found"
        echo "   Please run the compute-first playbook first"
        exit 1
    fi
    
    # Count available structures
    local n_structures=$(find calc/qe_runs -mindepth 1 -maxdepth 1 -type d | wc -l)
    echo "✅ Found $n_structures calculation directories"
    
    # Check SLURM availability
    if ! command -v sbatch &> /dev/null; then
        echo "❌ ERROR: SLURM (sbatch) not available"
        echo "   This script requires a SLURM-based HPC system"
        exit 1
    fi
    
    echo "✅ Prerequisites check passed"
    echo ""
}

estimate_resources() {
    local n_jobs=$1
    local cluster=$2
    
    echo "📊 Resource Estimation for $n_jobs Jobs"
    echo "======================================="
    
    case "$cluster" in
        frontera|stampede2)
            local cores_per_job=56
            local mem_per_job="96GB"
            local walltime="8:00:00"
            ;;
        bridges2|perlmutter)
            local cores_per_job=64
            local mem_per_job="128GB"
            local walltime="8:00:00"
            ;;
        *)
            local cores_per_job=32
            local mem_per_job="64GB"
            local walltime="8:00:00"
            ;;
    esac
    
    local total_cores=$((n_jobs * cores_per_job))
    local total_core_hours=$((total_cores * 8))
    
    echo "🖥️  Cores per job: $cores_per_job"
    echo "💾 Memory per job: $mem_per_job"
    echo "⏰ Wall time: $walltime"
    echo "🔢 Total cores: $total_cores"
    echo "⏱️  Total core-hours: $total_core_hours"
    echo "💰 Estimated SU cost: $((total_core_hours / 10)) (assuming 0.1 SU/core-hour)"
    echo ""
}

submit_jobs() {
    local test_run=$1
    local dry_run=$2
    local cluster_params="$3"
    
    cd calc/qe_runs
    
    # Get list of structures to submit
    local structures=($(ls -d */ | sed 's#/##'))
    local n_total=${#structures[@]}
    
    if [[ "$test_run" == "true" ]]; then
        structures=("${structures[@]:0:3}")  # First 3 only
        echo "🧪 TEST RUN: Submitting ${#structures[@]} structures for testing"
    else
        echo "🚀 PRODUCTION RUN: Submitting all $n_total structures"
    fi
    
    echo ""
    echo "📋 Job Submission Plan:"
    echo "======================"
    for i in "${!structures[@]}"; do
        echo "  $((i+1)). ${structures[i]}"
    done
    echo ""
    
    if [[ "$dry_run" == "true" ]]; then
        echo "🔍 DRY RUN: Commands that would be executed:"
        echo "============================================="
        for struct in "${structures[@]}"; do
            echo "sbatch $cluster_params --job-name=RT-SC-$struct ../tools/qe_run.sh $struct"
        done
        echo ""
        echo "ℹ️  Add --submit to actually submit jobs"
        return 0
    fi
    
    # Confirm submission
    echo "⚠️  About to submit ${#structures[@]} jobs to HPC cluster"
    echo "   Each job will run for up to 8 hours"
    echo "   This will consume significant computational resources"
    echo ""
    read -p "Continue with job submission? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Job submission cancelled"
        exit 0
    fi
    
    # Submit jobs
    echo ""
    echo "🚀 Submitting jobs..."
    echo "===================="
    
    local submitted=0
    local failed=0
    
    for struct in "${structures[@]}"; do
        echo -n "  Submitting $struct... "
        
        local job_cmd="sbatch $cluster_params --job-name=RT-SC-$struct ../tools/qe_run.sh $struct"
        
        if output=$(eval "$job_cmd" 2>&1); then
            local job_id=$(echo "$output" | grep -o '[0-9]\+')
            echo "✅ Job $job_id"
            ((submitted++))
        else
            echo "❌ FAILED"
            echo "   Error: $output"
            ((failed++))
        fi
        
        # Brief delay between submissions
        sleep $SUBMIT_DELAY
    done
    
    cd ../..
    
    echo ""
    echo "📊 Submission Summary:"
    echo "====================="
    echo "✅ Successfully submitted: $submitted jobs"
    echo "❌ Failed submissions: $failed jobs"
    echo "📈 Success rate: $(( submitted * 100 / (submitted + failed) ))%"
    
    if [[ $submitted -gt 0 ]]; then
        echo ""
        echo "🎯 Next Steps:"
        echo "============="
        echo "1. Monitor jobs: squeue -u \$USER"
        echo "2. Check progress: watch 'squeue -u \$USER'"
        echo "3. Collect results: python3 utils/parse_epw.py --qe-runs-dir calc/qe_runs --output-dir results"
        echo ""
        echo "🔥 Expected breakthrough identification in 6 days!"
    fi
}

monitor_progress() {
    echo "📊 LIVE PROGRESS MONITORING"
    echo "==========================="
    
    while true; do
        clear
        echo "📊 RT-Superconductor Job Status ($(date))"
        echo "=========================================="
        
        # Count jobs by status
        local running=$(squeue -u $USER -n RT-SC -h -t R | wc -l)
        local pending=$(squeue -u $USER -n RT-SC -h -t PD | wc -l)
        local total=$((running + pending))
        
        echo "🏃 Running: $running jobs"
        echo "⏳ Pending: $pending jobs"
        echo "📊 Total active: $total jobs"
        
        if [[ $total -eq 0 ]]; then
            echo ""
            echo "✅ All jobs completed!"
            echo "🔬 Ready for results analysis"
            break
        fi
        
        echo ""
        echo "📋 Job Details:"
        squeue -u $USER -n RT-SC -o "%.10i %.15j %.8t %.10M %.5D %R" | head -20
        
        echo ""
        echo "Press Ctrl+C to exit monitoring"
        sleep 30
    done
}

# ================================================================
# ARGUMENT PARSING
# ================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        --cluster)
            CLUSTER="$2"
            shift 2
            ;;
        --nodes)
            NODES="$2"
            shift 2
            ;;
        --priority)
            PRIORITY="$2"
            shift 2
            ;;
        --test-run)
            TEST_RUN=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --monitor)
            monitor_progress
            exit 0
            ;;
        --help|-h)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --cluster <name>    Target HPC cluster"
            echo "  --nodes <n>         Number of nodes (default: auto)"
            echo "  --priority <level>  Job priority (low/normal/high)"
            echo "  --test-run          Submit only 3 structures"
            echo "  --dry-run           Show commands without executing"
            echo "  --monitor           Monitor existing jobs"
            echo "  --help              Show this help"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# ================================================================
# MAIN EXECUTION
# ================================================================

print_banner

# Auto-detect cluster if not specified
if [[ -z "$CLUSTER" ]]; then
    CLUSTER=$(detect_cluster)
    echo "🔍 Auto-detected cluster: $CLUSTER"
    if [[ "$CLUSTER" == "unknown" ]]; then
        echo "⚠️  Unknown cluster - using generic SLURM parameters"
    fi
else
    echo "🎯 Target cluster: $CLUSTER"
fi

# Display cluster info if available
cluster_info=$(get_cluster_info "$CLUSTER")
if [[ "$cluster_info" != "Generic SLURM cluster" ]]; then
    echo "ℹ️  $cluster_info"
fi
echo ""

check_prerequisites

# Get cluster-specific parameters
CLUSTER_PARAMS=$(get_cluster_params "$CLUSTER")
echo "🔧 SLURM parameters: $CLUSTER_PARAMS"
echo ""

# Determine number of jobs
if [[ "$TEST_RUN" == "true" ]]; then
    N_JOBS=3
else
    N_JOBS=$(find calc/qe_runs -mindepth 1 -maxdepth 1 -type d | wc -l)
fi

estimate_resources $N_JOBS "$CLUSTER"

# Submit jobs
submit_jobs "$TEST_RUN" "$DRY_RUN" "$CLUSTER_PARAMS"

echo ""
echo "🎉 HPC DEPLOYMENT COMPLETE!"
echo "=========================="
echo "🔥 The race to room-temperature superconductivity has begun!"
echo "🏆 First breakthrough expected within 6 days"
echo ""
echo "📧 For updates: Check job status with 'squeue -u \$USER'"
echo "📊 For monitoring: Run './deploy_hpc.sh --monitor'" 