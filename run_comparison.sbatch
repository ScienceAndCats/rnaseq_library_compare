#!/usr/bin/env bash
#SBATCH --job-name=rnaseq_compare
#SBATCH --account=gts-mwhiteley3-biocluster
#SBATCH --partition=cpu-small
#SBATCH --qos=inferno
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=rnaseq_compare_%j.log
#SBATCH --export=ALL

set -euo pipefail

# Slurm executes a spool COPY of this script. Do NOT locate Python beside
# BASH_SOURCE here. submit.sh exports the real repo; direct sbatch uses the
# directory from which the job was submitted.
REPO_DIR="${RNA_COMPARE_REPO_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
PYTHON="${RNA_COMPARE_PYTHON:-python}"
CONFIG="${1:-config.yaml}"
if [[ $# -gt 0 ]]; then shift; fi

if [[ ! -f "$REPO_DIR/compare_rnaseq.py" ]]; then
    echo "ERROR: Missing $REPO_DIR/compare_rnaseq.py" >&2
    echo 'Submit from the repository directory, or use /path/to/repo/submit.sh.' >&2
    exit 1
fi
# Relative config paths are relative to the submitter's directory. Paths INSIDE
# the YAML are independently resolved by Python relative to the YAML file.
if [[ "$CONFIG" != /* ]]; then
    CONFIG="${SLURM_SUBMIT_DIR:-$PWD}/$CONFIG"
fi
if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config not found: $CONFIG" >&2
    exit 1
fi
if ! command -v "$PYTHON" >/dev/null 2>&1; then
    echo 'ERROR: Python is unavailable. Activate rnaseq-library-compare before submitting.' >&2
    exit 1
fi

# Let PyDESeq2 use the requested processes; prevent each BLAS worker spawning
# another full complement of CPU threads.
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
export MPLCONFIGDIR="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/rnaseq_compare_mpl_${SLURM_JOB_ID:-$$}"
mkdir -p "$MPLCONFIGDIR"
trap 'rm -rf -- "$MPLCONFIGDIR"' EXIT

printf 'Repository: %s\nConfig: %s\nPython: %s\n' "$REPO_DIR" "$CONFIG" "$PYTHON"
"$PYTHON" -c 'import sys, numpy, pandas, scipy, matplotlib, yaml; print("Interpreter:", sys.executable)'
"$PYTHON" "$REPO_DIR/compare_rnaseq.py" "$CONFIG" \
    --threads "${SLURM_CPUS_PER_TASK:-4}" "$@"
