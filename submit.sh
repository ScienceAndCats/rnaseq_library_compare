#!/usr/bin/env bash
# Use from any directory: bash /path/to/repo/submit.sh /path/to/config.yaml
set -euo pipefail
REPO_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
CONFIG="${1:-$REPO_DIR/config.yaml}"
if [[ $# -gt 0 ]]; then shift; fi
if [[ "$CONFIG" != /* ]]; then CONFIG="$PWD/$CONFIG"; fi
if [[ ! -f "$CONFIG" ]]; then echo "ERROR: Config not found: $CONFIG" >&2; exit 1; fi
command -v sbatch >/dev/null || { echo 'ERROR: sbatch is not on PATH.' >&2; exit 1; }
command -v python >/dev/null || { echo 'ERROR: Activate rnaseq-library-compare first.' >&2; exit 1; }
export RNA_COMPARE_REPO_DIR="$REPO_DIR"
export RNA_COMPARE_PYTHON="$(python -c 'import sys; print(sys.executable)')"
# The environment and absolute interpreter are captured at submission, not found
# by guessing where Slurm copied the batch script.
exec sbatch "$REPO_DIR/run_comparison.sbatch" "$CONFIG" "$@"
