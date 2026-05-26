#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

conda env create -f environment.yml

cat <<'EOF'

Conda environment created.

Next:
  conda activate pymol-siteview
  pymol-siteview --help

Example:
  pymol-siteview your_structure.pdb -o siteview_out

EOF
