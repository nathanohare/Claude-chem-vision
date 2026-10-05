#!/usr/bin/env bash
# Installs the chem-vision plugin's dependencies.
#   ./setup.sh            core: RDKit, PyMuPDF, Pillow, OSRA
#   ./setup.sh --ml       also MolScribe and DECIMER (large; downloads model weights on first use)
set -e
cd "$(dirname "$0")"
python3 -m pip install -r requirements.txt
if ! command -v osra >/dev/null; then
  if command -v apt-get >/dev/null; then sudo apt-get install -y osra || true
  elif command -v brew >/dev/null; then echo "OSRA: install via conda (conda install -c conda-forge osra) or build from source"
  elif command -v conda >/dev/null; then conda install -y -c conda-forge osra || true
  fi
fi
if [ "$1" = "--ml" ]; then
  python3 -m pip install molscribe huggingface_hub decimer
fi
python3 skills/chem-vision/scripts/chemvision.py engines
