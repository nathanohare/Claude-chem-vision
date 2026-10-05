#!/usr/bin/env bash
# Installs the chem-vision plugin's dependencies.
#   ./setup.sh            core: RDKit, PyMuPDF, Pillow, OSRA
#   ./setup.sh --ml       also MolScribe and DECIMER (large; downloads model weights on first use)
# Each package is installed separately, so one blocked or failed package does
# not stop the others. The final `engines` report says what is available.
cd "$(dirname "$0")"
failed=()
for pkg in "pillow>=10" "rdkit>=2024.3" "pymupdf>=1.24"; do
  if ! python3 -m pip install "$pkg" >/tmp/chemvision_pip.log 2>&1 && \
     ! python3 -m pip install --break-system-packages "$pkg" >>/tmp/chemvision_pip.log 2>&1; then
    failed+=("$pkg")
    if grep -qiE "No matching distribution|403|Forbidden" /tmp/chemvision_pip.log; then
      echo "chem-vision setup: $pkg could not be downloaded. The package index appears to be blocked by this"
      echo "  environment's network allowlist (retrying with other commands will not help)."
    else
      echo "chem-vision setup: $pkg failed to install (see /tmp/chemvision_pip.log)."
    fi
  fi
done
if ! command -v osra >/dev/null; then
  if command -v apt-get >/dev/null; then sudo apt-get install -y osra >/dev/null 2>&1 || true
  elif command -v conda >/dev/null; then conda install -y -c conda-forge osra || true
  elif command -v brew >/dev/null; then echo "OSRA: install via conda (conda install -c conda-forge osra) or build from source"
  fi
fi
if [ "$1" = "--ml" ]; then
  python3 -m pip install molscribe huggingface_hub decimer || failed+=("ml extras")
fi
if [ ${#failed[@]} -gt 0 ]; then
  echo "chem-vision setup: not installed: ${failed[*]}. The plugin still runs; see 'rdkit_missing' /"
  echo "  'pymupdf_note' in the report below for what is unavailable and the fallbacks."
fi
python3 skills/chem-vision/scripts/chemvision.py engines
