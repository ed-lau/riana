#!/usr/bin/env bash
# Build the standalone Riana GUI bundle with PyInstaller.
#
# Prereqs (in the active environment):
#     pip install -e ".[gui,packaging]"
# Usage (from anywhere):
#     packaging/build_app.sh
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root, regardless of caller CWD

# CFBundleShortVersionString must be numeric X.Y.Z, so strip any PEP 440
# suffix (e.g. "1.2.0.dev0" -> "1.2.0"); the .app still reports the exact
# build via the provenance header its outputs carry.
RIANA_VERSION="$(python -c 'import re, riana; m = re.match(r"\d+\.\d+\.\d+", riana.__version__); print(m.group(0) if m else "0.0.0")')"
export RIANA_VERSION
echo "Building Riana ${RIANA_VERSION} bundle with PyInstaller $(python -c 'import PyInstaller; print(PyInstaller.__version__)')..."

python -m PyInstaller packaging/riana-gui.spec --noconfirm --clean

echo
echo "Built:"
echo "  dist/riana-gui/         one-dir bundle (run dist/riana-gui/riana-gui)"
if [[ "$(uname)" == "Darwin" ]]; then
  echo "  dist/Riana.app          double-clickable macOS app"
fi
echo
echo "Smoke test (no display needed):"
echo "  dist/riana-gui/riana-gui --version"
