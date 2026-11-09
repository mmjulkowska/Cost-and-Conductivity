#!/bin/bash
# Install MECHA (openalea.mecha, MECHA-dev "develop") + GRANAP and the Python analysis packages, conda env mecha_dev (macOS, Apple Silicon).
# Pinned to the exact commits used for the Nov 2026 re-run, plus a one-line fix for a merge bug in develop.
# Run from the repository:  bash env/install_mecha_mac.sh      (needs conda and git; installs into ~/root_models)
set -e
DEST="$HOME/root_models"; ENV=mecha_dev
MECHA_COMMIT=945b2be429ab0493346cc8d10339f116ee8de385   # develop, 2026-10-26 (includes barrier 10 + 11)
GRANAP_COMMIT=248adcafcebc040bde9321509b299ad2811e67cf  # GRANAP-dev develop at the same time
HERE="$(cd "$(dirname "$0")" && pwd)"
LOG="$HERE/install_mecha_log.txt"; exec > >(tee "$LOG") 2>&1
echo "== $(date)"
for c in /opt/anaconda3/bin/conda /opt/homebrew/Caskroom/miniforge/base/bin/conda "$HOME/miniforge3/bin/conda"; do
  [ -x "$c" ] && eval "$("$c" shell.bash hook)" && break; done
command -v conda >/dev/null || { echo "conda not found"; exit 1; }
mkdir -p "$DEST" && cd "$DEST"
[ -d MECHA-dev ]  || git clone https://github.com/openalea-incubator/MECHA-dev.git
[ -d GRANAP-dev ] || git clone https://github.com/openalea-incubator/GRANAP-dev.git
git -C MECHA-dev fetch -q origin && git -C MECHA-dev checkout -q -f "$MECHA_COMMIT"
git -C GRANAP-dev fetch -q origin && git -C GRANAP-dev checkout -q -f "$GRANAP_COMMIT"
# Fix: in this commit a helper `_contact_nodes` was accidentally merged under the name `_apply_soil_boundary`
python3 - "$DEST/MECHA-dev/src/openalea/mecha/utils/hydraulic_solver.py" <<'PY'
import sys; p=sys.argv[1]; s=open(p).read()
sig="    def _apply_soil_boundary(self, x_contact, height, thickness, kw, barrier, boundary, rhs_s, rhs_C):\n"
if "def _contact_nodes(" in s: print("patch already applied")
else:
    i=s.find(sig); j=s.find(sig,i+1); assert i>0 and j>i, "unexpected file layout"
    open(p,"w").write(s[:i]+"    def _contact_nodes(self, x_contact):\n"+s[i+len(sig):]); print("patched hydraulic_solver.py")
PY
# native arm64 env (no Rosetta needed, unlike classic MECHA)
conda env list | grep -q "^$ENV " || CONDA_SUBDIR=osx-arm64 conda create -y -n $ENV -c conda-forge --override-channels \
  python=3.12 numpy scipy networkx lxml geopandas pandas matplotlib shapely pydantic opencv munkres scikit-learn roifile pip git \
  statsmodels openpyxl scikit-image pillow
conda run -n $ENV python -m pip install -q --no-deps -e "$DEST/GRANAP-dev" -e "$DEST/MECHA-dev"
echo "== self-test (needs 1_granar outputs; expect MP1_cell1 with B6: Kx 2.7612e-02, kr 5.9217e-05)"
A="${COSTCOND_OUT:-$HERE/../outputs}/granar/anatomies/MP1_cell1_root.xml"
if [ -f "$A" ]; then
  cd "$HERE/../2_mecha" && conda run -n $ENV python -c "
from run_one import run
r=run('$A',[6],1000.)[0]
print('Kx %.4e  kr %.4e' % (r['Kx'], r['kr']))
print('SELF-TEST', 'OK' if abs(r['Kx']/2.7612e-02-1)<1e-3 and abs(r['kr']/5.9217e-05-1)<1e-3 else 'MISMATCH')"
else echo "   (skipped: build the GRANAR anatomies first, then re-run this script)"; fi
echo "== done. Use:  conda activate $ENV ; bash reproduce/run_all.sh"
