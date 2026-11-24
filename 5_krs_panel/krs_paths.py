"""Locations for step 5 (Krs of the traced panel). Builds on the repository's paths.py.

Inputs (not in git; see README.md, section Inputs) are read from KRS_IN = $KRS_DATA or data/krs/.
Segment conductivities come from step 2 outputs (paths.COND, outputs/mecha/anatomy_by_segment.csv,
outputs/mecha/checks/bracket_laterals_summary.csv). Outputs go to KRS_OUT = outputs/krs/.
Every script in this folder changes into KRS_OUT, so outputs are written there under their own names.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import paths

KRS_IN = os.environ.get('KRS_DATA', os.path.join(paths.DATA, 'krs'))
KRS_OUT = os.path.join(paths.OUT, 'krs')
os.makedirs(KRS_OUT, exist_ok=True)

COND = paths.COND
ANAT_SEG = os.path.join(paths.MECHA, 'anatomy_by_segment.csv')
BRACKET = os.path.join(paths.CHECKS, 'bracket_laterals_summary.csv')


def inp(name):
    p = os.path.join(KRS_IN, name)
    if not os.path.exists(p):
        sys.exit(f'Missing input {p}. See 5_krs_panel/README.md, section Inputs.')
    return p


def to_out():
    os.chdir(KRS_OUT)
