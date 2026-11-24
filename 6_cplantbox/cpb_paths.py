"""Locations for step 6 (CPlantBox root architectures). Builds on the repository's paths.py and step 5's krs_paths.py.

Inputs (not in git) are the step-5 inputs in data/krs/ ($KRS_DATA): 01_panel_root_geometry.csv, pooled_plants_by_day.csv,
pooled_laterals.csv, timeseries_roots.csv (made by 5_krs_panel/prep/), plus the PimpiBig2 day-5 SmartRoot RSML files
($PIMPIBIG2_RSML, default ~/mnt/Pimp_RSA/SmartRooted/pimpiBIG2; Dropbox Tomato_Pimpinellifolium/Pimp_RSA/SmartRooted/pimpiBIG2).
Segment conductivities come from step 2 (outputs/mecha/). Every script changes into CPB_OUT = outputs/cplantbox/,
so intermediate tables, parameters, simulations and figures are written and read there under their own names.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
for _p in (ROOT, os.path.join(ROOT, '5_krs_panel'), HERE):
    if _p not in sys.path: sys.path.insert(0, _p)
import paths, krs_paths as KP

CPB_OUT = os.environ.get('CPB_OUT', os.path.join(paths.OUT, 'cplantbox'))
os.makedirs(CPB_OUT, exist_ok=True)
RSML_DIR = os.path.expanduser(os.environ.get('PIMPIBIG2_RSML', '~/mnt/Pimp_RSA/SmartRooted/pimpiBIG2'))
RSML_SELECTED = os.path.join(HERE, 'inputs', 'rsml_selected.csv')   # the 72 PimpiBig2 day-5 tracings used (B 14, C 28, dT 30)
PARAMS_REF = os.path.join(HERE, 'params')                            # fitted parameters of the November 2026 runs
COND, ANAT_SEG, BRACKET = KP.COND, KP.ANAT_SEG, KP.BRACKET


def inp(name):
    return KP.inp(name)


def rsml(name):
    p = os.path.join(RSML_DIR, name)
    if not os.path.exists(p):
        sys.exit(f'Missing RSML {p}. Set PIMPIBIG2_RSML to the PimpiBig2 SmartRoot folder (see 6_cplantbox/README.md).')
    return p


os.chdir(CPB_OUT)
