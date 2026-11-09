"""Where everything lives. Every script imports this file, so data and output locations are set in one place.

Inputs that are NOT in this repository (download them from Box, see data/README.md):
  DATA/AllDatammAndpixelsNov2026_2.xlsx   Emily Miller's measurement sheet
  DATA/Cross Sections/                     Emily Miller's cross-section photographs (only needed for the figure icons)
Set COSTCOND_DATA to use another folder, COSTCOND_OUT to write outputs elsewhere, ROOTSLICE_HOME for the RootSlice build.
"""
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('COSTCOND_DATA', os.path.join(ROOT, 'data'))
OUT = os.environ.get('COSTCOND_OUT', os.path.join(ROOT, 'outputs'))
ROOTSLICE = os.path.expanduser(os.environ.get('ROOTSLICE_HOME', '~/root_models/RootSlice'))

# ---- raw inputs (Box)
SHEET = os.path.join(DATA, 'AllDatammAndpixelsNov2026_2.xlsx')
IMAGES = os.path.join(DATA, 'Cross Sections')

# ---- model inputs kept in this repository
REPLICATES = os.path.join(ROOT, '1_granar', 'v5_replicates.csv')
GRANAR_PARAMS = os.path.join(ROOT, '1_granar', 'granar_params')
SEEDS = os.path.join(ROOT, '1_granar', 'seed_overrides.txt')
RS_TEMPLATE = os.path.join(ROOT, '3_rootslice', 'template_rootslice.xml')
RS_CONFIGS = os.path.join(ROOT, '3_rootslice', 'configs')
RS_CONFIGS_SS = os.path.join(ROOT, '3_rootslice', 'configs_smallstele')
RS_VALUES = os.path.join(ROOT, '3_rootslice', 'configs_values.csv')
ICON_SECTIONS = os.path.join(ROOT, 'figures', 'icon_sections.json')

# ---- outputs (not tracked; archived on Box)
ANAT = os.path.join(OUT, 'granar', 'anatomies')
GRANAR_SUMMARY = os.path.join(OUT, 'granar', 'granar_anatomy_summary.csv')
MECHA = os.path.join(OUT, 'mecha')
MECHA_RAW = os.path.join(MECHA, 'mecha_raw.csv')
COND = os.path.join(MECHA, 'conductivities_per_anatomy.csv')
STATS = os.path.join(MECHA, 'stats')
MC = os.path.join(MECHA, 'sensitivity_mc')
CHECKS = os.path.join(MECHA, 'checks')
RS_OUT = os.path.join(OUT, 'rootslice')
RS_VOL = os.path.join(RS_OUT, 'volume_analysis_configs.csv')
RS_VOL_SS = os.path.join(RS_OUT, 'volume_analysis_smallstele.csv')
COSTS = os.path.join(OUT, 'costs')
FIG = os.path.join(OUT, 'figures')
ICONS = os.path.join(FIG, 'icons')
DELIV = os.path.join(OUT, 'deliverables')

for _d in (os.path.dirname(GRANAR_SUMMARY), ANAT, MECHA, STATS, MC, CHECKS, RS_OUT, COSTS, FIG, ICONS, DELIV):
    os.makedirs(_d, exist_ok=True)
