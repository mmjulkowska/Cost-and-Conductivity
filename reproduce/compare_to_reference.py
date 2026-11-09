"""Compare every table in <outputs>/ with the frozen reference (download it from Box, see data/README.md).
Reference folder: $COSTCOND_REF (default <repo>/reference). Writes <outputs>/comparison_report.md.
Numbers must agree to a relative difference of 1e-6; text columns must be identical."""
import os, sys, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from paths import ROOT, OUT
REF = os.environ.get('COSTCOND_REF', os.path.join(ROOT, 'reference'))
if not os.path.isdir(REF): sys.exit(f'Reference folder not found: {REF}. Download it from Box and set COSTCOND_REF.')
TOL = 1e-6; rows = []
for ref in sorted(glob.glob(os.path.join(REF, '**', '*.csv'), recursive=True)):
    rel = os.path.relpath(ref, REF); new = os.path.join(OUT, rel)
    if not os.path.exists(new): rows.append((rel, 'missing in outputs', '')); continue
    a, b = pd.read_csv(new), pd.read_csv(ref)
    if a.shape != b.shape or list(a.columns) != list(b.columns): rows.append((rel, 'DIFFERENT SHAPE', f'{a.shape} vs {b.shape}')); continue
    worst, txt = 0.0, []
    for c in b.columns:
        if pd.api.types.is_numeric_dtype(b[c]) and pd.api.types.is_numeric_dtype(a[c]):
            x, y = a[c].to_numpy(float), b[c].to_numpy(float)
            if not np.array_equal(np.isnan(x), np.isnan(y)): worst = np.inf; continue
            m = ~np.isnan(y); d = np.abs(x[m] - y[m]) / np.maximum(np.abs(y[m]), 1e-30)
            if len(d): worst = max(worst, float(d.max()))
        elif not (a[c].astype(str) == b[c].astype(str)).all(): txt.append(c)
    ok = worst <= TOL and not txt
    rows.append((rel, 'identical' if ok else 'DIFFERENT', f'max rel diff {worst:.1e}' + (f'; text differs in {txt}' if txt else '')))
R = pd.DataFrame(rows, columns=['file', 'result', 'detail'])
n_bad = (~R.result.eq('identical')).sum()
with open(os.path.join(OUT, 'comparison_report.md'), 'w') as f:
    f.write(f'# Comparison with reference\n\nReference: `{REF}`\n\n{len(R)} tables, {len(R) - n_bad} identical (relative tolerance {TOL}), {n_bad} not.\n\n')
    f.write('| file | result | detail |\n|---|---|---|\n' + ''.join(f'| {r.file} | {r.result} | {r.detail} |\n' for r in R.itertuples()))
print(R.result.value_counts().to_string()); print('report:', os.path.join(OUT, 'comparison_report.md'))
sys.exit(1 if n_bad else 0)
