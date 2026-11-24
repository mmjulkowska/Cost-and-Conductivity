"""EMMAX phenotype files for the Krs traits, in the format of the Pimp_Pareto_GWAS_for_Ruth pipeline.

Accession link: the label prefix in the tracings (e.g. '194' in 194_1_C) is the panel accession_code; the panel
geometry table gives its TGRC code, which is the FID/IID in ALLCHR.filt.prune.tfam order (LB/LC re-sequencing
replicates stay NA, as in v1_prep_phenotypes.R). Values: experiment-adjusted accession means from
krs_heritability.py; INT = Blom rank-based inverse normal over the 265-long vector, as in v1.
Writes outputs/krs/krs_gwas/ (copy this folder to the server): run_krs_gwas.sh, pheno/emmax/trait.<name>.INT.pheno, krs_gwas/pheno/raw/, analysis_set.txt, alignment_report.txt
"""
import pandas as pd, numpy as np, os, sys
from scipy.stats import norm
TRAITS = (sys.argv[1] if len(sys.argv) > 1 else 'logKrs.C,logKrs.S,logKrs.SminusC,KrsResid.C,KrsResid.S').split(',')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import krs_paths as KP
key = pd.read_csv(KP.inp('sample_key.csv'), dtype=str)
tfam = pd.read_csv(KP.inp('ALLCHR.filt.prune.tfam'), sep=r'\s+', header=None, dtype=str).iloc[:, :2]; tfam.columns = ['FID', 'IID']
assert (key.IID.values == tfam.IID.values).all()
KP.to_out()
A = pd.read_csv('krs_accession_values.csv', dtype={'accession_code': str}).set_index('accession_code')
A.index = A.index.str.zfill(3)
out = 'krs_gwas/pheno'; os.makedirs(f'{out}/emmax', exist_ok=True); os.makedirs(f'{out}/raw', exist_ok=True)


def INT(x):
    x = x.copy(); ok = x.notna(); r = x[ok].rank(method='average')
    x[ok] = norm.ppf((r - 3 / 8) / (ok.sum() + 1 / 4)); return x


# accession code -> TGRC code from the panel geometry table (240 accessions), then to the .tfam IID.
# (sample_key.genotype_code only covers the 199 accessions in Arjun's trait file, so it is NOT used.)
geo = pd.read_csv(KP.inp('01_panel_root_geometry.csv'), dtype=str)[['accession_code', 'tgrc_code']].drop_duplicates()
geo['code'] = geo.accession_code.str.zfill(3); assert geo.code.is_unique
la2code = dict(zip(geo.tgrc_code, geo.code))
code = key.IID.map(la2code).where(key.is_replicate != 'True')
rep = []
anyp = pd.Series(False, index=key.index)
for t in TRAITS:
    v = code.map(A[t]); anyp |= v.notna()
    pd.DataFrame({'FID': key.FID, 'IID': key.IID, 'v': v}).to_csv(f'{out}/raw/trait.{t}.RAW.pheno', sep='\t', header=False, index=False, na_rep='NA')
    pd.DataFrame({'FID': key.FID, 'IID': key.IID, 'v': INT(v)}).to_csv(f'{out}/emmax/trait.{t}.INT.pheno', sep='\t', header=False, index=False, na_rep='NA', float_format='%.10g')
    rep.append(f'{t:28s} accessions with a value: {A[t].notna().sum():4d}   in the genotype panel: {v.notna().sum():4d}')
key.loc[anyp, ['FID', 'IID']].to_csv(f'{out}/analysis_set.txt', sep='\t', header=False, index=False)
miss = sorted(set(A.index) - set(code.dropna()))
rep += ['', f'samples in .tfam: {len(key)}  (LB/LC replicates kept NA: {(key.is_replicate == "True").sum()})',
        f'accessions in the analysis set (any trait): {anyp.sum()}',
        f'traced accession codes not in the genotype panel ({len(miss)}): {" ".join(miss)}']
import shutil; shutil.copy(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_krs_gwas.sh'), 'krs_gwas/')
open(f'{out}/alignment_report.txt', 'w').write('\n'.join(rep) + '\n'); print('\n'.join(rep))
