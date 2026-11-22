# -*- coding: utf-8 -*-
"""scDRS specificity control: project OA, FinnGen arthrosis and FinnGen T2D onto GSE324993."""
import sys, io, os, re, gzip
import numpy as np
import pandas as pd
import anndata as ad
import scdrs
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"
H5 = os.path.join(D, 'GSE324993.h5ad')

GENE_SETS = {
    'OA (GCST007090)': os.path.join(D, 'GCST007090_gene_z.tsv'),
    'FinnGen arthrosis': os.path.join(D, 'finngen_M13_ARTHROSIS_gene_z.tsv'),
    'FinnGen T2D': os.path.join(D, 'finngen_T2D_gene_z.tsv'),
}

print('reading obs...', flush=True)
_obs = ad.read_h5ad(H5, backed='r').obs
print('cells:', len(_obs), flush=True)
states = np.array([re.sub(r'\(.*\)', '', str(s)).strip() for s in _obs['cell_type'].values])
donors = np.array([s.split('_')[0] for s in _obs.index])
print('states:', {s: int((states == s).sum()) for s in sorted(set(states))}, flush=True)

summary = []
for label, path in GENE_SETS.items():
    if not os.path.exists(path):
        print('missing', path, flush=True)
        continue
    z = pd.read_csv(path, sep='\t').dropna()
    z = z[~z.GENE.duplicated()]
    print('\n=== %s: %d genes ===' % (label, len(z)), flush=True)

    b = ad.read_h5ad(H5)
    scdrs.preprocess(b, n_mean_bin=20, n_var_bin=20, copy=False)
    sc = scdrs.score_cell(b, z.GENE.values, gene_weight=z.Z.values.astype(float), n_ctrl=1000)
    print('scored', sc.shape, flush=True)
    s = sc['norm_score'].values

    for st in sorted(set(states)):
        summary.append(dict(trait=label, state=st, n=int((states == st).sum()),
                            mean=float(s[states == st].mean()),
                            median=float(np.median(s[states == st]))))

    # donor-level contrast for the two compartments of interest
    for focus in ('RegC', 'RepC'):
        rows = []
        for dn in sorted(set(donors)):
            m = (donors == dn)
            f = m & (states == focus)
            o = m & (states != focus)
            if f.sum() >= 5 and o.sum() >= 5:
                rows.append((dn, s[f].mean(), s[o].mean()))
        if len(rows) >= 3:
            r = pd.DataFrame(rows, columns=['donor', 'focus', 'other'])
            r['delta'] = r.focus - r.other
            w = stats.wilcoxon(r.focus, r.other)
            print('  %s vs others: %d/%d positive | mean delta %+.3f | Wilcoxon P=%.3g'
                  % (focus, int((r.delta > 0).sum()), len(r), r.delta.mean(), w.pvalue), flush=True)

    del b

df = pd.DataFrame(summary)
piv = df.pivot(index='state', columns='trait', values='mean')
print('\n=== state means by trait ===')
print(piv.round(3).to_string())
piv.to_csv(os.path.join(T, 'scdrs_specificity_control.csv'))
print('\nwritten scdrs_specificity_control.csv')

print('\n=== rank of each state per trait ===')
print(piv.rank(ascending=False).astype(int).to_string())
