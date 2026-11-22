# -*- coding: utf-8 -*-
"""scDRS sensitivity with absolute-value weights (standard practice) vs signed weights."""
import sys, io, os, re
import numpy as np
import pandas as pd
import anndata as ad
import scdrs
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"
H5 = os.path.join(D, 'GSE324993.h5ad')

z = pd.read_csv(os.path.join(D, 'GCST007090_gene_z.tsv'), sep='\t').dropna()
z = z[~z.GENE.duplicated()]
print('genes:', len(z), '| negative z:', int((z.Z < 0).sum()), flush=True)

obs = ad.read_h5ad(H5, backed='r').obs
states = np.array([re.sub(r'\(.*\)', '', str(s)).strip() for s in obs['cell_type'].values])
donors = np.array([s.split('_')[0] for s in obs.index])
print('cells:', len(states), flush=True)

rows = []
for tag, w in (('abs_z', np.abs(z.Z.values.astype(float))),
               ('signed_z_excl_VDR', z.Z.values.astype(float))):
    G = z.GENE.values
    if tag == 'signed_z_excl_VDR':
        keep = G != 'VDR'
        G, w2 = G[keep], w[keep]
    else:
        w2 = w
    print('\n=== %s: %d genes ===' % (tag, len(G)), flush=True)
    b = ad.read_h5ad(H5)
    scdrs.preprocess(b, n_mean_bin=20, n_var_bin=20, copy=False)
    sc = scdrs.score_cell(b, G, gene_weight=w2, n_ctrl=1000)
    s = sc['norm_score'].values
    for st in sorted(set(states)):
        rows.append(dict(run=tag, state=st, n=int((states == st).sum()),
                         mean=float(s[states == st].mean())))
    for focus in ('RegC', 'RepC'):
        rr = []
        for dn in sorted(set(donors)):
            m = donors == dn
            f, o = m & (states == focus), m & (states != focus)
            if f.sum() >= 5 and o.sum() >= 5:
                rr.append((dn, s[f].mean(), s[o].mean()))
        if len(rr) >= 3:
            r = pd.DataFrame(rr, columns=['donor', 'focus', 'other'])
            r['delta'] = r.focus - r.other
            wt = stats.wilcoxon(r.focus, r.other)
            print('  %s: %d/%d positive | delta %+.3f | P=%.3g'
                  % (focus, int((r.delta > 0).sum()), len(r), r.delta.mean(), wt.pvalue), flush=True)
    del b

df = pd.DataFrame(rows)
piv = df.pivot(index='state', columns='run', values='mean')
print('\n=== state means ===')
print(piv.round(3).to_string())
piv.to_csv(os.path.join(T, 'scdrs_weight_sensitivity.csv'))
print('\nwritten scdrs_weight_sensitivity.csv')
print('\nrank comparison:')
print(piv.rank(ascending=False).astype(int).to_string())
