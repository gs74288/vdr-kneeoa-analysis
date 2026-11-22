# -*- coding: utf-8 -*-
"""Non-circular test: are canonical VDR target genes expressed where VDR is detected?"""
import sys, io, os
import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\tools\scenic"
OUT = r"C:\Users\liuha\OpenClawVDR\repo\tables"

X = scipy.io.mmread(os.path.join(D, 'counts.mtx')).tocsr()
genes = [g for g in open(os.path.join(D, 'genes.txt')).read().split('\n') if g]
cells = [c for c in open(os.path.join(D, 'cells.txt')).read().split('\n') if c]
meta = pd.read_csv(os.path.join(D, 'cellmeta.csv')).set_index('cell').loc[cells]
gi = {g: i for i, g in enumerate(genes)}

TARGETS = ['CYP24A1', 'CYP27B1', 'TRPV6', 'SPP1', 'BGLAP', 'CD14', 'DEFB4A', 'CAMP',
           'LRP5', 'RUNX2', 'VDR']
present = [t for t in TARGETS if t in gi]
print('targets present:', present)

tot = np.asarray(X.sum(axis=0)).ravel()
tot[tot == 0] = 1
rows = {}
for t in present:
    v = np.asarray(X[gi[t], :].todense()).ravel()
    rows[t] = v * 1e4 / tot          # CP10K
df = pd.DataFrame(rows, index=cells)
df['cluster'] = meta.cluster.values
df['donor'] = meta.donor.values
df['vdr_det'] = (df['VDR'] > 0).astype(int)

print('\n=== detection rate by target gene ===')
for t in present:
    print('  %-9s %.2f%% of cells' % (t, 100 * (df[t] > 0).mean()))

print('\n=== CYP24A1 and CYP27B1 by cluster ===')
g = df.groupby('cluster').agg(
    n=('VDR', 'size'),
    vdr_pct=('vdr_det', lambda x: 100 * x.mean()),
    cyp24a1_pct=('CYP24A1', lambda x: 100 * (x > 0).mean()) if 'CYP24A1' in df else ('CYP24A1', 'size'),
    cyp24a1_mean=('CYP24A1', 'mean') if 'CYP24A1' in df else ('CYP24A1', 'size'),
)
print(g.sort_values('vdr_pct', ascending=False).round(2).to_string())

if 'CYP24A1' in df:
    print('\n=== CYP24A1 in VDR+ vs VDR- cells (within cluster) ===')
    for cl, sub in df.groupby('cluster'):
        a = sub.loc[sub.vdr_det == 1, 'CYP24A1']
        b = sub.loc[sub.vdr_det == 0, 'CYP24A1']
        if len(a) > 20 and len(b) > 20:
            t = stats.ttest_ind(a, b, equal_var=False)
            print('  %-4s n+=%4d n-=%5d  %.3f vs %.3f  P=%.2e'
                  % (cl, len(a), len(b), a.mean(), b.mean(), t.pvalue))
    r = stats.spearmanr(df['VDR'], df['CYP24A1'])
    print('\nSpearman(VDR, CYP24A1) across all cells: rho=%.3f  P=%.2e' % (r.statistic, r.pvalue))
