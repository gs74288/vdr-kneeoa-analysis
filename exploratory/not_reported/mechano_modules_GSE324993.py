# -*- coding: utf-8 -*-
"""Mechanobiology modules across GSE324993: state enrichment, correlation with scDRS, VDR adjustment."""
import sys, io, os, re
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

MODULES = {
    'YAP_TAZ_targets': ['CCN2', 'CTGF', 'CCN1', 'CYR61', 'ANKRD1', 'AMOTL2', 'TEAD1', 'TEAD4'],
    'mechanochannels': ['PIEZO1', 'PIEZO2', 'TRPV4', 'TRPC1', 'TRPC6'],
    'integrin_FAK_ROCK': ['ITGA5', 'ITGB1', 'PTK2', 'ROCK1', 'ROCK2', 'ACTA2', 'TAGLN', 'MYL9'],
    'ECM_stiffness': ['LOX', 'LOXL2', 'FN1', 'COL1A1', 'COL3A1', 'TIMP1', 'SPP1'],
}

a = ad.read_h5ad(os.path.join(D, 'GSE324993.h5ad'))
print('loaded', a.shape, flush=True)
states = np.array([re.sub(r'\(.*\)', '', str(s)).strip() for s in a.obs['cell_type'].values])
donors = np.array([s.split('_')[0] for s in a.obs_names])

sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)
for k, g in MODULES.items():
    present = [x for x in g if x in a.var_names]
    print('%-18s %d/%d genes' % (k, len(present), len(g)), flush=True)
    sc.tl.score_genes(a, present, score_name='sc_' + k)

vdr = np.asarray((a[:, 'VDR'].X > 0).todense()).ravel() if hasattr(a[:, 'VDR'].X, 'todense') else np.asarray(a[:, 'VDR'].X > 0).ravel()
df = pd.DataFrame({'cell': a.obs_names, 'donor': donors, 'state': states, 'vdr': vdr})
for k in MODULES:
    df['sc_' + k] = a.obs['sc_' + k].values

scd = pd.read_csv(os.path.join(D, 'scdrs_results_cell_scores.csv')).set_index('cell')
df = df.merge(scd[['norm_score']], left_on='cell', right_index=True, how='left')
print('cells:', len(df), '| with scDRS:', int(df.norm_score.notna().sum()))

mech = ['sc_YAP_TAZ_targets', 'sc_mechanochannels', 'sc_integrin_FAK_ROCK', 'sc_ECM_stiffness']
df['mechano_mean'] = df[mech].mean(1)

print('\n=== mean module score by state (z across states) ===')
g = df.groupby('state')[mech + ['mechano_mean', 'norm_score']].mean()
z = (g - g.mean()) / g.std()
print(z.round(2).to_string())

print('\n=== donor-state level: correlation with scDRS ===')
ds = df.groupby(['donor', 'state'])[mech + ['mechano_mean', 'norm_score']].mean().dropna()
for m in mech + ['mechano_mean']:
    r = stats.spearmanr(ds[m], ds.norm_score)
    print('  %-24s rho=%+.3f P=%.4g (n=%d)' % (m, r.statistic, r.pvalue, len(ds)))

print('\n=== within RegC: does VDR association survive mechano adjustment? ===')
reg = df[df.state == 'RegC'].dropna(subset=['norm_score'])
rows = []
for d in sorted(reg.donor.unique()):
    s = reg[reg.donor == d]
    if s.vdr.sum() < 5:
        continue
    # simple model: score ~ vdr + mechano (within donor)
    X = np.column_stack([np.ones(len(s)), s.vdr.astype(float), s.mechano_mean])
    y = s.norm_score.values
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    rows.append((d, len(s), int(s.vdr.sum()), beta[1], beta[2]))
r = pd.DataFrame(rows, columns=['donor', 'n', 'n_vdr', 'beta_vdr', 'beta_mech'])
print(r.round(3).to_string(index=False))
print('\nmean beta_vdr %+.3f | median %+.3f | negative in %d/%d donors'
      % (r.beta_vdr.mean(), r.beta_vdr.median(), int((r.beta_vdr < 0).sum()), len(r)))
print('mean beta_mech %+.3f | positive in %d/%d'
      % (r.beta_mech.mean(), int((r.beta_mech > 0).sum()), len(r)))

df.to_csv(os.path.join(T, 'GSE324993_mechano_scores.csv'), index=False)
print('\nwritten GSE324993_mechano_scores.csv')
