# -*- coding: utf-8 -*-
"""Depth control: does the VDR+ excess within RegC survive adjustment for detected-gene count?"""
import sys, io, os, re
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

a = ad.read_h5ad(os.path.join(D, 'GSE324993.h5ad'), backed='r')
genes = list(a.var_names)
n_cells, n_genes = a.shape
states = np.array([re.sub(r'\(.*\)', '', str(s)).strip() for s in a.obs['cell_type'].values])
donors = np.array([s.split('_')[0] for s in a.obs_names])
scores = pd.read_csv(os.path.join(D, 'scdrs_results_cell_scores.csv')) \
    .set_index('cell').loc[a.obs_names].norm_score.values

ndet = np.zeros(n_cells, dtype=np.int32)
vi = genes.index('VDR')
vdr = np.zeros(n_cells, dtype=bool)
CH = 1500
for start in range(0, n_genes, CH):
    end = min(start + CH, n_genes)
    Xc = a[:, start:end].to_memory().X
    Xc = Xc.toarray() if hasattr(Xc, 'toarray') else np.asarray(Xc)
    ndet += (Xc > 0).sum(1).astype(np.int32)
    if start <= vi < end:
        vdr = Xc[:, vi - start] > 0
    print('chunk %d-%d' % (start, end), flush=True)

df = pd.DataFrame({'donor': donors, 'state': states, 'score': scores,
                   'vdr': vdr, 'ndet': ndet})
df.to_csv(os.path.join(T, 'gse324993_cell_depth.csv'), index=False)
print('\nmedian detected genes by VDR status:')
print(df.groupby('vdr').ndet.median().to_string())
print('\nSpearman(Stringent) score vs ndet overall: rho=%.3f P=%.3g'
      % stats.spearmanr(df.score, df.ndet))

for label, sel in (('RegC', df.state == 'RegC'), ('RegC+RegC2', df.state.isin(['RegC', 'RegC2'])),
                   ('RepC', df.state == 'RepC')):
    s = df[sel].copy()
    s['lnd'] = np.log(s.ndet)
    # residualise score on log detected genes within this compartment
    s['res'] = s.score - np.polyval(np.polyfit(s.lnd, s.score, 1), s.lnd)
    rows = []
    for dn in sorted(s.donor.unique()):
        m = s.donor == dn
        p_adj, n_adj = s.res[m & s.vdr], s.res[m & ~s.vdr]
        p_raw, n_raw = s.score[m & s.vdr], s.score[m & ~s.vdr]
        if len(p_adj) >= 5 and len(n_adj) >= 5:
            rows.append((dn, len(p_adj), len(n_adj),
                         p_raw.mean() - n_raw.mean(), p_adj.mean() - n_adj.mean()))
    r = pd.DataFrame(rows, columns=['donor', 'pos', 'neg', 'delta_raw', 'delta_adj'])
    w_raw = stats.wilcoxon(r.delta_raw)
    w_adj = stats.wilcoxon(r.delta_adj)
    print('\n=== %s (n=%d cells; %d VDR+) ===' % (label, len(s), int(s.vdr.sum())))
    print(r.round(3).to_string(index=False))
    print('raw      : %d/%d positive | mean %+.3f | Wilcoxon P=%.4f'
          % (int((r.delta_raw > 0).sum()), len(r), r.delta_raw.mean(), w_raw.pvalue))
    print('adjusted : %d/%d positive | mean %+.3f | Wilcoxon P=%.4f'
          % (int((r.delta_adj > 0).sum()), len(r), r.delta_adj.mean(), w_adj.pvalue))
    ci = stats.t.interval(0.95, len(r) - 1, loc=r.delta_adj.mean(),
                          scale=stats.sem(r.delta_adj))
    print('adjusted 95%% CI %+.3f to %+.3f' % ci)
