# -*- coding: utf-8 -*-
"""Follow-up: CI for the within-RepC contrast, pooled RegC+RegC2 sensitivity, VDR share percentile."""
import sys, io, os, re
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

per = pd.read_csv(os.path.join(T, 'scdrs_vdr_within_gse324993.csv'))
print(per.round(3).to_string(index=False))
d = per.delta.values
ci = stats.t.interval(0.95, len(d) - 1, loc=d.mean(), scale=stats.sem(d))
print('\nRepC delta mean %+.3f  95%% CI %+.3f to %+.3f  n=%d' % (d.mean(), ci[0], ci[1], len(d)))

# share percentiles
for f, st in (('scdrs_contribution_gse324993.csv', 'RepC'),
              ('scdrs_contribution_gse104782.csv', 'RegC')):
    df = pd.read_csv(os.path.join(T, f))
    v = df[df.gene == 'VDR']
    if len(v):
        sh = v.share_pct.values[0]
        pct = 100.0 * (df.share_pct < sh).mean()
        z = v.z_expr.values[0]
        frac_small = 100.0 * (df.z_expr.abs() < abs(z)).mean()
        print('\n%s: VDR share %.3f%% (percentile %.1f of %d) | z %+.3f | %.0f%% of genes have smaller |z|'
              % (st, sh, pct, len(df), z, frac_small))
        top = df.nlargest(5, 'contribution')[['gene', 'z_expr', 'share_pct']]
        print('top positive contributors:')
        print(top.round(3).to_string(index=False))

# pooled RegC + RegC2 sensitivity in GSE324993
a = ad.read_h5ad(os.path.join(D, 'GSE324993.h5ad'), backed='r')
genes = list(a.var_names)
states = np.array([re.sub(r'\(.*\)', '', str(s)).strip() for s in a.obs['cell_type'].values])
donors = np.array([s.split('_')[0] for s in a.obs_names])
scores = pd.read_csv(os.path.join(D, 'scdrs_results_cell_scores.csv')) \
    .set_index('cell').loc[a.obs_names].norm_score.values
vi = genes.index('VDR')
vdr = np.zeros(len(states), dtype=bool)
for start in range(0, len(genes), 3000):
    if start <= vi < min(start + 3000, len(genes)):
        Xc = a[:, start:min(start + 3000, len(genes))].to_memory().X
        Xc = Xc.toarray() if hasattr(Xc, 'toarray') else np.asarray(Xc)
        vdr = Xc[:, vi - start] > 0
        break

for label, sel in (('RegC', states == 'RegC'),
                   ('RegC+RegC2', np.isin(states, ['RegC', 'RegC2']))):
    rows = []
    for dn in sorted(set(donors[sel])):
        m = sel & (donors == dn)
        p, n = scores[m & vdr], scores[m & ~vdr]
        if len(p) >= 5 and len(n) >= 5:
            rows.append((dn, len(p), len(n), p.mean(), n.mean(), p.mean() - n.mean()))
    r = pd.DataFrame(rows, columns=['donor', 'pos', 'neg', 'mp', 'mn', 'delta'])
    w = stats.wilcoxon(r.mp, r.mn)
    c = stats.t.interval(0.95, len(r) - 1, loc=r.delta.mean(), scale=stats.sem(r.delta))
    print('\n%s: %d donors | positive %d | delta %+.3f (%.3f to %.3f) | Wilcoxon P=%.3f'
          % (label, len(r), int((r.delta > 0).sum()), r.delta.mean(), c[0], c[1], w.pvalue))
    print(r.round(3).to_string(index=False))

print('\nCOL2A1 in GSE324993 var:', 'COL2A1' in set(genes))
