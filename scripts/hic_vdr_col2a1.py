# -*- coding: utf-8 -*-
"""Assess VDR-COL2A1 contact strength versus a distance-matched control block."""
import sys, io, os
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data\chromatin"


def load(name, cols=('x', 'y', 'v')):
    p = os.path.join(D, name)
    df = pd.read_csv(p, sep='\t', header=None, names=list(cols))
    return df


tg = load('hic_vdr_col2a1_oe.txt')
tg.columns = ['x', 'y', 'oe']
obs = load('hic_vdr_col2a1_obs.txt')
obs.columns = ['x', 'y', 'n']
ct = load('hic_ctrl_oe.txt')
ct.columns = ['x', 'y', 'oe']

for name, d in (('target', tg), ('control', ct)):
    d['dist'] = (d.y - d.x).abs() + 2500
print('target bins:', len(tg), '| control bins:', len(ct))

# distance-decay comparison
bins = np.arange(0, 400001, 25000)
for name, d in (('VDR x COL2A1', tg), ('control', ct)):
    d['db'] = pd.cut(d.dist, bins)
    g = d.groupby('db', observed=True).oe.agg(['median', 'count'])
    print('\n%s: median O/E by distance' % name)
    print(g.head(10).round(3).to_string())

# specific: VDR gene (47.84-47.95) x COL2A1 gene (47.97-48.01)
sel = tg[((tg.x >= 47830000) & (tg.x <= 47960000) &
          (tg.y >= 47960000) & (tg.y <= 48030000))]
print('\n--- VDR region x COL2A1 region (all bins) ---')
print('n bins %d | median O/E %.3f | mean %.3f | max %.3f | %%>1.5: %.1f%%'
      % (len(sel), sel.oe.median(), sel.oe.mean(), sel.oe.max(),
         100 * (sel.oe > 1.5).mean()))

# distance-only comparison, 100-200 kb band
band_t = tg[(tg.dist >= 100000) & (tg.dist <= 200000)].oe
band_c = ct[(ct.dist >= 100000) & (ct.dist <= 200000)].oe
from scipy import stats
print('\n100-200 kb band: target median %.3f (n=%d) | control median %.3f (n=%d) | MWU P=%.3g'
      % (band_t.median(), len(band_t), band_c.median(), len(band_c),
         stats.mannwhitneyu(band_t, band_c).pvalue))

# the closest bin pair to (47845054, 48009982)
tg['d_target'] = np.abs(tg.x - 47845054) + np.abs(tg.y - 48009982)
best = tg.nsmallest(5, 'd_target')
print('\n--- bins nearest the rs7975232 x rs11168351 pair ---')
print(best[['x', 'y', 'oe', 'dist']].to_string(index=False))
m = obs.assign(dist=(obs.y - obs.x).abs())
print('\nobserved counts overall: median %.1f max %.0f' % (m.n.median(), m.n.max()))
