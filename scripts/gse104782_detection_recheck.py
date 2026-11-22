# -*- coding: utf-8 -*-
"""Recompute GSE104782 VDR detection by state and donor from the raw count matrix."""
import sys, io, os
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

print('--- fig3_vdr_by_state.csv ---')
print(pd.read_csv(os.path.join(T, 'fig3_vdr_by_state.csv')).to_string(index=False))

cnt = pd.read_csv(os.path.join(D, 'GSE104782', 'GSE104782_allcells_UMI_count.txt.gz'),
                  sep='\t', index_col=0)
if cnt.shape[0] > cnt.shape[1]:
    cnt = cnt.T
cnt = cnt.groupby(level=0).sum()
print('\ncounts:', cnt.shape, '| VDR present:', 'VDR' in cnt.columns)
vdr = cnt['VDR']

lab = pd.read_excel(os.path.join(D, 'GSE104782',
                                 'GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx'))
lab = lab[['Cell', 'Cluster']].dropna()
d = pd.DataFrame({'cell': vdr.index, 'det': (vdr.values > 0).astype(int)}).merge(
    lab, left_on='cell', right_on='Cell', how='inner')
d['donor'] = d.cell.str.split('_').str[0]
print('cells:', len(d), '| detected:', int(d.det.sum()),
      '| %.2f%%' % (100 * d.det.mean()))

g = d.groupby('Cluster').det.agg(['size', 'sum', 'mean']).sort_values('mean', ascending=False)
g['pct'] = (100 * g['mean']).round(1)
print('\n=== detection by state ===')
print(g.to_string())

d['grp'] = np.where(d.Cluster.isin(['RegC', 'HomC']), 'RegHom', 'other')
per = d.groupby(['donor', 'grp']).det.agg(['size', 'sum', 'mean']).reset_index()
w = per.pivot(index='donor', columns='grp', values='mean')
n = per.pivot(index='donor', columns='grp', values='size')
out = pd.DataFrame({'RegHom': w.RegHom, 'other': w.other, 'n_RegHom': n.RegHom,
                    'delta_pp': (w.RegHom - w.other) * 100}).dropna()
print('\n=== per donor, RegC/HomC vs other ===')
print(out.round(4).to_string())

w_test = stats.wilcoxon(out.RegHom, out.other)
print('\npositive donors: %d/%d' % (int((out.delta_pp > 0).sum()), len(out)))
print('mean delta   %.2f pp | median delta %.2f pp' % (out.delta_pp.mean(), out.delta_pp.median()))
print('paired Wilcoxon P = %.4f' % w_test.pvalue)
ci = stats.t.interval(0.95, len(out) - 1, loc=out.delta_pp.mean(), scale=stats.sem(out.delta_pp))
print('mean 95%% CI %.2f to %.2f pp' % ci)

out.to_csv(os.path.join(T, 'gse104782_detection_per_donor.csv'))
print('\nwritten gse104782_detection_per_donor.csv')
