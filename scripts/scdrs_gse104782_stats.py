# -*- coding: utf-8 -*-
"""Cell-level and per-state donor comparisons for GSE104782 scDRS."""
import sys, io, os
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"
f = os.path.join(T, 'scdrs_GSE104782_cell_scores.csv')
print('exists:', os.path.exists(f))
d = pd.read_csv(f)
print('cols:', list(d.columns), '| rows:', len(d))

lab = pd.read_excel(r"C:\Users\liuha\OpenClawVDR\data\GSE104782\GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx")
lab = lab[['Cell', 'Cluster']].dropna()
d = d.merge(lab, left_on=d.columns[0], right_on='Cell', how='inner')
d['donor'] = d.Cell.str.split('_').str[0]

col = 'norm_score'
a = d.loc[d.Cluster == 'RegC', col]
b = d.loc[d.Cluster != 'RegC', col]
u = stats.mannwhitneyu(a, b, alternative='greater')
print('\ncell-level RegC vs rest: n=%d vs %d | mean %.3f vs %.3f | MWU one-sided P = %.3g' % (
    len(a), len(b), a.mean(), b.mean(), u.pvalue))

print('\nper-state, cell level:')
for st, g in d.groupby('Cluster'):
    o = d.loc[d.Cluster != st, col]
    p = stats.mannwhitneyu(g[col], o, alternative='greater').pvalue
    print('  %-7s n=%4d mean %+.3f  vs rest P=%.3g' % (st, len(g), g[col].mean(), p))

print('\nper-state donor-level mean scores (top 3 states):')
wide = d.groupby(['donor', 'Cluster'])[col].mean().unstack()
print(wide[['RegC', 'ProC', 'FC', 'preHTC']].round(3).to_string())

print('\nnumber of donors where RegC > ProC:', int((wide.RegC > wide.ProC).sum()), '/', wide.RegC.notna().sum())

summ = d.groupby('Cluster')[col].agg(['count', 'mean', 'median', 'std']).sort_values('mean', ascending=False)
summ.to_csv(os.path.join(T, 'scDRS_GSE104782_by_state.csv'))
print('\nrewrote by_state table')
