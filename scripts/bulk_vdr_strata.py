# -*- coding: utf-8 -*-
"""Bulk: do VDR-high OA samples show a stronger matrix-degradation programme?"""
import sys, io, os
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
BULK = r"C:\Users\liuha\OpenClawVDR\bulk"

oa = pd.read_csv(os.path.join(BULK, 'GSE114007_OA_counts.txt.gz'), sep='\t', index_col=0)
oa = oa.groupby(level=0).mean()
print('OA samples:', oa.shape[1])

TARGETS = ['MMP13', 'ADAMTS5', 'COL10A1', 'ACAN', 'COL2A1', 'MMP3', 'ADAMTS4',
           'COL1A1', 'SOX9', 'CYP24A1']
present = [t for t in TARGETS if t in oa.index]
print('genes present:', present)

vdr = oa.loc['VDR'].astype(float)
med = vdr.median()
grp = np.where(vdr > med, 'VDR-high', 'VDR-low')
print('\nVDR median %.3f | high n=%d low n=%d' % (med, (grp == 'VDR-high').sum(),
                                                  (grp == 'VDR-low').sum()))

rows = []
for g in present:
    a = oa.loc[g].astype(float)[grp == 'VDR-high']
    b = oa.loc[g].astype(float)[grp == 'VDR-low']
    t = stats.ttest_ind(a, b, equal_var=False)
    rows.append(dict(gene=g, high=a.mean(), low=b.mean(), diff=a.mean() - b.mean(),
                     p=t.pvalue))
res = pd.DataFrame(rows).sort_values('p')
print('\n=== VDR-high vs VDR-low OA samples ===')
print(res.round(3).to_string(index=False))

print('\n=== Spearman(VDR, gene) across 22 OA samples ===')
cor = []
for g in present:
    r = stats.spearmanr(vdr, oa.loc[g].astype(float))
    cor.append(dict(gene=g, rho=r.statistic, p=r.pvalue))
print(pd.DataFrame(cor).sort_values('p').round(3).to_string(index=False))
