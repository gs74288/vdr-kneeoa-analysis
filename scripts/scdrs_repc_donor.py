# -*- coding: utf-8 -*-
"""Donor-level scDRS: is the OA genetic signal concentrated in RepC?"""
import sys, io, gzip, re
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
SC = r"C:\Users\liuha\OpenClawVDR\data\scdrs_results_cell_scores.csv"
AN = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\data\GSE324993_cell_annotation.txt.gz"

rows, started = [], False
with gzip.open(AN, 'rt', errors='replace') as fh:
    for ln in fh:
        ln = ln.rstrip('\n')
        if ln.startswith('Cell_Label'):
            started = True
            continue
        if not started or not ln.strip():
            continue
        p = ln.split('\t')
        if len(p) < 2:
            continue
        rows.append((p[0].strip(), re.sub(r'\(.*\)', '', p[1].strip()).strip()))
ann = pd.DataFrame(rows, columns=['cell', 'state'])
ann['donor'] = ann.cell.str.split('_').str[0]

d = pd.read_csv(SC).merge(ann, on='cell', how='inner')
print('cells with scores and labels:', len(d))

per = d.groupby(['donor', 'state']).norm_score.agg(['mean', 'size']).reset_index()
wide = per.pivot(index='donor', columns='state', values='mean')
cnt = per.pivot(index='donor', columns='state', values='size')
print('\nmean normalised score by donor and state:')
print(wide.round(3).to_string())

rep = wide['RepC']
other = d[~d.state.isin(['RepC'])].groupby('donor').norm_score.mean()
print('\nRepC vs all other cells, by donor:')
out = pd.DataFrame({'RepC': rep, 'other': other, 'delta': rep - other,
                    'n_RepC': cnt['RepC']}).dropna()
print(out.round(3).to_string())

w = stats.wilcoxon(out.RepC, out.other)
t = stats.ttest_rel(out.RepC, out.other)
print('\npositive donors: %d/%d' % (int((out.delta > 0).sum()), len(out)))
print('mean delta %.3f (SD %.3f), median %.3f' % (out.delta.mean(), out.delta.std(ddof=1),
                                                  out.delta.median()))
print('paired Wilcoxon P = %.4f | paired t P = %.4f' % (w.pvalue, t.pvalue))
ci = stats.t.interval(0.95, len(out) - 1, loc=out.delta.mean(), scale=stats.sem(out.delta))
print('mean delta 95%% CI %.3f to %.3f' % ci)

out.to_csv(r"C:\Users\liuha\OpenClawVDR\repo\tables\scDRS_RepC_per_donor.csv")
print('\nwritten scDRS_RepC_per_donor.csv')
