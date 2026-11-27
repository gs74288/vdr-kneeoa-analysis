# -*- coding: utf-8 -*-
"""Normal-only correlation for the calcitriol-response score."""
import sys, io, os, re
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\bulk"
CORE = ['GPX3', 'IGFBP3', 'MGP', 'PLPP3', 'PTGIS', 'TXNRD1']


def load(fn):
    d = pd.read_csv(os.path.join(D, fn), sep='\t', index_col=0)
    keep = [c for c in d.columns if not re.search(r'^(average|max)', str(c).strip(), re.I)]
    return d.loc[:, keep].groupby(level=0).mean()


oa, no = load('GSE114007_OA_counts.txt.gz'), load('GSE114007_normal_counts.txt.gz')
X = pd.concat([oa, no], axis=1)
vdr = X.loc['VDR'].astype(float)
present = [g for g in CORE if g in X.index]
Z = X.loc[present].astype(float)
Z = Z.sub(Z.mean(1), axis=0).div(Z.std(1).replace(0, np.nan), axis=0)
score = Z.mean(0)

for label, cols in (('all', list(X.columns)), ('OA', list(oa.columns)), ('normal', list(no.columns))):
    r = stats.spearmanr(vdr[cols], score[cols])
    print('%-8s n=%2d  rho=%+.3f  P=%.4g' % (label, len(cols), r.statistic, r.pvalue))
