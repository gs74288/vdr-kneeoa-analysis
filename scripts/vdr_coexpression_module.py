# -*- coding: utf-8 -*-
"""VDR co-expression module in GSE152805 (explicitly NOT full SCENIC).

Rationale: VDR transcript detection is sparse, so we ask whether a VDR co-expression
module (the genes whose expression tracks VDR across cells) is enriched in regulatory
states and in VDR-detected cells.
"""
import sys, io, os
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\tools\scenic"
OUT = r"C:\Users\liuha\OpenClawVDR\repo\tables"
os.makedirs(OUT, exist_ok=True)

X = scipy.io.mmread(os.path.join(D, 'counts.mtx')).tocsr()
genes = [g for g in open(os.path.join(D, 'genes.txt')).read().split('\n') if g]
cells = [c for c in open(os.path.join(D, 'cells.txt')).read().split('\n') if c]
meta = pd.read_csv(os.path.join(D, 'cellmeta.csv')).set_index('cell').loc[cells]
print('matrix %d genes x %d cells' % X.shape, flush=True)

n_expr = np.asarray((X > 0).sum(axis=1)).ravel()
keep = n_expr >= 50
X = X[keep]
genes = [g for g, k in zip(genes, keep) if k]
print('genes kept:', len(genes), flush=True)

tot = np.asarray(X.sum(axis=0)).ravel()
tot[tot == 0] = 1
Xn = X.multiply(1e4 / tot).tocsc()
Xn.data = np.log1p(Xn.data)
M = Xn.T.tocsr()                       # cells x genes
print('normalised', flush=True)

gi = {g: i for i, g in enumerate(genes)}
vdr_i = gi['VDR']
vdr = np.asarray(M[:, vdr_i].todense()).ravel()
det = vdr > 0
print('VDR detected: %d / %d (%.2f%%)' % (det.sum(), len(det), 100 * det.mean()), flush=True)

# Spearman correlation of VDR with every gene
rank = lambda a: stats.rankdata(a)
rv = rank(vdr)
Dense = M.toarray().astype(np.float32)
print('dense:', Dense.shape, flush=True)
R = np.apply_along_axis(stats.rankdata, 0, Dense)
Rc = R - R.mean(0)
Rv = rv - rv.mean()
denom = np.sqrt((Rc ** 2).sum(0) * (Rv ** 2).sum())
rho = np.where(denom > 0, (Rc * Rv[:, None]).sum(0) / denom, 0.0)

res = pd.DataFrame({'gene': genes, 'rho': rho})
res = res[res.gene != 'VDR'].sort_values('rho', ascending=False)
res.to_csv(os.path.join(OUT, 'vdr_coexpression_ranking.csv'), index=False)
top = res.head(200).gene.tolist()
print('\ntop 15 co-expressed genes:', ', '.join(res.head(15).gene), flush=True)
print('top 200 span rho %.3f to %.3f' % (res.head(200).rho.max(), res.head(200).rho.min()))

# module score
cols = [gi[g] for g in top]
Z = (Dense[:, cols] - Dense[:, cols].mean(0)) / (Dense[:, cols].std(0) + 1e-9)
score = Z.mean(1)

meta['vdr_score'] = score
meta['vdr_det'] = det.astype(int)
meta.to_csv(os.path.join(OUT, 'vdr_module_scores.csv'))

print('\n=== module score by cluster ===')
print(meta.groupby('cluster').vdr_score.agg(['mean', 'size']).sort_values('mean', ascending=False).round(3).to_string())

print('\n=== VDR+ vs VDR- within cluster ===')
for cl, g in meta.groupby('cluster'):
    a = g.loc[g.vdr_det == 1, 'vdr_score']
    b = g.loc[g.vdr_det == 0, 'vdr_score']
    if len(a) > 20 and len(b) > 20:
        t = stats.ttest_ind(a, b, equal_var=False)
        print('  %-6s n+=%4d n-=%5d  mean %.3f vs %.3f  P=%.2e'
              % (cl, len(a), len(b), a.mean(), b.mean(), t.pvalue))
