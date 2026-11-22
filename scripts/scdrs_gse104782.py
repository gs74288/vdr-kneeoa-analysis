# -*- coding: utf-8 -*-
"""Second-atlas scDRS: map the knee OA GWAS onto GSE104782 (1,600 cells, 10 donors)."""
import sys, io, os, gzip
import numpy as np
import pandas as pd
import anndata as ad
import scdrs
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
OUT = r"C:\Users\liuha\OpenClawVDR\repo\tables"

# ---------- expression ----------
print('loading counts...', flush=True)
cnt = pd.read_csv(os.path.join(D, 'GSE104782', 'GSE104782_allcells_UMI_count.txt.gz'),
                  sep='\t', index_col=0)
print('raw shape:', cnt.shape, flush=True)
if cnt.shape[0] > cnt.shape[1]:
    cnt = cnt.T
    print('transposed ->', cnt.shape, flush=True)
cnt = cnt.groupby(level=0).sum()
print('cells x genes:', cnt.shape, flush=True)

# ---------- labels ----------
lab = pd.read_excel(os.path.join(D, 'GSE104782',
                                 'GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx'))
lab = lab[['Cell', 'Cluster']].dropna()
lab['donor'] = lab.Cell.str.split('_').str[0]
print('labelled cells:', len(lab), '| states:', dict(lab.Cluster.value_counts()), flush=True)

common = [c for c in lab.Cell if c in cnt.index]
print('cells in both:', len(common), flush=True)
cnt = cnt.loc[common]
meta = lab.set_index('Cell').loc[common]

# ---------- gene set ----------
z = pd.read_csv(os.path.join(D, 'GCST007090_gene_z.tsv'), sep='\t').dropna()
z = z[~z.GENE.duplicated()]
gene_list = z.GENE.values
gene_weight = z.Z.values.astype(float)
print('gene set:', len(gene_list), flush=True)

# ---------- scDRS ----------
X = cnt.values.astype(np.float64)
print('raw max %.1f | median total %.0f' % (X.max(), np.median(X.sum(1))), flush=True)
tot = X.sum(axis=1, keepdims=True)
X = X / np.maximum(tot, 1) * 1e4
X = np.log1p(X)
print('after log-normalisation: max %.2f' % X.max(), flush=True)

adata = ad.AnnData(X=X.astype(np.float32),
                   obs=meta.copy(), var=pd.DataFrame(index=cnt.columns))
adata.var_names = cnt.columns.astype(str)
scdrs.preprocess(adata, n_mean_bin=20, n_var_bin=20, copy=False)
print('preprocessed', flush=True)

sc = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000)
sc.to_csv(os.path.join(OUT, 'scdrs_GSE104782_cell_scores.csv'))
print('scored:', sc.shape, flush=True)

adata.obs = adata.obs.join(sc)
col = 'norm_score'

g = adata.obs.groupby('Cluster')[col].agg(['count', 'mean', 'median']).sort_values('mean',
                                                                                   ascending=False)
print('\n=== GSE104782: mean normalised scDRS by state ===')
print(g.round(4).to_string(), flush=True)

top = g.index[0]
print('\nhighest state:', top, flush=True)
obs = adata.obs
obs['is_top'] = (obs.Cluster == top).astype(int)

per = obs.groupby('donor').apply(
    lambda x: pd.Series({'top': x.loc[x.is_top == 1, col].mean(),
                         'other': x.loc[x.is_top == 0, col].mean(),
                         'n_top': int(x.is_top.sum())}))
per['delta'] = per.top - per.other
per = per.dropna()
print('\n=== donor-level: %s vs other cells ===' % top)
print(per.round(3).to_string(), flush=True)

if len(per) >= 3:
    w = stats.wilcoxon(per.top, per.other)
    t = stats.ttest_rel(per.top, per.other)
    print('\ndonors positive: %d/%d' % (int((per.delta > 0).sum()), len(per)))
    print('mean delta %.3f (SD %.3f) | Wilcoxon P = %.4g | t P = %.4g' % (
        per.delta.mean(), per.delta.std(ddof=1), w.pvalue, t.pvalue), flush=True)
    ci = stats.t.interval(0.95, len(per) - 1, loc=per.delta.mean(), scale=stats.sem(per.delta))
    print('95%% CI %.3f to %.3f' % ci, flush=True)

per.to_csv(os.path.join(OUT, 'scDRS_GSE104782_top_state_per_donor.csv'))
g.to_csv(os.path.join(OUT, 'scDRS_GSE104782_by_state.csv'))
print('\nwritten tables', flush=True)
