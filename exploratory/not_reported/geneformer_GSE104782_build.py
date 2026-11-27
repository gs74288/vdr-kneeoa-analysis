# -*- coding: utf-8 -*-
"""GSE104782 validation pilot: raw UMI counts -> 1,000-cell stratified subset."""
import sys, io, os, gzip
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data\GSE104782"
OUT = r"C:\Users\liuha\OpenClawVDR\data\gf_104782"
os.makedirs(OUT, exist_ok=True)

# --- metadata: cell -> cluster
meta = pd.read_excel(os.path.join(D, 'GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx'))
meta = meta.dropna(subset=['Cluster'])
cell2clu = dict(zip(meta['Cell'].astype(str), meta['Cluster'].astype(str)))
print('cells with cluster label:', len(cell2clu), '| clusters:', dict(meta['Cluster'].value_counts()))

# --- counts
path = os.path.join(D, 'GSE104782_allcells_UMI_count.txt.gz')
with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
    header = fh.readline().rstrip('\n').split('\t')
    genes = []
    rows = []
    for ln in fh:
        p = ln.rstrip('\n').split('\t')
        genes.append(p[0])
        rows.append(np.array(p[1:], dtype=np.int32))
cells = header[1:]
X = np.vstack(rows)
print('matrix: %d genes x %d cells' % X.shape)
print('X min %d max %d | integer %s' % (X.min(), X.max(), bool(np.all(X == np.round(X)))))

cell_idx = [i for i, c in enumerate(cells) if c in cell2clu]
print('cells matched to metadata:', len(cell_idx))
lab = np.array([cell2clu[cells[i]] for i in cell_idx])

rng = np.random.RandomState(0)
pick = []
for st in sorted(set(lab)):
    pool = np.where(lab == st)[0]
    k = min(len(pool), 140)
    pick.extend(rng.choice(pool, k, replace=False))
pick = np.array(sorted(pick))
print('subset:', len(pick), '| per-state:', dict(pd.Series(lab[pick]).value_counts()))

Xs = X[:, np.array(cell_idx)[pick]]
a = ad.AnnData(X=sp.csr_matrix(Xs.T.astype(np.float32)))
a.obs_names = [cells[cell_idx[i]] for i in pick]
a.var_names = genes
a.var_names_make_unique()
a.obs['state'] = lab[pick]
a.obs['donor'] = [c.split('_')[0] for c in a.obs_names]
print('adata:', a.shape, '| nnz %.1f%%' % (100 * a.X.nnz / (a.n_obs * a.n_vars)))
a.write(os.path.join(OUT, 'pilot.h5ad'))
print('written', os.path.join(OUT, 'pilot.h5ad'))
