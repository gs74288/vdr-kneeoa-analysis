# -*- coding: utf-8 -*-
"""pySCENIC on GSE152805: infer regulons, then score VDR activity per cell.

GRN inference is run on a random subsample of cells; AUCell scoring is applied to all cells.
"""
import sys, io, os, time
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\tools\scenic"
OUT = os.path.join(D, 'results')
os.makedirs(OUT, exist_ok=True)
N_SUB = 8000
SEED = 1

t0 = time.time()
print('loading matrix...', flush=True)
X = scipy.io.mmread(os.path.join(D, 'counts.mtx')).tocsr()   # genes x cells
genes = open(os.path.join(D, 'genes.txt')).read().split('\n')
genes = [g for g in genes if g]
cells = open(os.path.join(D, 'cells.txt')).read().split('\n')
cells = [c for c in cells if c]
print('matrix: %s genes x %s cells (%.0fs)' % (X.shape[0], X.shape[1], time.time() - t0), flush=True)

meta = pd.read_csv(os.path.join(D, 'cellmeta.csv')).set_index('cell')

# filter genes
n_cells_expr = np.asarray((X > 0).sum(axis=1)).ravel()
keep = n_cells_expr >= 30
X = X[keep]
genes = [g for g, k in zip(genes, keep) if k]
print('genes after filter:', len(genes), flush=True)

# normalise (CPM-ish) and log
tot = np.asarray(X.sum(axis=0)).ravel()
tot[tot == 0] = 1
Xn = X.multiply(1e4 / tot).tocsc()
Xn.data = np.log1p(Xn.data)
Xc = Xn.T.tocsr()          # cells x genes
print('normalised', flush=True)

rng = np.random.RandomState(SEED)
idx = rng.choice(Xc.shape[0], size=min(N_SUB, Xc.shape[0]), replace=False)
ex = pd.DataFrame(Xc[idx].toarray(), columns=genes)

tfs = [l.strip() for l in open(os.path.join(D, 'allTFs_hg38.txt')) if l.strip()]
tfs = [t for t in tfs if t in ex.columns]
print('TFs present:', len(tfs), flush=True)

from arboreto.algo import grnboost2
from distributed import Client, LocalCluster
cluster = LocalCluster(n_workers=4, threads_per_worker=1, processes=False,
                       memory_limit=None, dashboard_address=None)
client = Client(cluster)
print('running GRNBoost2 on %d cells (threads)...' % ex.shape[0], flush=True)
adj = grnboost2(ex, tf_names=tfs, client_or_address=client, verbose=False, seed=SEED)
adj.to_csv(os.path.join(OUT, 'adjacencies.tsv'), sep='\t', index=False)
print('adjacencies:', adj.shape, flush=True)

from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase

fdb = RankingDatabase(os.path.join(D, 'hg38_10k.feather'), 'hg38_10k')
modules = list(modules_from_adjacencies(adj, ex))
print('modules:', len(modules), flush=True)
motif_ann = os.path.join(D, 'motifs.tbl')
df = prune2df([fdb], modules, motif_ann, num_workers=4)
df.to_csv(os.path.join(OUT, 'regulons.csv'))
print('regulons:', df.shape, flush=True)

from pyscenic.aucell import aucell
regulons = {}
for name in df.index:
    if 'VDR' in str(name):
        targets = df.loc[name]
        hit = [g for g in targets.index if isinstance(targets[g], str) and targets[g]]
        if hit:
            regulons[str(name)] = hit
print('VDR regulons found:', list(regulons.keys()), flush=True)

if regulons:
    auc = aucell(Xc[idx], regulons, num_workers=4)
    auc.to_csv(os.path.join(OUT, 'aucell_subset.csv'))

import json
json.dump({k: len(v) for k, v in regulons.items()}, open(os.path.join(OUT, 'vdr_regulon_sizes.json'), 'w'))
print('DONE_SCENIC %.0fs' % (time.time() - t0), flush=True)
