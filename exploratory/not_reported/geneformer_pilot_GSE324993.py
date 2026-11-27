# -*- coding: utf-8 -*-
"""Pilot: build a 1,000-cell stratified subset of GSE324993 and tokenize with Geneformer V2."""
import sys, io, os, re, inspect
import numpy as np
import pandas as pd
import anndata as ad

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
OUT = r"C:\Users\liuha\OpenClawVDR\data\gf_pilot"
os.makedirs(OUT, exist_ok=True)

a = ad.read_h5ad(os.path.join(D, 'GSE324993.h5ad'))
states = np.array([re.sub(r'\(.*\)', '', str(s)).strip() for s in a.obs['cell_type'].values])
rng = np.random.RandomState(0)

idx = []
for st in sorted(set(states)):
    pool = np.where(states == st)[0]
    k = min(len(pool), 125)
    idx.extend(rng.choice(pool, k, replace=False))
idx = np.array(sorted(idx))
sub = a[idx].copy()
sub.obs['state'] = states[idx]
sub.obs['donor'] = [s.split('_')[0] for s in sub.obs_names]
print('subset:', sub.shape, '| states:', dict(pd.Series(states[idx]).value_counts()))

# Geneformer expects raw counts; keep X as is but report range
X = sub.X
print('X dtype', type(X), '| min %.2f max %.2f' % (X.min(), X.max()))
sub.write(os.path.join(OUT, 'pilot.h5ad'))
print('written pilot.h5ad')

import geneformer
print('\nTokenizer signature:', inspect.signature(geneformer.Tokenizer.__init__))
