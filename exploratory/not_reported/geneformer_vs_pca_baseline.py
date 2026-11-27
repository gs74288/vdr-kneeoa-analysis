# -*- coding: utf-8 -*-
"""Baseline comparison: does Geneformer beat plain PCA on log-normalized counts?"""
import sys, io
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

for tag, h5, emb in [
    ('GSE104782', r"C:\Users\liuha\OpenClawVDR\data\gf_104782\raw\pilot.h5ad",
     r"C:\Users\liuha\OpenClawVDR\data\gf_104782\embs\pilot.csv"),
]:
    a = ad.read_h5ad(h5)
    X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
    # library-size normalize to 10k + log1p (standard scRNA preprocessing)
    lib = X.sum(1, keepdims=True)
    lib[lib == 0] = 1
    Xn = np.log1p(X / lib * 1e4)
    # top-variable genes
    var = Xn.var(0)
    top = np.argsort(var)[-2000:]
    Xh = Xn[:, top]
    Z = PCA(n_components=50, random_state=0).fit_transform(Xh)
    y = a.obs['state'].astype(str).values
    acc = cross_val_score(KNeighborsClassifier(15), StandardScaler().fit_transform(Z), y, cv=5).mean()
    print('%s  PCA50 baseline: kNN %.3f | silhouette %.3f | chance %.3f'
          % (tag, acc, silhouette_score(Z, y), 1 / len(set(y))))

    d = pd.read_csv(emb)
    ec = [c for c in d.columns if c not in ('Unnamed: 0', 'state', 'donor')]
    E = StandardScaler().fit_transform(d[ec].values.astype(float))
    y2 = d['state'].astype(str).values
    print('%s  Geneformer  : kNN %.3f | silhouette %.3f'
          % (tag, cross_val_score(KNeighborsClassifier(15), E, y2, cv=5).mean(), silhouette_score(E, y2)))
