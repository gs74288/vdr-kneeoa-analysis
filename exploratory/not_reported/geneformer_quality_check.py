# -*- coding: utf-8 -*-
"""Quality check: do Geneformer embeddings separate chondrocyte states?"""
import sys, io, os
import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
P = r"C:\Users\liuha\OpenClawVDR\data\gf_pilot\embs\pilot.csv"
d = pd.read_csv(P)
print('shape:', d.shape)
meta = [c for c in d.columns if not c.replace('.', '').isdigit()]
print('non-numeric cols:', meta[:5])

emb_cols = [c for c in d.columns if c not in meta]
X = d[emb_cols].values.astype(float)
print('embeddings:', X.shape)

y = d['state'].astype(str).values
Xs = StandardScaler().fit_transform(X)
knn = KNeighborsClassifier(n_neighbors=15)
acc = cross_val_score(knn, Xs, y, cv=5).mean()
print('5-fold kNN accuracy for state:', round(acc, 3))

sil = silhouette_score(Xs, y)
print('silhouette by state:', round(sil, 3))

print('\nper-state centroid distance check (mean pairwise):')
cats = sorted(set(y))
cent = {c: X[y == c].mean(0) for c in cats}
import itertools
rows = []
for a, b in itertools.combinations(cats, 2):
    rows.append((a, b, float(np.linalg.norm(cent[a] - cent[b]))))
rows.sort(key=lambda r: r[2])
print('  closest pairs:', [(a, b, round(v, 1)) for a, b, v in rows[:5]])
print('  farthest pairs:', [(a, b, round(v, 1)) for a, b, v in rows[-3:]])
