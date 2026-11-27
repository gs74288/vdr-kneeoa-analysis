# -*- coding: utf-8 -*-
"""Mechano modules in GSE104782: state pattern and correlation with scDRS."""
import sys, io, os
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

MODULES = {
    'YAP_TAZ_targets': ['CCN2', 'CTGF', 'CCN1', 'CYR61', 'ANKRD1', 'AMOTL2', 'TEAD1', 'TEAD4'],
    'mechanochannels': ['PIEZO1', 'PIEZO2', 'TRPV4', 'TRPC1', 'TRPC6'],
    'integrin_FAK_ROCK': ['ITGA5', 'ITGB1', 'PTK2', 'ROCK1', 'ROCK2', 'ACTA2', 'TAGLN', 'MYL9'],
    'ECM_stiffness': ['LOX', 'LOXL2', 'FN1', 'COL1A1', 'COL3A1', 'TIMP1', 'SPP1'],
}

cnt = pd.read_csv(os.path.join(D, 'GSE104782', 'GSE104782_allcells_UMI_count.txt.gz'),
                  sep='\t', index_col=0)
if cnt.shape[0] > cnt.shape[1]:
    cnt = cnt.T
cnt = cnt.groupby(level=0).sum()
lab = pd.read_excel(os.path.join(D, 'GSE104782',
                                 'GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx'))
lab = lab[['Cell', 'Cluster']].dropna()
cells = [c for c in lab.Cell if c in cnt.index]
X = cnt.loc[cells]
X = np.log1p(X.div(X.sum(1), axis=0) * 1e4)
states = lab.set_index('Cell').loc[cells, 'Cluster'].values
donors = np.array([c.split('_')[0] for c in cells])

sc = pd.read_csv(os.path.join(T, 'scdrs_GSE104782_cell_scores.csv')).rename(columns={'Cell': 'cell'})
sc = sc.set_index('cell').loc[cells]

rows = {}
for k, gs in MODULES.items():
    present = [g for g in gs if g in X.columns]
    Z = X[present]
    Z = Z.sub(Z.mean(1), axis=0).div(Z.std(1).replace(0, np.nan), axis=0)
    rows['sc_' + k] = Z.mean(1).values
df = pd.DataFrame(rows, index=cells)
df['state'] = states
df['donor'] = donors
df['norm_score'] = sc['norm_score'].values
df['mechano_mean'] = df[[c for c in df.columns if c.startswith('sc_')]].mean(1)

mech = ['sc_' + k for k in MODULES]
print('=== mean module by state (z across states) ===')
g = df.groupby('state')[mech + ['mechano_mean', 'norm_score']].mean()
print(((g - g.mean()) / g.std()).round(2).to_string())

print('\n=== donor-state correlation with scDRS ===')
ds = df.groupby(['donor', 'state'])[mech + ['mechano_mean', 'norm_score']].mean().dropna()
for m in mech + ['mechano_mean']:
    r = stats.spearmanr(ds[m], ds.norm_score)
    print('  %-24s rho=%+.3f P=%.4g (n=%d)' % (m, r.statistic, r.pvalue, len(ds)))

df.to_csv(os.path.join(T, 'GSE104782_mechano_scores.csv'))
print('\nwritten GSE104782_mechano_scores.csv')
