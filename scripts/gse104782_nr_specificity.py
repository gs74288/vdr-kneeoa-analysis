# -*- coding: utf-8 -*-
"""GSE104782: are other low-abundance nuclear receptors also enriched in RegC?
   And does the RegC excess survive depth stratification?"""
import sys, io, os
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

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
states = lab.set_index('Cell').loc[cells, 'Cluster'].values
ngene = (X > 0).sum(1).values

NRS = ['VDR', 'AR', 'PGR', 'ESR1', 'ESR2', 'RORA', 'RORB', 'RORC', 'NR4A1', 'NR4A2',
       'NR1H3', 'NR1H2', 'NR2F1', 'NR2F2', 'PPARG', 'PPARD', 'THRA', 'THRB', 'RARA', 'RARB',
       'NR2C1', 'NR2C2', 'NR5A1', 'NR5A2', 'NR0B1', 'NR0B2', 'NR1D1', 'NR1D2', 'NR1I2', 'NR1I3',
       'ESRRA', 'ESRRB', 'ESRRG', 'HNF4A', 'HNF4G', 'NR6A1', 'NR3C1', 'NR3C2', 'RXRA', 'RXRB',
       'RXRG', 'NR2E1', 'NR2E3', 'NR2F6', 'NR4A3', 'NR1H4']
present = [g for g in NRS if g in X.columns]
rows = []
reg = states == 'RegC'
for g in present:
    d = (X[g].values > 0)
    p_reg = d[reg].mean()
    p_oth = d[~reg].mean()
    if d.sum() >= 5:
        tab = [[d[reg].sum(), reg.sum() - d[reg].sum()],
               [d[~reg].sum(), (~reg).sum() - d[~reg].sum()]]
        orr, p = stats.fisher_exact(tab)
        rows.append(dict(gene=g, detected=int(d.sum()), pct_reg=100 * p_reg, pct_other=100 * p_oth,
                         ratio=(p_reg / p_oth) if p_oth > 0 else np.nan, OR=orr, p=p))
res = pd.DataFrame(rows).sort_values('p')
print('=== nuclear receptors detected in GSE104782 (>=5 cells) ===')
print(res.round(3).to_string(index=False))
print('\nreceptors enriched in RegC at P<0.05:', int((res.p < 0.05).sum()), 'of', len(res))
print('VDR rank by P:', int((res.p < res.loc[res.gene == "VDR", "p"].iloc[0]).sum()) + 1)

# depth stratification
print('\n=== RegC excess by detected-gene quartile ===')
q = pd.qcut(ngene, 4, labels=False)
out = []
for k in range(4):
    m = q == k
    a = (X['VDR'].values > 0)[m & reg]
    b = (X['VDR'].values > 0)[m & ~reg]
    if len(a) >= 5 and len(b) >= 5:
        out.append((k, int(m.sum()), 100 * a.mean(), 100 * b.mean(), 100 * (a.mean() - b.mean())))
print(pd.DataFrame(out, columns=['quartile', 'cells', 'RegC %', 'other %', 'difference pp']).round(1).to_string(index=False))
