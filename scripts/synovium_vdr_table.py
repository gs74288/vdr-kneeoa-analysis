# -*- coding: utf-8 -*-
"""Tabulate synovium VDR detection with per-sample Fisher tests."""
import sys, io, os
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

rows = [
    ('Infla1', 'GSM8655601', 'macrophage-like', 4645, 1946),
    ('Infla1', 'GSM8655601', 'fibroblast-like', 2076, 370),
    ('Infla1', 'GSM8655601', 'other', 291, 14),
    ('Infla2', 'GSM8655602', 'macrophage-like', 5404, 1136),
    ('Infla2', 'GSM8655602', 'fibroblast-like', 4375, 325),
    ('Infla2', 'GSM8655602', 'other', 912, 19),
]
df = pd.DataFrame(rows, columns=['sample', 'gsm', 'class', 'n', 'detected'])
df['pct'] = (100 * df.detected / df.n).round(1)

fisher = []
for s in df['sample'].unique():
    d = df[df['sample'] == s].set_index('class')
    tab = [[d.loc['macrophage-like', 'detected'], d.loc['macrophage-like', 'n'] - d.loc['macrophage-like', 'detected']],
           [d.loc['fibroblast-like', 'detected'], d.loc['fibroblast-like', 'n'] - d.loc['fibroblast-like', 'detected']]]
    orr, p = stats.fisher_exact(tab)
    ratio = d.loc['macrophage-like', 'pct'] / d.loc['fibroblast-like', 'pct']
    fisher.append((s, round(ratio, 2), float('%.3g' % orr), float('%.3g' % p)))
    print('%s: ratio %.2f | OR %.2f | Fisher P=%.3g' % (s, ratio, orr, p))

print()
print(df.to_string(index=False))
df.to_csv(os.path.join(T, 'synovium_vdr_detection.csv'), index=False)
pd.DataFrame(fisher, columns=['sample', 'ratio', 'OR', 'fisher_p']).to_csv(
    os.path.join(T, 'synovium_vdr_tests.csv'), index=False)
print('\nwritten synovium tables')
