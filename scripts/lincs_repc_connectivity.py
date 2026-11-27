# -*- coding: utf-8 -*-
"""LINCS connectivity for a RepC signature: calcitriol vs ROCK/YAP/STAT3 inhibitors."""
import sys, io, os
import numpy as np
import pandas as pd
import gseapy as gp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

pb = pd.read_csv(os.path.join(T, 'pseudobulk_gse324993_donor_state.csv'), index_col=0)
states = pd.Series([c.split('|')[1] for c in pb.columns], index=pb.columns)
rep = pb.loc[:, states == 'RepC'].mean(1)
oth = pb.loc[:, states != 'RepC'].mean(1)
frac = (pb > 0.1).mean(1)
d = (rep - oth)
tab = pd.DataFrame({'diff': d, 'frac': frac})
tab = tab[tab.frac > 0.5].sort_values('diff', ascending=False)
up = tab.head(150).index.tolist()
print('RepC up genes:', len(up))
print('  ', ', '.join(up[:20]))

res = gp.enrichr(gene_list=up, gene_sets='LINCS_L1000_Chem_Pert_Consensus_Sigs',
                 outdir=None, no_plot=True).results
print('\nLINCS hits:', len(res))
res['Term_l'] = res['Term'].str.lower()

WANT = ['calcitriol', 'alfacalcidol', 'vitamin d', 'y-27632', 'y27632', 'fasudil',
        'verteporfin', 'stattic', 'rock', 'yap', 'fak', 'stat3']
print('\n=== compounds of interest ===')
hit = res[res.Term_l.apply(lambda t: any(w in t for w in WANT))]
if len(hit):
    print(hit.sort_values('Adjusted P-value')[['Term', 'Adjusted P-value', 'Combined Score']]
          .head(20).to_string(index=False))
else:
    print('  none among the top-ranked terms')

print('\n=== top 15 overall ===')
print(res.sort_values('Adjusted P-value')[['Term', 'Adjusted P-value', 'Combined Score']]
      .head(15).to_string(index=False))
res.to_csv(os.path.join(T, 'lincs_RepC_signature.csv'), index=False)
