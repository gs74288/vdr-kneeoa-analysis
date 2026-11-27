# -*- coding: utf-8 -*-
"""MOA family enrichment with a tie-free ranking (significance-based)."""
import sys, io, os
import pandas as pd
import numpy as np
import gseapy as gp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

res = pd.read_csv(os.path.join(T, 'lincs_RepC_signature.csv')).dropna(subset=['Adjusted P-value'])
res['Term'] = res['Term'].astype(str)
res['score'] = -np.log10(res['Adjusted P-value'].clip(lower=1e-300)) + \
               res['Combined Score'] / res['Combined Score'].max()
res = res.drop_duplicates('Term')
print('terms:', len(res), '| unique scores:', res.score.nunique())

FAMILIES = {
    'VDR_agonist': r'calcitriol|ercalcitriol|alfacalcidol|paricalcitol|doxercalciferol|vitamin d',
    'ROCK_inhibitor': r'y-?27632|fasudil|gsk429286|rockout|sr3677|thiazovivin',
    'SRC_inhibitor': r'dasatinib|saracatinib|bosutinib|pp2|pp1|kx2-391|azd0530',
    'YAP_TEAD_inhibitor': r'verteporfin|super-tdu|\bca3\b|yap',
    'TGFb_inhibitor': r'sb-?525334|sb-?431542|ly2157299|galunisertib',
    'MEK_inhibitor': r'trametinib|pd-?0325901|u0126|selumetinib',
}
sets = {}
for fam, pat in FAMILIES.items():
    for direction in ('Up', 'Down'):
        names = res[res.Term.str.contains(pat, case=False, regex=True) &
                    res.Term.str.endswith(direction)]['Term'].tolist()
        if names:
            sets['%s_%s' % (fam, direction)] = names

rnk = res.set_index('Term')['score'].sort_values(ascending=False)
print('families:', {k: len(v) for k, v in sets.items()})

pre = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=2000, seed=0,
                 outdir=None, min_size=1, max_size=500, threads=4)
out = pre.res2d[['Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val', 'Lead_genes']]
print('\n=== MOA family enrichment (significance-ranked) ===')
print(out.sort_values('NES', ascending=False).round(3).to_string(index=False))
out.to_csv(os.path.join(T, 'lincs_MOA_families_ranked.csv'), index=False)
print('\nwritten lincs_MOA_families_ranked.csv')
