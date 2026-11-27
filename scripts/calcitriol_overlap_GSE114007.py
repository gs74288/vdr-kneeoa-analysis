# -*- coding: utf-8 -*-
"""Step 1: Calcitriol-overlap genes -> test against VDR in GSE114007; check VDR targeting."""
import sys, io, os, re
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"
D = r"C:\Users\liuha\OpenClawVDR\bulk"

res = pd.read_csv(os.path.join(T, 'lincs_RepC_signature.csv'))
res['Term'] = res['Term'].astype(str)

def genes_of(term):
    r = res[res.Term == term]
    if not len(r):
        return set()
    g = str(r.iloc[0]['Genes']) if 'Genes' in r.columns else ''
    return set(x.strip() for x in g.split(';') if x.strip())

vdr_up = genes_of('Calcitriol Up')
erc = genes_of('Ercalcitriol Up')
alf = genes_of('Alfacalcidol Up')
print('Calcitriol Up overlap genes:', len(vdr_up))
print('Ercalcitriol Up:', len(erc), '| Alfacalcidol Up:', len(alf))
core = vdr_up & erc
print('intersection Calcitriol & Ercalcitriol:', len(core))
print('  ', ', '.join(sorted(core))[:300])

# CollecTRI: are these direct VDR targets?
net = pd.read_csv(r"C:\Users\liuha\OpenClawVDR\data\collectri_human.csv")
tgt = set(net.loc[net.source == 'VDR', 'target'])
bg = set(net.target)
def frac_hit(s, ref):
    s = s & ref
    return len(s)
print('\nVDR-target enrichment among Calcitriol-Up genes:')
for label, s in (('calcitriol', vdr_up), ('core', core)):
    k, n = frac_hit(s, tgt), len(s)
    if n:
        p = stats.hypergeom.sf(k - 1, 20000, len(tgt), n)
        print('  %-12s %d/%d are VDR targets | hypergeom P=%.3g' % (label, k, n, p))

# GSE114007 test
def load(fn):
    d = pd.read_csv(os.path.join(D, fn), sep='\t', index_col=0)
    keep = [c for c in d.columns if not re.search(r'^(average|max)', str(c).strip(), re.I)]
    return d.loc[:, keep].groupby(level=0).mean()

oa, no = load('GSE114007_OA_counts.txt.gz'), load('GSE114007_normal_counts.txt.gz')
X = pd.concat([oa, no], axis=1)
vdr = X.loc['VDR'].astype(float)

for label, gs in (('calcitriol-core', sorted(core)),):
    present = [g for g in gs if g in X.index]
    if len(present) < 3:
        print('\n%s: too few genes in bulk (%d)' % (label, len(present)))
        continue
    Z = X.loc[present].astype(float)
    Z = Z.sub(Z.mean(1), axis=0).div(Z.std(1).replace(0, np.nan), axis=0)
    score = Z.mean(0)
    r_all = stats.spearmanr(vdr, score)
    r_oa = stats.spearmanr(vdr[oa.columns], score[oa.columns])
    print('\n%s: %d genes' % (label, len(present)))
    print('  vs VDR, all  : rho=%+.3f P=%.4g' % (r_all.statistic, r_all.pvalue))
    print('  vs VDR, OA   : rho=%+.3f P=%.4g' % (r_oa.statistic, r_oa.pvalue))
