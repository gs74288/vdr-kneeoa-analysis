# -*- coding: utf-8 -*-
"""Figure 7: MOA family enrichment and transcriptional decoupling in bulk cartilage."""
import sys, io, os, re
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7,
                     'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6})

T = r"C:\Users\liuha\OpenClawVDR\repo\tables"
D = r"C:\Users\liuha\OpenClawVDR\bulk"
OUT = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\Figures"

CORE = ['GPX3', 'IGFBP3', 'MGP', 'PLPP3', 'PTGIS', 'TXNRD1']


def load(fn):
    d = pd.read_csv(os.path.join(D, fn), sep='\t', index_col=0)
    keep = [c for c in d.columns if not re.search(r'^(average|max)', str(c).strip(), re.I)]
    return d.loc[:, keep].groupby(level=0).mean()


oa, no = load('GSE114007_OA_counts.txt.gz'), load('GSE114007_normal_counts.txt.gz')
X = pd.concat([oa, no], axis=1)
vdr = X.loc['VDR'].astype(float)

present = [g for g in CORE if g in X.index]
Z = X.loc[present].astype(float)
Z = Z.sub(Z.mean(1), axis=0).div(Z.std(1).replace(0, np.nan), axis=0)
score = Z.mean(0)

r_all = stats.spearmanr(vdr, score)
r_oa = stats.spearmanr(vdr[oa.columns], score[oa.columns])
print('core genes used:', present)
print('all: rho %.3f P %.4g | OA: rho %.3f P %.4g' % (r_all.statistic, r_all.pvalue,
                                                      r_oa.statistic, r_oa.pvalue))

moa = pd.read_csv(os.path.join(T, 'lincs_MOA_families_ranked.csv'))
sel = ['VDR_agonist_Up', 'SRC_inhibitor_Down', 'TGFb_inhibitor_Down', 'MEK_inhibitor_Up',
       'ROCK_inhibitor_Down', 'YAP_TEAD_inhibitor_Down']
m = moa[moa.Term.isin(sel)].copy().sort_values('NES')

fig = plt.figure(figsize=(6.7, 2.6))
axA = fig.add_axes([0.30, 0.22, 0.26, 0.68])
axB = fig.add_axes([0.68, 0.22, 0.13, 0.68])
axC = fig.add_axes([0.86, 0.22, 0.13, 0.68])

# A: NES bars
colors = ['#B03A2E' if f < 0.05 else '#B0B0B0' for f in m['FDR q-val']]
axA.barh(range(len(m)), m['NES'], color=colors, height=0.62)
axA.set_yticks(range(len(m)))
axA.set_yticklabels([t.replace('_', ' ') for t in m['Term']], fontsize=6)
axA.axvline(0, color='k', lw=0.5)
axA.set_xlabel('Normalised enrichment score')
for i, (nes, fdr) in enumerate(zip(m['NES'], m['FDR q-val'])):
    axA.text(nes + 0.03, i, 'FDR %.3f' % fdr, va='center', fontsize=5.6)
axA.set_xlim(0, max(m['NES']) * 1.45)
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)
axA.set_title('MOA family enrichment', pad=4)

# B: scatter VDR vs core score (OA)
x = vdr[oa.columns].values
y = score[oa.columns].values
axB.scatter(x, y, s=12, color='#2C3E50', linewidths=0, alpha=0.9)
b, a = np.polyfit(x, y, 1)
xs = np.linspace(x.min(), x.max(), 30)
axB.plot(xs, a + b * xs, color='#B03A2E', lw=0.9)
axB.set_xlabel('VDR (log2)', fontsize=6.5)
axB.set_ylabel('Calcitriol-response score', fontsize=6.5)
axB.text(0.04, 0.06, 'rho = %.2f\nP = %.3f' % (r_oa.statistic, r_oa.pvalue),
         transform=axB.transAxes, fontsize=6)
axB.tick_params(labelsize=5.8)
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)

# C: per-gene rho with VDR
rhos = []
for g in present:
    r = stats.spearmanr(vdr[oa.columns], X.loc[g, oa.columns].astype(float))
    rhos.append((g, r.statistic))
rhos.sort(key=lambda t: t[1])
names = [t[0] for t in rhos]
vals = [t[1] for t in rhos]
axC.barh(range(len(vals)), vals, color='#7F8C8D', height=0.62)
axC.set_yticks(range(len(vals)))
axC.set_yticklabels(names, fontsize=6)
axC.axvline(0, color='k', lw=0.5)
axC.set_xlabel('Spearman rho with VDR', fontsize=6.5)
axC.tick_params(labelsize=5.8)
axC.spines['top'].set_visible(False)
axC.spines['right'].set_visible(False)

for ax, letter, dx in ((axA, 'A', -0.42), (axB, 'B', -0.30), (axC, 'C', -0.30)):
    ax.text(dx, 1.05, letter, transform=ax.transAxes, fontsize=10, fontweight='bold')

png = os.path.join(OUT, 'Figure_7.png')
tif = os.path.join(OUT, 'Figure_7.tif')
fig.savefig(png, dpi=300, facecolor='white')
fig.savefig(tif, dpi=300)
plt.close(fig)

from PIL import Image
im = Image.open(tif).convert('RGB')
w, h = im.size
if w != 2010:
    im = im.resize((2010, int(h * 2010 / w)), Image.LANCZOS)
im.save(tif, format='TIFF', dpi=(300, 300), compression='tiff_lzw')
print('written Figure_7.tif', im.size, '%.2f MB' % (os.path.getsize(tif) / 1e6))
print('per-gene rho:', [(n, round(v, 2)) for n, v in rhos])
