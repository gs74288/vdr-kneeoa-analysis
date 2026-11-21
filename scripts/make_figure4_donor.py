# -*- coding: utf-8 -*-
"""Figure 4 as three panels A-C (panel D duplicated panel C)."""
import sys, io, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7,
                     'axes.linewidth': 0.7, 'axes.labelsize': 7.2})

TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"
OUT = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\Figures"

ab = pd.read_csv(os.path.join(TAB, 'GSE324993_VDR_per_donor.csv'))
eff = pd.read_csv(os.path.join(TAB, 'fig4_gse152805_within_state.csv'))

fig = plt.figure(figsize=(6.70, 2.90))
axA = fig.add_axes([0.105, 0.200, 0.215, 0.690])
axB = fig.add_axes([0.375, 0.200, 0.195, 0.690])
axC = fig.add_axes([0.665, 0.200, 0.300, 0.690])

C1, C2 = '#4E79A7', '#F28E2B'

# A
for i, r in ab.iterrows():
    axA.plot([0, 1], [r.other_frac, r.reg_frac], color='0.78', lw=0.7, zorder=1)
axA.scatter(np.zeros(len(ab)), ab.other_frac, s=22, c=C1, zorder=3)
axA.scatter(np.ones(len(ab)), ab.reg_frac, s=22, c=C2, zorder=3)
axA.set_xticks([0, 1]); axA.set_xticklabels(['other', 'RegC'], fontsize=6.4)
axA.set_xlim(-0.35, 1.35); axA.set_ylim(0.08, 0.28)
axA.set_ylabel('VDR-detected cells (%)', fontsize=6.8)
axA.set_yticks([0.1, 0.15, 0.2, 0.25])
axA.set_yticklabels(['10', '15', '20', '25'], fontsize=6.2)
axA.tick_params(labelsize=6.2)
for sp in ('top', 'right'):
    axA.spines[sp].set_visible(False)

# B
yB = np.arange(len(ab))
cols = [C2 if v > 0 else C1 for v in ab.delta]
axB.axvline(0, color='0.45', lw=0.7, ls=(0, (4, 3)))
axB.hlines(yB, 0, ab.delta.values * 100, color=cols, lw=1.0)
axB.scatter(ab.delta.values * 100, yB, s=20, c=cols, zorder=3)
axB.set_yticks(yB)
axB.set_yticklabels(['%s (%s)' % (s, g) for s, g in zip(ab['sample'], ab['group'])],
                    fontsize=5.8)
axB.invert_yaxis()
axB.set_xlabel('difference (pp)', fontsize=6.6)
axB.tick_params(labelsize=5.8)
for sp in ('top', 'right'):
    axB.spines[sp].set_visible(False)

# C: signature shifts (three signatures, six states)
pv = eff.pivot(index='state', columns='signature', values='diff')
order = list(pv.index)
xx = np.arange(len(order))
wd = 0.26
handles = []
for i, (k, c) in enumerate([('regulatory', '#B2182B'), ('catabolic', '#F28E2B'),
                            ('NLRP3', '#4E79A7')]):
    b = axC.bar(xx + (i - 1) * wd, pv[k].reindex(order).values, wd, color=c, label=k, alpha=0.9)
    handles.append(b)
axC.axhline(0, color='0.4', lw=0.7)
axC.set_xticks(xx); axC.set_xticklabels(order, fontsize=5.8, rotation=35, ha='right')
axC.set_ylabel('VDR+ minus VDR- (z units)', fontsize=6.8)
axC.tick_params(labelsize=6.0)
axC.legend(fontsize=5.4, frameon=False, loc='upper left', handlelength=1.0,
           handletextpad=0.3, borderpad=0.1, ncol=1)
for sp in ('top', 'right'):
    axC.spines[sp].set_visible(False)

for letter, ax, dx in (('A', axA, -0.22), ('B', axB, -0.40), ('C', axC, -0.20)):
    ax.text(dx, 1.04, letter, transform=ax.transAxes, fontsize=11, fontweight='bold',
            ha='left', va='bottom')

png = os.path.join(OUT, 'Figure_4_3p.png')
tif = os.path.join(OUT, 'Figure_4_3p.tif')
fig.savefig(png, dpi=300)
fig.savefig(tif, dpi=300)
plt.close(fig)

from PIL import Image
im = Image.open(tif)
if im.mode != 'RGB':
    rgb = im.convert('RGB'); im.close(); rgb.save(tif, dpi=(300, 300), compression='tiff_lzw'); rgb.close()
print('written:', tif, Image.open(tif).size)
