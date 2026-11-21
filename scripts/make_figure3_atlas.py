# -*- coding: utf-8 -*-
"""Figure 3 rebuilt from GSE104782 (UMAP recomputed, published cluster labels)."""
import sys, io, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7,
                     'axes.linewidth': 0.7, 'axes.labelsize': 7.2})

TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"
OUT = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\Figures"

d = pd.read_csv(os.path.join(TAB, 'fig3_umap.csv'))
f = pd.read_csv(os.path.join(TAB, 'fig3_vdr_by_state.csv'), index_col=0)

STATES = ['HomC', 'RegC', 'preHTC', 'HTC', 'ProC', 'FC', 'EC']
COL = {'HomC': '#4E79A7', 'RegC': '#F28E2B', 'preHTC': '#59A14F', 'HTC': '#B07AA1',
       'ProC': '#76B7B2', 'FC': '#EDC948', 'EC': '#E15759'}

fig = plt.figure(figsize=(6.70, 5.95))
axA1 = fig.add_axes([0.055, 0.560, 0.415, 0.390])
axA2 = fig.add_axes([0.545, 0.560, 0.415, 0.390])
axB1 = fig.add_axes([0.055, 0.085, 0.415, 0.390])
axB2 = fig.add_axes([0.585, 0.085, 0.375, 0.390])

# ---- A1: states ----
for s in STATES:
    m = d.Cluster == s
    axA1.scatter(d.UMAP1[m], d.UMAP2[m], s=2.2, c=COL[s], alpha=0.85, linewidths=0,
                 label=s)
axA1.set_xlabel('UMAP 1', fontsize=6.8)
axA1.set_ylabel('UMAP 2', fontsize=6.8)
axA1.set_xticks([]); axA1.set_yticks([])
axA1.legend(markerscale=3.0, fontsize=5.2, frameon=False, loc='lower left',
            handletextpad=0.2, labelspacing=0.18, borderpad=0.15)

# ---- A2: donors ----
donors = sorted(d.donor.unique())
dcol = dict(zip(donors, ['#4E79A7', '#F28E2B', '#59A14F', '#B07AA1', '#76B7B2', '#EDC948']))
for dnr in donors:
    m = d.donor == dnr
    axA2.scatter(d.UMAP1[m], d.UMAP2[m], s=2.2, c=dcol.get(dnr, '0.5'), alpha=0.85,
                 linewidths=0, label=dnr)
axA2.set_xlabel('UMAP 1', fontsize=6.8)
axA2.set_ylabel('UMAP 2', fontsize=6.8)
axA2.set_xticks([]); axA2.set_yticks([])
axA2.legend(markerscale=3.0, fontsize=5.2, frameon=False, loc='lower left',
            handletextpad=0.2, labelspacing=0.18, borderpad=0.15)

# ---- B1: VDR detection overlay ----
neg = d.VDR_detected == 0
axB1.scatter(d.UMAP1[neg], d.UMAP2[neg], s=2.0, c='0.85', alpha=0.8, linewidths=0,
             label='VDR not detected')
pos = ~neg
axB1.scatter(d.UMAP1[pos], d.UMAP2[pos], s=6.0, c='#B2182B', alpha=0.95, linewidths=0,
             label='VDR detected (n=%d)' % pos.sum())
axB1.set_xlabel('UMAP 1', fontsize=6.8)
axB1.set_ylabel('UMAP 2', fontsize=6.8)
axB1.set_xticks([]); axB1.set_yticks([])
axB1.legend(markerscale=2.4, fontsize=5.8, frameon=False, loc='upper left',
            handletextpad=0.3)

# ---- B2: detection fraction by state ----
f = f.reindex(STATES)
y = np.arange(len(f))
axB2.barh(y, f.frac * 100, height=0.66, color=[COL[s] for s in f.index], alpha=0.9)
for i, (s, r) in enumerate(f.iterrows()):
    axB2.text(r.frac * 100 + 0.3, i, '%d/%d' % (r.det, r.n), va='center', fontsize=5.4)
axB2.set_yticks(y)
axB2.set_yticklabels(f.index, fontsize=6.2)
axB2.invert_yaxis()
axB2.set_xlabel('VDR-detected cells (%)', fontsize=6.8)
axB2.set_xlim(0, 19)
axB2.tick_params(labelsize=6.2)
for sp in ('top', 'right'):
    axB2.spines[sp].set_visible(False)

for letter, ax, dx, dy in (('A', axA1, -0.10, 1.03), ('B', axB1, -0.10, 1.03)):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=11, fontweight='bold',
            ha='left', va='bottom')

os.makedirs(OUT, exist_ok=True)
png = os.path.join(OUT, 'Figure_3_new.png')
tif = os.path.join(OUT, 'Figure_3_new.tif')
fig.savefig(png, dpi=300)
fig.savefig(tif, dpi=300)
plt.close(fig)

from PIL import Image
im = Image.open(tif)
if im.mode != 'RGB':
    rgb = im.convert('RGB'); im.close(); rgb.save(tif, dpi=(300, 300), compression='tiff_lzw'); rgb.close()
print('written:', tif, Image.open(tif).size)
