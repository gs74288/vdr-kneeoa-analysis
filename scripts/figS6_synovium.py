# -*- coding: utf-8 -*-
"""Supplementary figure S6: paired synovium-cartilage VDR detection by donor."""
import sys, io, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7,
                     'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6})
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"
OUT = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\Figures"

syn = pd.read_csv(os.path.join(T, 'GSE152805_synovium_clusters.csv'))
car = pd.read_csv(os.path.join(T, 'GSE152805_VDR_per_donor_state.csv'))

mac = syn[syn.cell_class == 'macrophage'].set_index('donor')['pct']
fib = syn[syn.cell_class == 'fibroblast'].set_index('donor')['pct']
donors = [d for d in ['113', '116', '118'] if d in mac.index and d in fib.index]

car_d = car.groupby('donor').apply(lambda g: 100 * g.detected.sum() / g.n.sum())

fig, (axA, axB) = plt.subplots(1, 2, figsize=(6.7, 2.9))

# A: paired macrophage vs fibroblast
x = np.arange(len(donors))
for i, d in enumerate(donors):
    axA.plot([0, 1], [fib[d], mac[d]], color='#B0B0B0', lw=0.9, zorder=1)
    axA.scatter([0], [fib[d]], s=16, color='#7F8C8D', zorder=3, linewidths=0)
    axA.scatter([1], [mac[d]], s=16, color='#B03A2E', zorder=3, linewidths=0)
    axA.text(1.06, mac[d], d, fontsize=6, va='center', color='#B03A2E')
axA.set_xticks([0, 1])
axA.set_xticklabels(['Fibroblast', 'Macrophage'])
axA.set_ylabel('VDR detected (% of cells)')
axA.set_xlim(-0.35, 1.35)
axA.set_title('Synovium, same donors', pad=4)
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)

# B: cartilage overall vs synovial macrophage
w = 0.36
for i, d in enumerate(donors):
    axB.bar(i - w / 2, car_d.get(d, np.nan), width=w, color='#2C3E50', alpha=0.75)
    axB.bar(i + w / 2, mac[d], width=w, color='#B03A2E', alpha=0.75)
axB.set_xticks(range(len(donors)))
axB.set_xticklabels(donors)
axB.set_xlabel('Donor')
axB.set_ylabel('VDR detected (% of cells)')
axB.set_title('Cartilage vs synovial macrophages', pad=4)
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)
from matplotlib.patches import Patch
axB.legend(handles=[Patch(color='#2C3E50', alpha=0.75, label='Cartilage (all states)'),
                    Patch(color='#B03A2E', alpha=0.75, label='Synovial macrophages')],
           frameon=False, fontsize=6, loc='upper left')

for ax, letter in ((axA, 'A'), (axB, 'B')):
    ax.text(-0.16, 1.06, letter, transform=ax.transAxes, fontsize=10, fontweight='bold')

fig.tight_layout()
png = os.path.join(OUT, 'S6_Fig.png')
tif = os.path.join(OUT, 'S6_Fig.tif')
fig.savefig(png, dpi=300, facecolor='white')
fig.savefig(tif, dpi=300)
plt.close(fig)

from PIL import Image
im = Image.open(tif).convert('RGB')
w0, h0 = im.size
if w0 != 2010:
    im = im.resize((2010, int(h0 * 2010 / w0)), Image.LANCZOS)
im.save(tif, format='TIFF', dpi=(300, 300), compression='tiff_lzw')
print('written S6_Fig.tif', im.size, '%.2f MB' % (os.path.getsize(tif) / 1e6))
print('\ncartilage overall per donor:', {k: round(v, 1) for k, v in car_d.items()})
print('synovial macrophage:', {k: round(v, 1) for k, v in mac.items()})
