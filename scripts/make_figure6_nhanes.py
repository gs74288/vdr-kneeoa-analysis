# -*- coding: utf-8 -*-
"""Figure 6 rebuilt from the NHANES design file (no in-panel titles, journal width)."""
import sys, io, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7,
                     'axes.linewidth': 0.7, 'xtick.major.width': 0.7,
                     'ytick.major.width': 0.7, 'axes.labelsize': 7.5})

TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables\fig6_data"
OUT = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\Figures"

A = pd.read_csv(os.path.join(TAB, 'panelA_mean_by_grade.csv'))
B = pd.read_csv(os.path.join(TAB, 'panelB_progressive_OR.csv'))

fig = plt.figure(figsize=(6.70, 3.05))
axA = fig.add_axes([0.085, 0.175, 0.375, 0.700])
axB = fig.add_axes([0.615, 0.175, 0.360, 0.700])

# ---------- panel A ----------
C = '#2166AC'
axA.errorbar(A.grade, A['mean'], yerr=1.96 * A.se, fmt='o', ms=3.6,
             color=C, ecolor=C, elinewidth=0.8, capsize=2.0, capthick=0.8, zorder=3)
axA.plot(A.grade, A['mean'], color=C, lw=0.7, alpha=0.55, zorder=2)
axA.set_xticks(A.grade)
axA.set_xticklabels(['%d\n(n=%d)' % (g, n) for g, n in zip(A.grade, A.n)], fontsize=6.2)
axA.set_xlabel('Kellgren-Lawrence grade', fontsize=7.2)
axA.set_ylabel('Serum 25(OH)D (nmol/L)', fontsize=7.2)
axA.set_ylim(40, 68)
axA.tick_params(labelsize=6.4)
for sp in ('top', 'right'):
    axA.spines[sp].set_visible(False)

# ---------- panel B ----------
B = B.iloc[::-1].reset_index(drop=True)      # crude at the bottom of the axis -> top after invert
y = np.arange(len(B))
axB.axvline(1.0, color='0.45', ls=(0, (4, 3)), lw=0.7, zorder=1)
axB.errorbar(B.OR, y, xerr=[B.OR - B.lo, B.hi - B.OR], fmt='o', ms=4.0,
             color='#B2182B', ecolor='#B2182B', elinewidth=0.9, capsize=2.2,
             capthick=0.9, zorder=3)
axB.set_yticks(y)
axB.set_yticklabels(B.model, fontsize=6.4)
axB.invert_yaxis()
axB.set_xscale('log')
axB.set_xlim(0.62, 1.45)
axB.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4])
axB.set_xticklabels(['0.7', '0.8', '0.9', '1.0', '1.1', '1.2', '1.3', '1.4'], fontsize=6.2)
axB.set_xlabel('OR per 25 nmol/L (95% CI)', fontsize=7.2)
axB.tick_params(labelsize=6.4)
for sp in ('top', 'right'):
    axB.spines[sp].set_visible(False)

for letter, ax, dx in (('A', axA, -0.19), ('B', axB, -0.30)):
    ax.text(dx, 1.06, letter, transform=ax.transAxes, fontsize=11,
            fontweight='bold', ha='left', va='bottom')

os.makedirs(OUT, exist_ok=True)
png = os.path.join(OUT, 'Figure_6_new.png')
tif = os.path.join(OUT, 'Figure_6_new.tif')
fig.savefig(png, dpi=300)
fig.savefig(tif, dpi=300)
plt.close(fig)

from PIL import Image
im = Image.open(tif)
if im.mode != 'RGB':
    rgb = im.convert('RGB')
    im.close()
    rgb.save(tif, dpi=(300, 300), compression='tiff_lzw')
    rgb.close()
im2 = Image.open(tif)
print('written:', tif, im2.size, im2.mode, round(os.path.getsize(tif) / 1e6, 2), 'MB')
print('\npanel A range: %.1f - %.1f nmol/L' % (A['mean'].min(), A['mean'].max()))
print('panel B ORs  : %.3f -> %.3f' % (B.OR.iloc[-1], B.OR.iloc[0]))
