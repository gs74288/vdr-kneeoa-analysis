# -*- coding: utf-8 -*-
"""Figure 2 as three panels: scatter, forest, leave-one-out (funnel dropped: 3 variants)."""
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

d = pd.read_csv(os.path.join(TAB, 'cisMR_GCST90566800_instruments.csv'))
d = d.rename(columns={'snp': 'SNP'})
wald = d.bo.values / d.bx.values
wse = np.abs(d.seo.values / d.bx.values)
w = 1 / wse ** 2
b_ivw = (w * wald).sum() / w.sum()
se_fixed = np.sqrt(1 / w.sum())
q = (w * (wald - b_ivw) ** 2).sum()
dfq = len(d) - 1
phi = max(1.0, q / dfq)
se_rand = se_fixed * np.sqrt(phi)

fig = plt.figure(figsize=(6.70, 2.75))
axA = fig.add_axes([0.105, 0.215, 0.235, 0.660])
axB = fig.add_axes([0.395, 0.215, 0.265, 0.660])
axC = fig.add_axes([0.735, 0.215, 0.225, 0.660])
CR = '#B2182B'

# A: scatter
axA.axhline(0, color='0.78', lw=0.5, zorder=1)
axA.axvline(0, color='0.78', lw=0.5, zorder=1)
axA.errorbar(d.bx, d.bo, yerr=1.96 * d.seo, xerr=1.96 * d.sex, fmt='o', ms=4.0,
             color=CR, ecolor=CR, elinewidth=0.8, capsize=2.0, zorder=3)
xs = np.linspace(d.bx.min() * 1.25, d.bx.max() * 1.25, 40)
axA.plot(xs, b_ivw * xs, color='0.2', lw=1.0, zorder=2)
axA.set_xlabel('VDR eQTL effect (per allele)', fontsize=7.0)
axA.set_ylabel('Knee OA log OR (per allele)', fontsize=7.0)
axA.tick_params(labelsize=6.2)
for sp in ('top', 'right'):
    axA.spines[sp].set_visible(False)
axA.text(0.04, 0.95, 'IVW slope %.4f' % b_ivw, transform=axA.transAxes, fontsize=5.6,
         va='top')

# B: forest
labs = list(d.SNP) + ['IVW']
ors = np.exp(list(wald) + [b_ivw])
los = np.exp(list(wald - 1.96 * wse) + [b_ivw - 1.96 * se_rand])
his = np.exp(list(wald + 1.96 * wse) + [b_ivw + 1.96 * se_rand])
y = np.arange(len(labs))
axB.axvline(1, color='0.45', ls=(0, (4, 3)), lw=0.7, zorder=1)
for i in range(len(labs)):
    c = '0.15' if i == len(labs) - 1 else CR
    axB.plot([los[i], his[i]], [y[i], y[i]], color=c, lw=1.0, zorder=3)
    axB.plot(ors[i], y[i], 'o', ms=4.4, color=c, zorder=4)
axB.set_yticks(y)
axB.set_yticklabels(labs, fontsize=6.0)
axB.invert_yaxis()
axB.set_xscale('log')
axB.set_xlabel('Causal OR per SD of expression (95% CI)', fontsize=7.0)
axB.tick_params(labelsize=6.2)
for sp in ('top', 'right'):
    axB.spines[sp].set_visible(False)

# C: leave-one-out
loo = []
for i in range(len(d)):
    m = np.ones(len(d), bool); m[i] = False
    wi = 1 / wse[m] ** 2
    bi = (wi * wald[m]).sum() / wi.sum()
    si = np.sqrt(1 / wi.sum())
    qi = (wi * (wald[m] - bi) ** 2).sum()
    ph = max(1.0, qi / max(int(m.sum()) - 1, 1))
    loo.append((d.SNP.iloc[i], bi, si * np.sqrt(ph)))
labsC = ['omit ' + x[0] for x in loo] + ['all']
orsC = [np.exp(x[1]) for x in loo] + [np.exp(b_ivw)]
losC = [np.exp(x[1] - 1.96 * x[2]) for x in loo] + [np.exp(b_ivw - 1.96 * se_rand)]
hisC = [np.exp(x[1] + 1.96 * x[2]) for x in loo] + [np.exp(b_ivw + 1.96 * se_rand)]
yC = np.arange(len(labsC))
axC.axvline(1, color='0.45', ls=(0, (4, 3)), lw=0.7, zorder=1)
for i in range(len(labsC)):
    c = '0.15' if i == len(labsC) - 1 else '#2166AC'
    axC.plot([losC[i], hisC[i]], [yC[i], yC[i]], color=c, lw=1.0, zorder=3)
    axC.plot(orsC[i], yC[i], 'o', ms=4.4, color=c, zorder=4)
axC.set_yticks(yC)
axC.set_yticklabels(labsC, fontsize=6.0)
axC.invert_yaxis()
axC.set_xlabel('IVW OR after excluding one variant', fontsize=7.0)
axC.tick_params(labelsize=6.2)
for sp in ('top', 'right'):
    axC.spines[sp].set_visible(False)

for letter, ax, dx in (('A', axA, -0.20), ('B', axB, -0.30), ('C', axC, -0.34)):
    ax.text(dx, 1.04, letter, transform=ax.transAxes, fontsize=11, fontweight='bold',
            ha='left', va='bottom')

png = os.path.join(OUT, 'Figure_2.png')
tif = os.path.join(OUT, 'Figure_2.tif')
fig.savefig(png, dpi=300)
fig.savefig(tif, dpi=300)
plt.close(fig)

from PIL import Image
im = Image.open(tif)
if im.mode != 'RGB':
    rgb = im.convert('RGB'); im.close(); rgb.save(tif, dpi=(300, 300), compression='tiff_lzw'); rgb.close()
print('written:', tif, Image.open(tif).size)
