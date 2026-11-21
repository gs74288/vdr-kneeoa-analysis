# -*- coding: utf-8 -*-
"""Figure 5 rebuilt from the raw GSE114007 count matrices."""
import sys, io, os, gzip
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7,
                     'axes.linewidth': 0.7, 'xtick.major.width': 0.7,
                     'ytick.major.width': 0.7, 'axes.labelsize': 7.5})

BULK = r"C:\Users\liuha\OpenClawVDR\bulk"
OUT = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\Figures"

AXIS = ["NLRP3", "PYCARD", "CASP1", "IL1B", "IL18", "GSDMD", "TXNIP", "AIM2", "NLRC4", "CASP8"]

oa = pd.read_csv(os.path.join(BULK, 'GSE114007_OA_counts.txt.gz'), sep='\t', index_col=0)
nm = pd.read_csv(os.path.join(BULK, 'GSE114007_normal_counts.txt.gz'), sep='\t', index_col=0)
print('OA', oa.shape, '| Normal', nm.shape)

X = pd.concat([oa, nm], axis=1)
X = X.groupby(level=0).mean()
groups = np.array(['OA'] * oa.shape[1] + ['Normal'] * nm.shape[1])
print('combined', X.shape, '| groups', np.unique(groups, return_counts=True))

vdr = X.loc['VDR'].values.astype(float)
g = [x for x in AXIS if x in X.index]
Z = X.loc[g].astype(float)
Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1, ddof=0).replace(0, 1), axis=0)
axis_score = Z.mean(axis=0).values
print('axis genes used:', g)

m = groups == 'OA'
p_vdr = stats.mannwhitneyu(vdr[m], vdr[~m]).pvalue
p_axis = stats.mannwhitneyu(axis_score[m], axis_score[~m]).pvalue
r, pr = stats.spearmanr(vdr[m], axis_score[m])
print('\nVDR   OA %.3f vs Normal %.3f  P=%.4g' % (vdr[m].mean(), vdr[~m].mean(), p_vdr))
print('AXIS  OA %.3f vs Normal %.3f  P=%.4g' % (axis_score[m].mean(), axis_score[~m].mean(), p_axis))
print('Spearman(VDR, axis) in OA: r=%.3f P=%.4g' % (r, pr))

fig = plt.figure(figsize=(6.70, 2.65))
axA = fig.add_axes([0.075, 0.185, 0.240, 0.680])
axB = fig.add_axes([0.395, 0.185, 0.240, 0.680])
axC = fig.add_axes([0.720, 0.185, 0.255, 0.680])


def violin(ax, data_by_group, ylab, tag_p):
    vals = [data_by_group['Normal'], data_by_group['OA']]
    parts = ax.violinplot(vals, positions=[1, 2], widths=0.72, showextrema=False)
    for pc, c in zip(parts['bodies'], ('#9E9E9E', '#B2182B')):
        pc.set_facecolor(c); pc.set_edgecolor('none'); pc.set_alpha(0.55)
    for i, v in enumerate(vals):
        xj = np.random.RandomState(1).normal(i + 1, 0.045, len(v))
        ax.scatter(xj, v, s=5, c=('#6E6E6E', '#8C1A22')[i], alpha=0.85, linewidths=0, zorder=3)
        ax.plot([i + 0.78, i + 1.22], [np.median(v)] * 2, color='k', lw=0.9, zorder=4)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Normal\n(n=%d)' % len(vals[0]), 'OA\n(n=%d)' % len(vals[1])], fontsize=6.2)
    ax.set_ylabel(ylab, fontsize=7.2)
    ax.set_xlim(0.55, 2.45)
    ax.tick_params(labelsize=6.4)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.text(0.97, 0.96, tag_p, transform=ax.transAxes, ha='right', va='top', fontsize=6.0)


def fmt_p(p):
    return 'P < 0.001' if p < 0.001 else 'P = %.3f' % p


violin(axA, {'Normal': vdr[~m], 'OA': vdr[m]}, 'VDR (log2 expression)',
       fmt_p(p_vdr))
violin(axB, {'Normal': axis_score[~m], 'OA': axis_score[m]}, 'NLRP3-axis score (mean z)',
       fmt_p(p_axis))

axC.scatter(vdr[m], axis_score[m], s=16, c='#B2182B', alpha=0.85, linewidths=0, zorder=3)
b1, b0 = np.polyfit(vdr[m], axis_score[m], 1)
xs = np.linspace(vdr[m].min(), vdr[m].max(), 50)
axC.plot(xs, b0 + b1 * xs, color='0.25', lw=0.8, zorder=2)
axC.set_xlabel('VDR (log2 expression)', fontsize=7.2)
axC.set_ylabel('NLRP3-axis score', fontsize=7.2)
axC.tick_params(labelsize=6.4)
for sp in ('top', 'right'):
    axC.spines[sp].set_visible(False)
axC.text(0.04, 0.96, 'r = %.2f\n%s' % (r, fmt_p(pr)), transform=axC.transAxes,
         ha='left', va='top', fontsize=6.0)

for letter, ax, dx in (('A', axA, -0.28), ('B', axB, -0.26), ('C', axC, -0.24)):
    ax.text(dx, 1.05, letter, transform=ax.transAxes, fontsize=11, fontweight='bold',
            ha='left', va='bottom')

os.makedirs(OUT, exist_ok=True)
png = os.path.join(OUT, 'Figure_5_new.png')
tif = os.path.join(OUT, 'Figure_5_new.tif')
fig.savefig(png, dpi=300)
fig.savefig(tif, dpi=300)
plt.close(fig)

from PIL import Image
im = Image.open(tif)
if im.mode != 'RGB':
    rgb = im.convert('RGB'); im.close(); rgb.save(tif, dpi=(300, 300), compression='tiff_lzw'); rgb.close()
print('\nwritten:', tif, Image.open(tif).size)
