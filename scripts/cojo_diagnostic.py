# -*- coding: utf-8 -*-
"""Recompute conditional statistics with a variance floor and stratified reporting."""
import sys, io, os
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
OUT = r"C:\Users\liuha\OpenClawVDR\repo\tables\cojo"
d = pd.read_csv(os.path.join(OUT, 'conditional_oa.csv'))

FLOOR = 1e-3
r = d.r_to_cond.values
r2 = np.clip(1 - r ** 2, FLOOR, None)
se_cond = d.se2_orig.values * np.sqrt(r2)
b_cond = d.b2_orig.values - (r * d.se2_orig.values / d.se2_orig.values[84]) * d.b2_orig.values[84]
z = b_cond / se_cond
p = 2 * (1 - stats.norm.cdf(np.abs(z)))
d['b2_cond'] = b_cond
d['se2_cond'] = se_cond
d['p2_cond'] = p
d['abs_r'] = np.abs(r)
d.to_csv(os.path.join(OUT, 'conditional_oa_fixed.csv'), index=False)

print('=== conditioning on rs11168351 ===')
for lab, m in (('all variants', np.ones(len(d), bool)),
               ('|r| < 0.8', d.abs_r < 0.8),
               ('|r| < 0.5', d.abs_r < 0.5),
               ('VDR gene body', (d.pos >= 47_841_537) & (d.pos <= 47_943_048))):
    s = d[m]
    if s.empty:
        continue
    print('  %-14s n=%4d  min conditional p = %.2e  (orig %.2e)'
          % (lab, len(s), s.p2_cond.min(), s.p2_orig.min()))

print('\nstrongest residual signals (|r|<0.8):')
s = d[d.abs_r < 0.8].nsmallest(6, 'p2_cond')
print(s[['snp', 'pos', 'p2_orig', 'p2_cond', 'abs_r']].to_string(index=False))

print('\nVDR-body variants after conditioning:')
v = d[(d.pos >= 47_841_537) & (d.pos <= 47_943_048)].nsmallest(5, 'p2_cond')
print(v[['snp', 'pos', 'p2_orig', 'p2_cond', 'abs_r']].to_string(index=False))
