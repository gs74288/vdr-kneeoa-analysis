# -*- coding: utf-8 -*-
"""Summarise the GCTA-COJO conditional results."""
import sys, io, os
import pandas as pd
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
OUT = r"C:\Users\liuha\OpenClawVDR\tools\cojo"
d = pd.read_csv(os.path.join(OUT, 'cond_result.cma.cojo'), sep='\t')
d.columns = [c.strip() for c in d.columns]
print('columns:', list(d.columns))
print('rows:', len(d))

VDR = (47_841_537 + 393_783, 47_943_048 + 393_783)
COL = (47_972_967 + 393_783, 48_004_554 + 393_783)

d['bp'] = pd.to_numeric(d.bp, errors='coerce')
d['p'] = pd.to_numeric(d.p, errors='coerce')
d['pC'] = pd.to_numeric(d.pC, errors='coerce')

print('\n=== unconditional vs conditional ===')
print('  all:      min p = %.2e   min pC = %.2e' % (d.p.min(), d.pC.min()))
for lab, (lo, hi) in (('VDR body', VDR), ('COL2A1 body', COL)):
    s = d[(d.bp >= lo) & (d.bp <= hi)]
    if s.empty:
        print('  %-12s no variants' % lab); continue
    print('  %-12s n=%3d  min p = %.2e  min pC = %.2e  |  pC<5e-8: %d'
          % (lab, len(s), s.p.min(), s.pC.min(), int((s.pC < 5e-8).sum())))

print('\n=== strongest residual signals after conditioning ===')
top = d.nsmallest(8, 'pC')[['SNP', 'bp', 'refA', 'b', 'p', 'bC', 'pC']]
print(top.to_string(index=False))

# how many genome-wide significant variants remain?
print('\nvariants with pC < 5e-8: %d (of %d)' % (int((d.pC < 5e-8).sum()), len(d)))
print('variants with p  < 5e-8: %d' % int((d.p < 5e-8).sum()))
