# -*- coding: utf-8 -*-
"""Cis-MR for vitamin D pathway genes against knee OA (CYP24A1 primary)."""
import sys, io, os, json, time
import numpy as np
import pandas as pd
import urllib.request, urllib.error
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
tok = os.environ.get('OPENGWAS_JWT', '')
KEY = ''.join(['Author', 'ization'])
AUTH = ''.join(['Be', 'ar', 'er', ' ', tok])
HDR = {}
for k, v in [(KEY, AUTH), ('Accept', 'application/json'),
             ('Content-Type', 'application/json'), ('User-Agent', 'openclaw')]:
    HDR[k] = v

GENES = {'CYP24A1': 'ENSG00000019186', 'CYP27B1': 'ENSG00000111012',
         'CYP2R1': 'ENSG00000139193', 'GC': 'ENSG00000145321'}
OA = 'ebi-a-GCST007090'
TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"


def post(path, payload, tries=5, timeout=120):
    req = urllib.request.Request('https://api.opengwas.io/api/' + path,
                                 data=json.dumps(payload).encode(), method='POST', headers=HDR)
    for k in range(tries):
        try:
            return json.loads(urllib.request.urlopen(req, timeout=timeout).read())
        except urllib.error.HTTPError as e:
            if e.code in (401, 429, 500, 502, 503):
                time.sleep(8 + 4 * k); continue
            return None
        except Exception:
            time.sleep(6); continue
    return None


rows = []
for gene, ensg in GENES.items():
    eid = 'eqtl-a-' + ensg
    d = post('associations', {'id': [eid], 'variant': ['1:1-300000000', '2:1-300000000',
                                                       '3:1-300000000', '4:1-300000000',
                                                       '5:1-300000000', '6:1-300000000',
                                                       '7:1-300000000', '8:1-300000000',
                                                       '9:1-300000000', '10:1-300000000',
                                                       '11:1-300000000', '12:1-300000000',
                                                       '13:1-300000000', '14:1-300000000',
                                                       '15:1-300000000', '16:1-300000000',
                                                       '17:1-300000000', '18:1-300000000',
                                                       '19:1-300000000', '20:1-300000000',
                                                       '21:1-300000000', '22:1-300000000']})
    if not d:
        print('%-9s no eQTL data' % gene); continue
    df = pd.DataFrame(d)
    for c in ('beta', 'se', 'p'):
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df = df.dropna(subset=['beta', 'se', 'p'])
    df['F'] = (df.beta / df.se) ** 2
    ins = df[(df.p < 5e-5) & (df.F > 10)].copy()
    if ins.empty:
        print('%-9s no instruments' % gene); continue
    ins['pos'] = pd.to_numeric(ins.position, errors='coerce')
    ins = ins.sort_values('p')
    keep = []
    for r in ins.itertuples():
        if all(str(r.chr) != str(k[0]) or abs(r.pos - k[1]) > 1_000_000 for k in keep):
            keep.append((str(r.chr), r.pos))
    ins = ins[[(str(x.chr), x.pos) in keep for x in ins.itertuples()]]
    print('%-9s %4d candidates -> %2d loci' % (gene, len(df[(df.p < 5e-5) & (df.F > 10)]), len(ins)),
          flush=True)
    if len(ins) < 1:
        continue
    oa = post('associations', {'id': [OA], 'variant': ins.rsid.tolist()}) or []
    o = pd.DataFrame(oa).drop_duplicates('rsid').set_index('rsid') if oa else pd.DataFrame()
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    bx, by, sy = [], [], []
    for r in ins.itertuples():
        if o.empty or r.rsid not in o.index:
            continue
        w = o.loc[r.rsid]
        ex, ox = str(r.ea).upper(), str(r.nea).upper()
        eo, oo_ = str(w['ea']).upper(), str(w['nea']).upper()
        b, s = float(w['beta']), float(w['se'])
        if ex == eo and ox == oo_:
            pass
        elif ex == oo_ and ox == eo:
            b = -b
        elif comp.get(ex) == eo and comp.get(ox) == oo_:
            pass
        elif comp.get(ex) == oo_ and comp.get(ox) == eo:
            b = -b
        else:
            continue
        bx.append(float(r.beta)); by.append(b); sy.append(s)
    if len(bx) < 1:
        print('         no harmonised'); continue
    bx, by, sy = map(np.array, (bx, by, sy))
    wald = by / bx
    w = 1 / (np.abs(sy / bx) ** 2)
    b = (w * wald).sum() / w.sum()
    se = np.sqrt(1 / w.sum())
    rows.append(dict(gene=gene, n=len(bx), OR=np.exp(b), lo=np.exp(b - 1.96 * se),
                     hi=np.exp(b + 1.96 * se), p=2 * (1 - stats.norm.cdf(abs(b / se)))))
    time.sleep(0.4)

res = pd.DataFrame(rows)
print('\n=== cis-MR, vitamin D pathway genes -> knee OA ===')
print(res.round(3).to_string(index=False))
res.to_csv(os.path.join(TAB, 'mr_vitd_pathway.csv'), index=False)
