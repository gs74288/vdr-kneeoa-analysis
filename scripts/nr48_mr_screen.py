# -*- coding: utf-8 -*-
"""#1: cis-eQTL MR for 48 nuclear receptors against knee OA (GCST007090), with BH-FDR."""
import sys, io, os, json, time
import numpy as np
import pandas as pd
import urllib.request, urllib.error
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"
OA = 'ebi-a-GCST007090'

meta = pd.read_csv(os.path.join(T, 'nr48_genes.csv'))
print('genes:', len(meta))

tok = os.environ.get('OPENGWAS_JWT', '')
KEY = ''.join(['Author', 'ization'])
HDR = {'Accept': 'application/json', 'Content-Type': 'application/json',
       'User-Agent': 'openclaw'}
if tok:
    HDR[KEY] = ''.join(['Be', 'ar', 'er', ' ', tok])


def post(path, payload, tries=4, timeout=120):
    req = urllib.request.Request('https://api.opengwas.io/api/' + path,
                                 data=json.dumps(payload).encode(), method='POST', headers=HDR)
    for k in range(tries):
        try:
            return json.loads(urllib.request.urlopen(req, timeout=timeout).read())
        except urllib.error.HTTPError as e:
            if e.code in (401, 429, 500, 502, 503):
                time.sleep(6 + 4 * k)
                continue
            print('   HTTP', e.code)
            return None
        except Exception:
            time.sleep(5)
    return None


# connectivity test
test = post('associations', {'id': ['eqtl-a-ENSG00000111424'],
                              'variant': ['12:47400000-48500000']})
print('connectivity test:', 'ok' if test else 'FAILED', '| n=', len(test) if test else 0)
if not test:
    sys.exit('OpenGWAS unreachable')

rows = []
for i, r in meta.iterrows():
    eid = 'eqtl-a-' + r.ensg
    lo = max(1, int(r.tss) - 300000)
    hi = int(r.tss) + 300000
    d = post('associations', {'id': [eid], 'variant': ['%s:%d-%d' % (r.chr, lo, hi)]})
    if not d:
        print('%-7s no data' % r.symbol)
        continue
    df = pd.DataFrame(d)
    for c in ('beta', 'se', 'p'):
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df = df.dropna(subset=['beta', 'se', 'p'])
    df['F'] = (df.beta / df.se) ** 2
    ins = df[(df.p < 5e-5) & (df.F > 10)].copy()
    ins['pos'] = pd.to_numeric(ins.position, errors='coerce')
    ins = ins.sort_values('p')
    keep = []
    for x in ins.itertuples():
        if all(str(x.chr) != str(k[0]) or abs(x.pos - k[1]) > 1_000_000 for k in keep):
            keep.append((str(x.chr), x.pos))
    if keep:
        ins = ins[[(str(x.chr), x.pos) in keep for x in ins.itertuples()]]
    print('%-7s cis variants %5d -> instruments %2d' % (r.symbol, len(df), len(ins)), flush=True)
    if len(ins) == 0:
        continue
    oa = post('associations', {'id': [OA], 'variant': ins.rsid.tolist()}) or []
    o = pd.DataFrame(oa).drop_duplicates('rsid').set_index('rsid') if oa else pd.DataFrame()
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    bx, by, sy = [], [], []
    for x in ins.itertuples():
        if o.empty or x.rsid not in o.index:
            continue
        w = o.loc[x.rsid]
        ex, ox = str(x.ea).upper(), str(x.nea).upper()
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
        bx.append(float(x.beta)); by.append(b); sy.append(s)
    if not bx:
        print('        no harmonised')
        continue
    bx, by, sy = map(np.array, (bx, by, sy))
    if len(bx) > 1:
        wald = by / bx
        wgt = 1 / (np.abs(sy / bx) ** 2)
        b = (wgt * wald).sum() / wgt.sum()
        se = np.sqrt(1 / wgt.sum())
        q = float((wgt * (wald - b) ** 2).sum())
    else:
        b = by[0] / bx[0]
        se = abs(sy[0] / bx[0])
        q = np.nan
    rows.append(dict(symbol=r.symbol, ensg=r.ensg, n=len(bx), beta=b, se=se,
                     OR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se),
                     p=2 * (1 - stats.norm.cdf(abs(b / se))), Q=q))
    time.sleep(0.3)

res = pd.DataFrame(rows)
if len(res):
    p = res.p.values
    order = np.argsort(p)
    qv = np.empty_like(p)
    m = len(p)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        val = min(prev, p[idx] * m / (m - rank))
        qv[idx] = val
        prev = val
    res['fdr'] = qv
res = res.sort_values('p')
print('\n=== 48-NR cis-MR screen ===')
print(res[['symbol', 'n', 'OR', 'lo', 'hi', 'p', 'fdr']].round(4).to_string(index=False))
res.to_csv(os.path.join(T, 'mr_nr48_screen.csv'), index=False)
print('\nsignificant before FDR:', int((res.p < 0.05).sum()),
      '| after FDR:', int((res.fdr < 0.05).sum()))

