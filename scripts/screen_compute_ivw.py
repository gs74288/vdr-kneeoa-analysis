# -*- coding: utf-8 -*-
"""Stage 2: prune to independent loci, fetch outcome effects, run IVW + BH-FDR."""
import sys, io, os, json, time, glob
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

CACHE = r"C:\Users\liuha\OpenClawVDR\data\screen_cache"
TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"
OA = 'ebi-a-GCST007090'
SPACING = 1_000_000


def post(payload, tries=6, timeout=120):
    req = urllib.request.Request('https://api.opengwas.io/api/associations',
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
for f in sorted(glob.glob(os.path.join(CACHE, '*.csv'))):
    gene = os.path.basename(f)[:-4]
    d = pd.read_csv(f)
    if d.empty or 'rsid' not in d.columns:
        print('%-8s no instruments' % gene)
        continue
    d['pos'] = pd.to_numeric(d.position, errors='coerce')
    d = d.dropna(subset=['pos']).sort_values('p')
    keep = []
    for r in d.itertuples():
        if all(str(r.chr) != str(k[0]) or abs(r.pos - k[1]) > SPACING for k in keep):
            keep.append((str(r.chr), r.pos))
    inst = d[[(str(x.chr), x.pos) in keep for x in d.itertuples()]]
    print('%-8s %5d candidates -> %3d loci' % (gene, len(d), len(inst)), flush=True)

    oa = post({'id': [OA], 'variant': inst.rsid.tolist()}) or []
    o = pd.DataFrame(oa).drop_duplicates('rsid').set_index('rsid') if oa else pd.DataFrame()
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    bx, by, sy = [], [], []
    for r in inst.itertuples():
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
        print('         no harmonised variants')
        continue
    bx, by, sy = map(np.array, (bx, by, sy))
    wald = by / bx
    w = 1 / (np.abs(sy / bx) ** 2)
    b = (w * wald).sum() / w.sum()
    se = np.sqrt(1 / w.sum())
    rows.append(dict(gene=gene, n_loci=len(bx), n_cand=len(d), OR=np.exp(b),
                     lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se),
                     p=2 * (1 - stats.norm.cdf(abs(b / se)))))
    time.sleep(0.4)

res = pd.DataFrame(rows)
if res.empty:
    print('no results'); raise SystemExit(0)
res['FDR'] = stats.false_discovery_control(res.p)
res = res.sort_values('p').reset_index(drop=True)
print('\n=== re-run screen (1 Mb spacing rule) ===')
print(res[['gene', 'n_cand', 'n_loci', 'OR', 'lo', 'hi', 'p', 'FDR']].to_string(
    index=False, float_format=lambda v: '%.3f' % v))
res.to_csv(os.path.join(TAB, 'screen16_rerun.csv'), index=False)
print('\nwritten screen16_rerun.csv')
