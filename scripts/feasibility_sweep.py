# -*- coding: utf-8 -*-
"""Feasibility sweep: decoupleR availability, cartilage ATAC/Hi-C data, Genebass phenotypes."""
import sys, io, os, json, importlib, subprocess, urllib.request, urllib.parse, urllib.error

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
UA = {'User-Agent': 'Mozilla/5.0 (research script)'}

print('=== A. local packages ===')
for m in ('decoupler', 'pyscenic', 'arboreto', 'requests', 'gseapy'):
    try:
        mod = importlib.import_module(m)
        print('  %-10s OK %s' % (m, getattr(mod, '__version__', '')))
    except Exception as e:
        print('  %-10s missing (%s)' % (m, type(e).__name__))

print('\n=== B. pip availability (decoupler) ===')
try:
    r = subprocess.run([sys.executable, '-m', 'pip', 'index', 'versions', 'decoupler'],
                       capture_output=True, text=True, timeout=90)
    out = (r.stdout or r.stderr).strip().splitlines()
    print('  rc=%d %s' % (r.returncode, out[:3]))
except Exception as e:
    print('  pip index failed:', e)


def eutils(term, retmax=6):
    q = urllib.parse.urlencode({'db': 'gds', 'term': term, 'retmax': retmax,
                                'retmode': 'json', 'sort': 'relevance'})
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + q
    try:
        d = json.loads(OPENER.open(urllib.request.Request(url, headers=UA), timeout=45).read())
        ids = d['esearchresult'].get('idlist', [])
        total = d['esearchresult'].get('count', '0')
    except Exception as e:
        print('  query failed:', e)
        return
    print('  total=%s' % total)
    if not ids:
        return
    q2 = urllib.parse.urlencode({'db': 'gds', 'id': ','.join(ids), 'retmode': 'json'})
    url2 = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?' + q2
    try:
        s = json.loads(OPENER.open(urllib.request.Request(url2, headers=UA), timeout=45).read())
        for uid in s['result']['uids']:
            rec = s['result'][uid]
            print('    %-10s %-6s n=%-6s %s' % (rec.get('accession'), rec.get('gdstype'),
                                                rec.get('n_samples'), rec.get('title')[:95]))
    except Exception as e:
        print('  summary failed:', e)


print('\n=== C. NCBI GEO: cartilage ATAC ===')
eutils('(cartilage[Title] OR chondrocyte[Title]) AND (ATAC[Title/Abstract] OR '
       '"chromatin accessibility"[Title/Abstract]) AND "Homo sapiens"[Organism]')

print('\n=== C2. NCBI GEO: cartilage Hi-C / PCHi-C ===')
eutils('(cartilage[Title] OR chondrocyte[Title]) AND ("Hi-C"[Title/Abstract] OR '
       '"capture Hi-C"[Title/Abstract] OR "PCHi-C"[Title/Abstract]) AND "Homo sapiens"[Organism]')

print('\n=== C3. NCBI GEO: cartilage H3K27ac ===')
eutils('(cartilage[Title] OR chondrocyte[Title]) AND ("H3K27ac"[Title/Abstract] OR '
       '"H3K27ac"[All Fields]) AND "Homo sapiens"[Organism]')

print('\n=== D. Genebass endpoints ===')
for url in ('https://genebass.org/api/phenotypes',
            'https://app.genebass.org/api/phenotypes',
            'https://storage.googleapis.com/genebass-results/manifest.csv',
            'https://genebass.org/'):
    try:
        r = OPENER.open(urllib.request.Request(url, headers=UA), timeout=25)
        body = r.read(400)
        print('  %-58s %s | %r' % (url, r.status, body[:120]))
    except urllib.error.HTTPError as e:
        print('  %-58s HTTP %s' % (url, e.code))
    except Exception as e:
        print('  %-58s %s' % (url, type(e).__name__))
