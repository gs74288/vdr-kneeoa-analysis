# -*- coding: utf-8 -*-
"""Patch geneformer to select CPU when CUDA is unavailable."""
import sys, io, os, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

SRC = r"C:\Users\liuha\OpenClawVDR\tools\geneformer\src\geneformer"
files = [f for f in os.listdir(SRC) if f.endswith('.py')]
changed = []
for f in files:
    p = os.path.join(SRC, f)
    s = open(p, encoding='utf-8').read()
    if '"cuda"' not in s:
        continue
    n = s.count('device="cuda"') + s.count("device='cuda'") + s.count('.to("cuda")') + s.count("map_location=\"cuda\"")
    s = s.replace('device="cuda"', 'device=DEVICE')
    s = s.replace("device='cuda'", 'device=DEVICE')
    s = s.replace('.to("cuda")', '.to(DEVICE)')
    s = s.replace('map_location="cuda"', 'map_location=DEVICE')
    if 'DEVICE = ' not in s:
        # insert after imports
        lines = s.split('\n')
        idx = max(i for i, ln in enumerate(lines[:40]) if ln.startswith(('import ', 'from ')))
        lines.insert(idx + 1, '')
        lines.insert(idx + 2, 'import torch as _torch')
        lines.insert(idx + 3, 'DEVICE = "cuda" if _torch.cuda.is_available() else "cpu"')
        s = '\n'.join(lines)
    open(p, 'w', encoding='utf-8', newline='').write(s)
    changed.append((f, n))

print('patched files:')
for f, n in changed:
    print('  %-32s %d replacements' % (f, n))
