# -*- coding: utf-8 -*-
"""Repair the broken DEVICE insertion in geneformer modules."""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

SRC = r"C:\Users\liuha\OpenClawVDR\tools\geneformer\src\geneformer"
TRIO = ['', 'import torch as _torch',
        'DEVICE = "cuda" if _torch.cuda.is_available() else "cpu"']
HEADER = ['import torch as _torch',
          'DEVICE = "cuda" if _torch.cuda.is_available() else "cpu"', '']

for f in sorted(os.listdir(SRC)):
    if not f.endswith('.py'):
        continue
    p = os.path.join(SRC, f)
    s = open(p, encoding='utf-8').read()
    if 'DEVICE' not in s:
        continue
    lines = s.split('\n')
    # remove any occurrence of the inserted trio (in any position)
    out, i = [], 0
    while i < len(lines):
        if lines[i:i + 3] == TRIO:
            i += 3
            continue
        out.append(lines[i])
        i += 1
    # drop a stray duplicated definition if present
    out = [ln for ln in out if ln != 'DEVICE = "cuda" if _torch.cuda.is_available() else "cpu"']
    out = [ln for ln in out if ln != 'import torch as _torch']
    # ensure the header is present at the very top
    while out and not out[0].strip():
        out.pop(0)
    out = HEADER + out
    open(p, 'w', encoding='utf-8', newline='').write('\n'.join(out))
    print('repaired', f)

# sanity: compile every module
import py_compile
bad = []
for f in sorted(os.listdir(SRC)):
    if not f.endswith('.py'):
        continue
    try:
        py_compile.compile(os.path.join(SRC, f), doraise=True)
    except Exception as e:
        bad.append((f, str(e)[:120]))
print('\ncompile errors:', bad if bad else 'none')
