# -*- coding: utf-8 -*-
"""Windows shim for the `tdigest` package (accumulation-tree is not buildable here).

Provides a TDigest class backed by exact numpy quantiles; the interface used by
geneformer.emb_extractor is update() / percentile().
"""
import io
import os
import sys

target = None
for p in sys.path:
    if p.endswith('site-packages'):
        target = p
        break
if target is None:
    import site
    target = site.getsitepackages()[0]

path = os.path.join(target, 'tdigest.py')
src = '''# -*- coding: utf-8 -*-
"""Minimal pure-Python stand-in for tdigest.TDigest (exact quantiles via numpy)."""
import numpy as np


class TDigest:
    def __init__(self, delta=0.01, K=25):
        self.delta = delta
        self.K = K
        self._values = []

    def update(self, x, w=1):
        try:
            self._values.append(float(x))
        except (TypeError, ValueError):
            pass

    def batch_update(self, xs, ws=None):
        for x in xs:
            self.update(x)

    def percentile(self, q):
        if not self._values:
            return float('nan')
        return float(np.percentile(self._values, q * 100.0))

    def quantile(self, q):
        return self.percentile(q)

    def cdf(self, x):
        if not self._values:
            return 0.0
        return float((np.array(self._values) <= x).mean())

    def trimed_mean(self):
        return float(np.mean(self._values)) if self._values else float('nan')

    def to_dict(self):
        return {'values': list(self._values)}
'''
with open(path, 'w', encoding='utf-8', newline='') as fh:
    fh.write(src)
print('shim written to', path)
