# Exploratory analyses not reported in the manuscript

This folder holds analyses that were run while preparing the manuscript but are
**not reported in the paper or its supplement**. They are kept for transparency,
not as supporting evidence. Nothing here should be cited as a result of the study.

## `mechano_modules_*.py` — mechanobiology module scoring

State-level scores for integrin/FAK, ROCK, YAP-TAZ and actomyosin gene modules in
two chondrocyte atlases, correlated against the scDRS disease-relevance score.

| Atlas | Peak state | Correlation with scDRS |
|---|---|---|
| GSE324993 | RepC (+2.04) | rho = +0.63, P = 1.6e-6 |
| GSE104782 | RegC (+1.76) | rho = +0.012, P = 0.93 |

The association did not reproduce across atlases, so it was excluded from the
manuscript. Tables: `tables/GSE*_mechano_scores.csv`.

## Geneformer in-silico perturbation (not completed)

Attempted zero-shot in-silico deletion/overexpression of VDR in chondrocytes with
Geneformer V2-104M, as an alternative to wet-lab perturbation.

Environment (Windows, CPU + RTX 5060):

- geneformer source from the HuggingFace model repository (the GitHub project is gone)
- `transformers==4.46` is required; 5.x removes `SpecialTokensMixin`
- `datasets==3.2.0`, `torch==2.11.0+cu128` (Blackwell GPU needs cu128 wheels)
- `tdigest` cannot be installed on Windows (its `accumulation-tree` dependency needs a
  C++ toolchain); `tdigest_shim.py` substitutes an exact numpy-quantile implementation
- the package hard-codes `device="cuda"`; `geneformer_patch_device.py` makes device
  selection adaptive and `geneformer_repair_patch.py` repairs the patched modules
- `forward_batch_size=4` is required on an 8 GB card; 16 saturates VRAM and runs ~10x slower

Why the analysis was stopped: the zero-shot `cls` embeddings did not separate
chondrocyte states, and did worse than a trivial baseline.

| Dataset | Cells | States | kNN accuracy (5-fold) | Silhouette | Chance |
|---|---|---|---|---|---|
| GSE324993 | 1,000 | 8 | 0.267 | -0.026 | 0.125 |
| GSE104782 | 882 | 7 | 0.521 | 0.002 | 0.143 |
| GSE104782, PCA50 on log-normalised HVGs | 882 | 7 | **0.721** | **0.093** | 0.143 |

Geneformer's `cls` embeddings are not intended for zero-shot cell-type
classification (fine-tuning is required), so this comparison is not a fair test of
the model in general. It does mean the embeddings could not be shown to represent
chondrocyte state here, and there is no chondrocyte VDR knockdown experiment to
validate a predicted perturbation shift against. The pipeline runs end to end and
can be reused with a fine-tuned model or a matched ground-truth signature.

Run order: `geneformer_pilot_GSE324993.py` → `geneformer_tokenize.py` →
`geneformer_embed.py` → `geneformer_quality_check.py` → `geneformer_vs_pca_baseline.py`
(the `GSE104782_*` scripts repeat the same steps on the second atlas).
