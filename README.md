# Hummell-2026

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23824188.svg)](https://doi.org/10.5281/zenodo.23824188)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Statistical analysis code for:

> Hummell, N. A., Jonnagaddala, H., Kalocsay, M. & Moore, T. I. Receptor-proximal phosphoform organization precedes β2 integrin activation. *bioRxiv* (2026). doi: [PREPRINT DOI]

## Overview

This repository contains the statistical engine used to map the cytoplasmic signaling territory of the β2 integrin CD18 (ITGB2) by APEX2 proximity labeling in Jurkat T cells. CD18-APEX2 is compared against a compartment-matched membrane reference (CAAX-APEX2) under three conditions: resting (BSA), immobilized ICAM-1, and immobilized ICAM-1 with phorbol ester (PMA). The proteome and phosphoproteome were measured from the same labeling reactions in a single 18-plex TMTpro experiment with three independent biological replicates.

The code implements:

- Median loading normalization of summed reporter-ion signal-to-noise
- A day-blocked cell-means design with per-channel quality weights
- limma-trend empirical Bayes moderated t statistics with Benjamini-Hochberg correction
- Equivalence testing (two one-sided tests) for compositional stability

The engine is a Python re-implementation of the limma workflow ([Smyth 2004](https://doi.org/10.2202/1544-6115.1027); [Ritchie et al. 2015](https://doi.org/10.1093/nar/gkv007)) and does not require R or Bioconductor.

## Files

| File | Contents |
|---|---|
| `limma_engine.py` | Median loading normalization; cell-means design with preparation-day blocking; arrayWeights-style channel-quality weights; weighted least squares; intensity-dependent variance trend with empirical Bayes moderated t statistics; Benjamini-Hochberg adjustment |
| `tost_equivalence.py` | Two one-sided tests (TOST) against a symmetric margin; minimum detectable effect; three-way classification (different, equivalent, inconclusive); margin sweep; smallest equivalence margin |

## Requirements

Python ≥ 3.10 with NumPy, SciPy, statsmodels and openpyxl (see `requirements.txt`). Tested with Python 3.13, NumPy 2.5, SciPy 1.18, statsmodels 0.15 and openpyxl 3.1.

```
pip install -r requirements.txt
```

## Usage

```python
import numpy as np
from limma_engine import (normalize_log2, keep_mask_all, cellmeans_design,
                          array_weights, fit, STANDARD_CONTRASTS)
from tost_equivalence import classify

# M: features x 18 matrix of summed reporter signal-to-noise values, channels in CHANS order
L, factors = normalize_log2(M)
intensity = np.nanmean(np.log2(np.where(M > 0, M, np.nan)), axis=1)
keep = keep_mask_all(M)
X, cols = cellmeans_design(block_reps=True)
w = array_weights(L, X, keep)
res = fit(L, X, cols, STANDARD_CONTRASTS(cols), keep, weights=w, intensity=intensity, trend=True)

R = res["CD18_PMA_vs_BSA"]
df_total = res["_meta"]["dfres"] + res["_meta"]["d0"]
labels, q_equivalence = classify(R["logFC"], R["se"], R["q"], df_total)
```

Loading factors are calculated from all positively quantified features. When refitting filtered quantification tables, apply the published loading factors rather than recomputing them on the filtered subset.

Contrast names: `enrich_<condition>` is CD18-APEX2 minus CAAX-APEX2 within a condition; `CD18_<a>_vs_<b>` and `CAAX_<a>_vs_<b>` are within-arm activation contrasts. `PMA` denotes ICAM-1 plus PMA.

## Data

Raw mass spectrometry data are deposited in PRIDE (PXD085432) and will be made public upon publication. Processed quantification tables will be released with the published article.

## Citation

If you use this code, please cite the article above and this archive:

> Hummell, N. A., Jonnagaddala, H., Kalocsay, M. & Moore, T. I. (2026). Hummell-2026: analysis code for "Receptor-proximal phosphoform organization precedes β2 integrin activation" (Version 1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.23824188

## License

Code is released under the [MIT License](LICENSE).

## Contact

Travis I. Moore, PhD
Department of Integrative Biology and Pharmacology, McGovern Medical School
The University of Texas Health Science Center at Houston (UTHealth Houston)
travis.i.moore@uth.tmc.edu
