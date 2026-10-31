# Hummell-2026

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23824188.svg)](https://doi.org/10.5281/zenodo.23824188)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Analysis code and processed data for:

> Hummell, N. A., Jonnagaddala, H., Kalocsay, M. & Moore, T. I. Receptor-proximal phosphoform organization precedes β2 integrin activation. *bioRxiv* (2026). doi: [PREPRINT DOI]

## Overview

This repository contains the statistical engine and processed quantification tables used to map the cytoplasmic signaling territory of the β2 integrin CD18 (ITGB2) by APEX2 proximity labeling in Jurkat T cells. CD18-APEX2 is compared against a compartment-matched membrane reference (CAAX-APEX2) under three conditions: resting (BSA), immobilized ICAM-1, and immobilized ICAM-1 with phorbol ester (PMA). The proteome and phosphoproteome were measured from the same labeling reactions in a single 18-plex TMTpro experiment with three independent biological replicates.

The code implements:

- Median loading normalization of summed reporter-ion signal-to-noise
- A day-blocked cell-means design with per-channel quality weights
- limma-trend empirical Bayes moderated t statistics with Benjamini-Hochberg correction
- Equivalence testing (two one-sided tests) for compositional stability
- The phosphosite selectivity index (PSI), the log2 CD18/CAAX enrichment of a phosphosite minus that of its parent protein
- Difference-of-differences contrasts for receptor-preferential activation responses

The engine is a Python re-implementation of the limma workflow ([Smyth 2004](https://doi.org/10.2202/1544-6115.1027); [Ritchie et al. 2015](https://doi.org/10.1093/nar/gkv007)) and does not require R or Bioconductor.

## Citation

If you use this code or data, please cite the article above and this archive:

> Hummell, N. A., Jonnagaddala, H., Kalocsay, M. & Moore, T. I. (2026). Hummell-2026: analysis code for "Receptor-proximal phosphoform organization precedes β2 integrin activation" (Version 1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.23824188

## License

- **Code** is released under the [MIT License](LICENSE).
- **Processed data** in `data/` are released under [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
- **Raw data** in PRIDE (PXD085432) are subject to the PRIDE repository terms.

## Contact

Travis I. Moore, PhD
Department of Integrative Biology and Pharmacology, McGovern Medical School
The University of Texas Health Science Center at Houston (UTHealth Houston)
travis.i.moore@uth.tmc.edu
