#!/usr/bin/env python3
"""
tost_equivalence: equivalence testing for compositional stability
=================================================================
Two one-sided tests (TOST) against a prespecified symmetric margin, as used in
Hummell et al. A non-significant differential test is absence of evidence; TOST
gives positive evidence that a change lies within +/- delta.

Note: for a limma-trend fit, df_total = residual df + prior df (d0).
"""
import numpy as np
from scipy import stats

from limma_engine import _bh


def tost(logFC, se, df_total, delta=0.5):
    """Two one-sided tests for equivalence within +/- delta.

    Returns the TOST p-value: max of the two one-sided p-values. Small p means
    the effect is significantly *inside* the bound, i.e. statistically
    equivalent to no change.
    """
    logFC = np.asarray(logFC, float)
    se = np.asarray(se, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p_lower = stats.t.sf((logFC + delta) / se, df=df_total)   # H0: beta <= -delta
        p_upper = stats.t.cdf((logFC - delta) / se, df=df_total)  # H0: beta >= +delta
    return np.maximum(p_lower, p_upper)


def min_detectable_effect(se, power=0.8, alpha=0.05):
    """Approximate two-sided MDE. ~2.8*se at 80% power, alpha 0.05."""
    z = stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power)
    return z * np.asarray(se, float)


def classify(logFC, se, q_diff, df_total, delta=0.5, qcut=0.05, effect=0.5):
    """Three-way classification used in the paper.
    'different'    : differential q < qcut and |logFC| >= effect
    'equivalent'   : BH-adjusted TOST q < qcut at +/- delta
    'inconclusive' : neither
    Returns (labels array, BH-adjusted TOST q)."""
    logFC = np.asarray(logFC, float)
    q_eq = _bh(tost(logFC, se, df_total, delta))
    diff = (np.asarray(q_diff, float) < qcut) & (np.abs(logFC) >= effect)
    lab = np.where(diff, "different", np.where(q_eq < qcut, "equivalent", "inconclusive"))
    return lab, q_eq


def margin_sweep(logFC, se, df_total, deltas=np.arange(0.1, 2.01, 0.05), qcut=0.05):
    """Fraction of features declared equivalent (BH-adjusted TOST q < qcut)
    across a range of symmetric margins."""
    return {round(float(d), 3): float(np.mean(_bh(tost(logFC, se, df_total, d)) < qcut))
            for d in deltas}


def smallest_margin(logFC, se, df_total, alpha=0.05):
    """Smallest symmetric margin at which the unadjusted TOST is significant at
    alpha: |logFC| + t(1 - alpha, df) * se (Supplementary Fig. 8b)."""
    return np.abs(np.asarray(logFC, float)) + stats.t.ppf(1 - alpha, df_total) * np.asarray(se, float)
