#!/usr/bin/env python3
"""
limma_engine: moderated linear modeling for multiplexed TMT proximity proteomics
=================================================================================
Python implementation of the limma-style workflow used in Hummell et al.,
"Receptor-proximal phosphoform organization precedes β2 integrin activation".
No R or Bioconductor installation is required.

Workflow
  1. Quantity = summed reporter signal-to-noise (sn_sum) per channel.
  2. Median loading normalization across the 18 channels, then log2.
  3. Weighted least-squares fit of a cell-means design (six bait-by-condition
     groups) with preparation day as a blocking factor.
  4. Per-channel quality weights estimated iteratively (arrayWeights-style).
  5. Intensity-dependent mean-variance trend (limma-trend) and empirical-Bayes
     moderated t statistics; Benjamini-Hochberg adjustment within each contrast.

Channel layout (C = CAAX-APEX2 reference, KO = CD18-APEX2 in CD18-knockout cells;
replicate index = preparation day):
  C_BSA   126 127n 127c    C_ICAM1  128n 128c 129n    C_PMA  129c 130n 130c
  KO_BSA  131n 131c 132n   KO_ICAM1 132c 133n 133c    KO_PMA 134n 134c 135n
"PMA" denotes ICAM-1 plus PMA.

"""
import numpy as np
from openpyxl import load_workbook
from scipy.special import polygamma, digamma
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

CHANS = ['126','127n','127c','128n','128c','129n','129c','130n','130c',
         '131n','131c','132n','132c','133n','133c','134n','134c','135n']
IDX = {c:i for i,c in enumerate(CHANS)}
GROUPS = ['C_BSA','C_BSA','C_BSA','C_ICAM1','C_ICAM1','C_ICAM1','C_PMA','C_PMA','C_PMA',
          'KO_BSA','KO_BSA','KO_BSA','KO_ICAM1','KO_ICAM1','KO_ICAM1','KO_PMA','KO_PMA','KO_PMA']
GL  = ['C_BSA','C_ICAM1','C_PMA','KO_BSA','KO_ICAM1','KO_PMA']
REP = [0,1,2]*6

COND_CH = {
  'C_BSA':['126','127n','127c'],   'C_ICAM1':['128n','128c','129n'],  'C_PMA':['129c','130n','130c'],
  'KO_BSA':['131n','131c','132n'], 'KO_ICAM1':['132c','133n','133c'], 'KO_PMA':['134n','134c','135n'],
}

def load_sn_sum(path, prefix, label_cols):
    wb = load_workbook(path, read_only=True); ws = wb[wb.sheetnames[0]]
    it = ws.iter_rows(values_only=True); h = [str(x) if x is not None else "" for x in next(it)]
    ci = [h.index(f"{prefix}{c}_sn_sum") for c in CHANS]
    li = {k:h.index(v) for k,v in label_cols.items() if v in h}
    labs, rows = [], []
    for r in it:
        rows.append([np.nan if r[c] is None else float(r[c]) for c in ci])
        labs.append({k:r[c] for k,c in li.items()})
    wb.close()
    return labs, np.array(rows, float)

def normalize_log2(M):
    med = np.array([np.median(M[:,j][M[:,j]>0]) if np.any(M[:,j]>0) else np.nan for j in range(M.shape[1])])
    fac = np.nanmean(med)/med
    with np.errstate(divide='ignore'):
        L = np.log2(M*fac)
    return L, fac

def keep_mask_all(M):
    return np.all(M>0, axis=1) & np.all(np.isfinite(M), axis=1)

def keep_mask(M, cols_needed):
    sub = M[:, [IDX[c] for c in cols_needed]]
    return np.all(sub>0, axis=1) & np.all(np.isfinite(sub), axis=1)

def cellmeans_design(block_reps=True):
    X = np.zeros((18, len(GL)))
    for j,g in enumerate(GROUPS): X[j, GL.index(g)] = 1.0
    cols = list(GL)
    if block_reps:
        for rr in (1,2):
            X = np.column_stack([X, np.array([1.0 if REP[j]==rr else 0.0 for j in range(18)])])
            cols.append(f"day{rr+1}")
    return X, cols

def contrast_vec(cols, plus, minus):
    c = np.zeros(len(cols)); c[cols.index(plus)] = 1.0; c[cols.index(minus)] = -1.0; return c

def STANDARD_CONTRASTS(cols):
    mk = lambda a,b: contrast_vec(cols,a,b)
    return {
      'enrich_BSA':   mk('KO_BSA','C_BSA'),
      'enrich_ICAM1': mk('KO_ICAM1','C_ICAM1'),
      'enrich_PMA':   mk('KO_PMA','C_PMA'),
      'CD18_ICAM1_vs_BSA': mk('KO_ICAM1','KO_BSA'),
      'CD18_PMA_vs_BSA':   mk('KO_PMA','KO_BSA'),
      'CD18_PMA_vs_ICAM1': mk('KO_PMA','KO_ICAM1'),
      'CAAX_ICAM1_vs_BSA': mk('C_ICAM1','C_BSA'),
      'CAAX_PMA_vs_BSA':   mk('C_PMA','C_BSA'),
      'CAAX_PMA_vs_ICAM1': mk('C_PMA','C_ICAM1'),
    }

def array_weights(L, X, keep, n_iter=6, damp=0.5, bounds=(0.1,5.0)):
    Y = L[keep]; n = 18; p = X.shape[1]; w = np.ones(n)
    for _ in range(n_iter):
        W = np.diag(w); inv = np.linalg.pinv(X.T@W@X)
        h = np.clip(np.diag(X@inv@X.T@W), 0, 0.999)
        coef = inv@X.T@W@Y.T
        resid = Y.T - X@coef
        s2 = np.sum(w[:,None]*resid**2, axis=0)/(n-p); s2 = np.where(s2>0, s2, np.nan)
        z = np.nanmedian((w[:,None]*resid**2)/s2[None,:], axis=1)/(1-h)
        w = w*np.power(np.where(z>0, 1.0/z, 1.0), damp)
        w = np.clip(w, *bounds); w = w/np.exp(np.mean(np.log(w)))
    return w

def _bh(p):
    p = np.asarray(p, float); n = len(p); ok = ~np.isnan(p); pv = p[ok]; nv = len(pv)
    if nv == 0: return np.full(n, np.nan)
    o = np.argsort(pv); q = pv[o]*nv/np.arange(1, nv+1)
    for i in range(nv-2, -1, -1): q[i] = min(q[i], q[i+1])
    q = np.clip(q, 0, 1); ov = np.empty(nv); ov[o] = q
    out = np.full(n, np.nan); out[ok] = ov; return out

def _trigamma_inv(val):
    x = 0.5 + 1.0/val
    for _ in range(60):
        t = polygamma(1, x); step = t*(1-t/val)/polygamma(2, x)
        x = x + step
        if abs(step) < 1e-10: break
    return x

def fit(L, X, cols, contrasts, keep, weights=None, intensity=None, trend=True):
    n = 18; p = X.shape[1]; G = L.shape[0]
    w = np.ones(n) if weights is None else weights
    W = np.diag(w); inv = np.linalg.pinv(X.T@W@X)
    Y = L[keep]
    coef = inv@X.T@W@Y.T
    resid = Y.T - X@coef
    dfres = n - p
    s2 = np.sum(w[:,None]*resid**2, axis=0)/dfres
    s2 = np.where(s2 > 1e-12, s2, np.nan)
    Abar = Y.mean(1) if intensity is None else intensity[keep]
    ok = np.isfinite(s2) & np.isfinite(Abar)
    if trend:
        lo = lowess(np.log(s2[ok]), Abar[ok], frac=0.4, return_sorted=True)
        s2_trend = np.exp(np.interp(Abar, lo[:,0], lo[:,1]))
    else:
        s2_trend = np.full_like(s2, np.nanmedian(s2))
    z = np.log(s2[ok]) - np.log(s2_trend[ok])
    e = z - digamma(dfres/2) + np.log(dfres/2)
    evar = np.var(e, ddof=1) - polygamma(1, dfres/2)
    d0 = 2*_trigamma_inv(evar) if evar > 0 else np.inf
    if np.isfinite(d0):
        s2_post = (d0*s2_trend + dfres*s2)/(d0+dfres); df_total = dfres + d0
    else:
        s2_post = s2_trend.copy(); df_total = np.full_like(s2, 1e6)
    def scat(a):
        o = np.full(G, np.nan); o[keep] = a; return o
    out = {}
    for name, c in contrasts.items():
        beta = c@coef; u = float(c@inv@c)
        se = np.sqrt(s2_post*u); t = beta/se
        pv = 2*stats.t.sf(np.abs(t), df=df_total)
        out[name] = dict(logFC=scat(beta), t=scat(t), p=scat(pv), q=scat(_bh(pv)), se=scat(se))
    out['_meta'] = dict(d0=float(d0), dfres=dfres, s2_0_med=float(np.nanmedian(s2_trend)), weights=w)
    return out