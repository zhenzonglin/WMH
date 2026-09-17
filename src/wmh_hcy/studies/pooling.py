"""Rubin scalar pooling and Li-Raghunathan-Rubin D1 multivariate pooling."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from ..common import DataError
from ..imputation import pool_scalar


def joint_test(estimates, covariances):
    q, u = np.asarray(estimates, float), np.asarray(covariances, float)
    m, k = q.shape
    if not np.isfinite(q).all() or not np.isfinite(u).all():
        raise DataError("Nonfinite input to MI joint test")
    if k == 1:
        result = pool_scalar(q[:, 0], u[:, 0, 0])
        return {"p": result["p"], "statistic": (result["estimate"]/result["se"])**2,
                    "df1": 1, "df2": result["df"], "method": "Rubin scalar Wald", "m": m}
    mean, within = q.mean(axis=0), u.mean(axis=0)
    if np.linalg.eigvalsh(within).min() <= 0:
        raise DataError("Singular within-imputation covariance in D1")
    between = np.cov(q, rowvar=False, ddof=1) if m > 1 else np.zeros((k, k))
    r = max(0.0, float((1 + 1/m)*np.trace(np.linalg.solve(within, between))/k))
    statistic = float(mean @ np.linalg.solve(within, mean) / (k*(1+r)))
    t = k*(m-1)
    if m == 1 or r < 1e-14:
        df2, p = np.inf, stats.chi2.sf(k*statistic, k)
    else:
        df2 = 4+(t-4)*(1+(1-2/t)/r)**2 if t > 4 else t*(1+1/k)*(1+1/r)**2/2
        p = stats.f.sf(statistic, k, df2)
    return {"p": float(p), "statistic": statistic, "df1": k, "df2": float(df2), "r": r,
                "method": "D1 pooled multivariate Wald", "m": m}


def terms_test(fits, terms):
    indices = [fits[0].terms.index(t) for t in terms]
    return joint_test([f.params[indices] for f in fits],
                      [f.covariance[np.ix_(indices, indices)] for f in fits])


def coefficients(fits, family):
    rows = []
    for j, term in enumerate(fits[0].terms):
        r = pool_scalar([f.params[j] for f in fits], [f.covariance[j, j] for f in fits])
        r.update(term=term, scale="difference" if family == "ols" else
                 "HR" if family == "cox" else "odds_ratio" if family == "binary" else "relative_probability_ratio")
        if family != "ols":
            r.update(ratio=float(np.exp(r["estimate"])), ratio_lower=float(np.exp(r["lower"])),
                     ratio_upper=float(np.exp(r["upper"])))
        rows.append(r)
    return pd.DataFrame(rows)


def contrast(fits, vectors):
    if isinstance(vectors, np.ndarray) and vectors.ndim == 1:
        vectors = [vectors]*len(fits)
    return pool_scalar([v @ f.params for v, f in zip(vectors, fits, strict=True)],
                       [v @ f.covariance @ v for v, f in zip(vectors, fits, strict=True)])
