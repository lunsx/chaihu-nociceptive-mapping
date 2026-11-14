"""Pure-Python two-sample MR estimators (Wald, IVW FE/RE, Egger, weighted median)
with Cochran's Q heterogeneity and the Egger intercept pleiotropy test."""
import numpy as np
from scipy import stats


def wald_ratio(bx, by, sx, sy):
    theta = by / bx
    se = np.sqrt((sy ** 2) / (bx ** 2) + (by ** 2 * sx ** 2) / (bx ** 4))
    return theta, se


def _p(theta, se):
    return 2 * stats.norm.sf(abs(theta / se))


def ivw_fixed(bx, by, sy):
    w = 1 / sy ** 2
    theta = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    se = 1 / np.sqrt(np.sum(w * bx ** 2))
    return theta, se


def ivw_random(bx, by, sy):
    # multiplicative random effects
    w = 1 / sy ** 2
    theta_f = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    q = np.sum(w * (by - theta_f * bx) ** 2)
    k = len(bx)
    phi = max(1.0, q / (k - 1))
    wr = w / phi
    theta = np.sum(wr * bx * by) / np.sum(wr * bx ** 2)
    se = 1 / np.sqrt(np.sum(wr * bx ** 2))
    return theta, se, q, phi


def egger(bx, by, sy):
    # sign alignment handled by caller; orient bx positive
    sign = np.sign(bx)
    sign[sign == 0] = 1
    bxo = bx * sign
    byo = by * sign
    w = 1 / sy ** 2
    sw = np.sum(w)
    sxw = np.sum(w * bxo)
    syw = np.sum(w * byo)
    sxx = np.sum(w * bxo ** 2)
    sxy = np.sum(w * bxo * byo)
    denom = sw * sxx - sxw ** 2
    intercept = (syw * sxx - sxw * sxy) / denom
    slope = (sw * sxy - sxw * syw) / denom
    # residual variance for SE
    fitted = intercept + slope * bxo
    rss = np.sum(w * (byo - fitted) ** 2)
    dfree = max(1, len(bx) - 2)
    sigma2 = rss / dfree
    se_int = np.sqrt(sigma2 * sxx / denom)
    se_slope = np.sqrt(sigma2 * sw / denom)
    return slope, se_slope, intercept, se_int


def weighted_median(bx, by, sx, sy, n_boot=10000, seed=42):
    ratio = by / bx
    se_ratio = np.sqrt((sy ** 2) / (bx ** 2) + (by ** 2 * sx ** 2) / (bx ** 4))
    weights = 1 / se_ratio ** 2
    order = np.argsort(ratio)
    rs = ratio[order]
    ws = weights[order]
    cdf = np.cumsum(ws) - 0.5 * ws
    cdf /= np.sum(ws)
    theta = np.interp(0.5, cdf, rs)

    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        bxs = np.random.normal(bx, sx)
        bys = np.random.normal(by, sy)
        rr = bys / bxs
        o = np.argsort(rr)
        rr, ww = rr[o], weights[o]
        c = np.cumsum(ww) - 0.5 * ww
        c /= np.sum(ww)
        boots[i] = np.interp(0.5, c, rr)
    se = np.std(boots, ddof=1)
    return theta, se


def cochran_q(bx, by, sy):
    w = 1 / sy ** 2
    theta = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    q = np.sum(w * (by - theta * bx) ** 2)
    p = stats.chi2.sf(q, len(bx) - 1)
    return q, p


def run_all(bx, by, sx, sy):
    bx = np.asarray(bx, float)
    by = np.asarray(by, float)
    sx = np.asarray(sx, float)
    sy = np.asarray(sy, float)
    out = {}
    if len(bx) == 1:
        t, se = wald_ratio(bx[0], by[0], sx[0], sy[0])
        out["Wald ratio"] = (t, se, _p(t, se))
        return out
    t, se = ivw_fixed(bx, by, sy)
    out["IVW fixed"] = (t, se, _p(t, se))
    tr, ser, q, phi = ivw_random(bx, by, sy)
    out["IVW random"] = (tr, ser, _p(tr, ser))
    slope, se_sl, intc, se_int = egger(bx, by, sy)
    out["MR-Egger"] = (slope, se_sl, _p(slope, se_sl))
    out["Egger intercept"] = (intc, se_int, _p(intc, se_int))
    tw, sew = weighted_median(bx, by, sx, sy)
    out["Weighted median"] = (tw, sew, _p(tw, sew))
    qval, qp = cochran_q(bx, by, sy)
    out["Cochran Q"] = (qval, None, qp)
    return out
