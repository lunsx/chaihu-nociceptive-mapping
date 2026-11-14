"""Approximate Bayes factor colocalization (coloc.abf style), single causal variant
assumption per trait. Log-ABF follows Wakefield (2009). Sums are kept in log space
with logsumexp so absolute evidence is preserved."""
import numpy as np
from scipy.special import logsumexp


def log_abf(beta, varbeta, w):
    beta = np.asarray(beta, float)
    v = np.asarray(varbeta, float)
    z2 = beta ** 2 / v
    return 0.5 * (np.log(v / (v + w)) + z2 * w / (v + w))


def coloc_abf(beta1, var1, beta2, var2, w1=0.04 ** 2, w2=0.04 ** 2,
              p1=1e-4, p2=1e-4, p12=1e-5):
    l1 = log_abf(beta1, var1, w1)
    l2 = log_abf(beta2, var2, w2)
    log_s1 = logsumexp(l1)
    log_s2 = logsumexp(l2)
    log_s12 = logsumexp(l1 + l2)

    h0 = 0.0
    h1 = np.log(p1) + log_s1
    h2 = np.log(p2) + log_s2
    # H3: distinct variants -> s1*s2 - s12
    gap = log_s12 - log_s1 - log_s2
    if gap < -1e-12:
        h3 = np.log(p1) + np.log(p2) + log_s1 + log_s2 + np.log1p(-np.exp(gap))
    else:
        h3 = -np.inf
    h4 = np.log(p12) + log_s12

    allh = np.array([h0, h1, h2, h3, h4])
    finite = allh[np.isfinite(allh)]
    pp = np.exp(allh - np.max(finite))
    pp[~np.isfinite(allh)] = 0.0
    pp /= pp.sum()
    keys = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    return dict(zip(keys, pp))
