import math


# log likelihood test (LLH) of Scherbaum et al. (2009).
def LLH(data, pred, s):
    from scipy.stats import norm

    f = 0
    N = len(data)
    for i in range(N):
        v = norm.pdf(data[i], pred[i], s[i])  # s is the sigma of gmpe
        l = math.log2(v)
        f = f + l
    llh = f / N
    return -llh


# chi-square Misfit (CHISQ-MF) test.
def CHIMF(data, pred_without_sigma, s):
    f = 0
    N = len(data)
    for i in range(N):
        sq = ((data[i] - pred_without_sigma[i]) / s[i]) ** 2
        f = f + sq
    CHI = f / N
    return CHI
