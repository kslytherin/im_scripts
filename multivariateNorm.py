"""
Likelihood of multivariate normal distribution
Written by Mak et al. (2017)
Mak, S., R. A. Clements, and D. Schorlemmer (2017). Empirical evaluation of hierachical ground motion models: Score uncertainty and model weighting. Bulletin of the Seismological Society of America. Vol 107 (2). DOI: 10.1785/0120160232.


Require packages:
 1. numpy
 2. scipy
 3. scikit-sparse
"""
import numpy as np
from scipy.linalg import solve
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix, diags
from scipy.stats import norm, multivariate_normal
from sksparse.cholmod import cholesky
#from scikit.sparse.cholmod import cholesky

def multiNorm_ll(pred, obs, bwEvtSig, wiEvtSig, evtNum, useScipy=False, sparse=True):
  # Equation 7 of the main article.
  #   (not multiplied by -1)
  #
  # Return the log-likelihood of multivariate normal distribution
  #
  # INPUT ARGUMENTS
  #
  #  pred
  #    Predicted values {ndarray of N members}
  #    in natural log scale
  #
  #  obs
  #    Observed values {ndarray of N members}
  #    in natural log scale
  #
  #  bwEvtSig
  #    between-event sigma
  #    float: Same value for all events and records
  #    ndarray of M members: Same for all records for the same event,
  #                          but different among events
  #    ndarray of N members: Different for each record
  #    in natural log scale
  #
  #  wiEvtSig
  #    float: Same value for all events and records
  #    ndarray of N members: Different for each record
  #    in natural log scale
  #
  # (ATTENTION: sigma is standard deviation, not variance)
  # (M is the number of earthquakes)
  # (N is the total number of records)
  #
  #  evtNum
  #    list
  #    Length is total number of events
  #    Each member is the number of records in that event
  #    e.g.: [10,13,5] is 3 events, with 10,13,5 records, respectively
  #
  #  useScipy:
  #    True: Use scipy multivariate_normal()
  #    False: Use self-defined functions (allows sparse matrix)
  #
  #  sparse:
  #     True: Optimized using sparse matrix
  #     Require useScipy=False

  #breakpoint()
  # Form variance-covariance matrix
  V = form_cov(bwEvtSig, wiEvtSig, evtNum, sparse)
  N = sum(evtNum)
  #breakpoint()
  # Call scipy.stats function
  if useScipy:
    if sparse:  
      return multivariate_normal(pred, V.todense()).logpdf(obs)
    else:
      return multivariate_normal(pred, V).logpdf(obs)

  # Self-defined formula
  b = obs - pred
  if sparse:
    factor = cholesky(V)
    # will convert csr to csv (please ignore warning)
    sign, logdetV = factor.slogdet()
    llh = -N*np.log(2*np.pi)/2 - logdetV/2 - b.T.dot( spsolve(V,b) )/2
  else:
    sign, logdetV = np.linalg.slogdet(V)
    llh = -N*np.log(2*np.pi)/2 - logdetV/2 - b.T.dot( solve(V,b) )/2

  return llh
######################## END

def form_cov(bwEvtSig, wiEvtSig, evtNum, sparse=True):
  # Return covariance matrix
  #
  # See multiNorm_ll() for input
  #
  # Optimized using sparse matrix

  # number of events
  M = len(evtNum)
  # number of records altogether
  N = sum(evtNum)

  if not hasattr(wiEvtSig, "__len__"):
    wiEvtSig = np.array([wiEvtSig]*N)  
  if not hasattr(bwEvtSig, "__len__"):
    bwEvtSig = np.array([bwEvtSig]*M)

  # breakpoint()
  # V = R + Z*G*Z'
  if len(bwEvtSig)==M:

    # Use sparse matrix
    if sparse:
      Z = lil_matrix((N,M))
      for i in range(M):
        beg = sum(evtNum[:i])
        end = sum(evtNum[:i+1])
        Z[beg:end, i] = 1
      Z = Z.tocsc()
      G = diags(bwEvtSig**2, offsets=0) # MxM matrix
      R = diags(wiEvtSig**2, offsets=0) # NxN matrix
      V = Z.dot(G).dot(Z.transpose()) + R
      breakpoint()
  
    # Use dense matrix
    else:
      Z = np.zeros((N,M))
      for i in range(M):
        beg = sum(evtNum[:i])
        end = sum(evtNum[:i+1])
        Z[beg:end, i] = 1
      G = np.diag(bwEvtSig**2) # MxM matrix
      R = np.diag(wiEvtSig**2) # NxN matrix
      V = Z.dot(G).dot(Z.T) + R

  # V = R + Zg * Zg'
  elif len(bwEvtSig)==N:

    # Use sparse matrix
    # Note that lil_matrix cannot slice axis 0
    # So construct Zg' instead of Zg
    if sparse:
      Zg_t = lil_matrix((M, N))
      for i in range(M):
        beg = sum(evtNum[:i])
        end = sum(evtNum[:i+1])
        Zg_t[i, beg:end] = bwEvtSig[beg:end]
      Zg_t = Zg_t.tocsc()
      R = diags(wiEvtSig**2, offsets=0) # NxN matrix
      V = Zg_t.transpose().dot(Zg_t) + R
      # V is now csr instead of csc, due to the transpose
  
    # Use dense matrix
    else:
      Zg = np.zeros((N,M))
      for i in range(M):
        beg = sum(evtNum[:i])
        end = sum(evtNum[:i+1])
        Zg[beg:end, i] = bwEvtSig[beg:end]
      R = np.diag(wiEvtSig**2) # NxN matrix
      V = Zg.dot(Zg.T) + R

  else:
    raise ValueError("Dimension of bwEvtSig not right: M=%s" %(bwEvtSig.shape))

  #breakpoint()
  return V
######################## END



if __name__=="__main__":
  # Sample run
  #
  # Between-event sigma and within-evet sigma given in Example A2
  # Set observed and predicted ground motions to zero
  #
  # Ignore warning

  bwEvtSig = np.array([0.30, 0.31, 0.32, 0.32, 0.33])
  wiEvtSig = np.array([0.45, 0.46, 0.47, 0.47, 0.48])
  evtNum = [2, 3, 2, 3, 3]
  #evtNum = [2, 3]
  pred = np.zeros(5) # natural log of predicted ground motion
  obs = np.zeros(5)  # natural log of observed gorund motion

  V = multiNorm_ll(pred, obs, bwEvtSig, wiEvtSig, evtNum)
  print("Mak log-likelihood is {}".format(V))
  #print -multiNorm_ll(pred, obs, bwEvtSig, wiEvtSig, evtNum)
