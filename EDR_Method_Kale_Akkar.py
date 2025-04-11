import numpy as np
import pandas as pd
from scipy.stats import norm

'''
Euclidean-Distance Based Ranking (EDR) Method
This script is used for performing EDR Method.
EDR Method is developed by Ozkan Kale and Sinan Akkar.

Format of Input File ("InputData.txt"):

Intensity Measure File -> 1st column: ln(Observed)
                          2nd column: ln(Estimated)
                          3rd column: Corresponding sigma values

This python code is written by Kyle Smith on 2024-08-12, and adopted
from the matlab code from Kale and Akkar.
'''


# def EDR_Method_Kale_Akkar(fname, fsep):
def EDR_Method_Kale_Akkar(obs, Y_est, sigma):

    # Input values
    dd = 0.01          # infinitesmall bandwidth
    x = 3             # the multiplier of sigma
   
    # Intensity Measure File
    #data = pd.read_csv(fname, sep=fsep, header=None)
    #N, col = np.shape(data)       # number of rows and columns in InputFile
    #breakpoint()
    N = len(obs)
    #obs = data.iloc[:, 0]              # observed data
    dmin = dd/2                # minimum discrete distance

    ## Calculation of the Kappa value
    mua = np.mean(obs)

    #Y_est = data.iloc[:, 1]
    muY = np.mean(Y_est)

    # Calculate the coefficients of Y=B1*A+B0 (Fitted line)
    x1 = (obs - np.full_like(obs, mua)) * (Y_est - np.full_like(obs, muY))
    y1 = (obs - np.full_like(obs, mua))**2
    B1 = np.sum(x1)/np.sum(y1)
    B0 = muY - B1*mua

    Yfit = np.full_like(obs, B0) + np.full_like(obs, B1)*obs
    Yc = Y_est - (Yfit - obs)    	# corrected estimations (Equation 10)

    DEorg = np.linalg.norm(obs-Y_est)          # Equation 9.b
    DEcor = np.linalg.norm(obs-Yc)         # Equation 9.c

    Kappa = DEorg/DEcor        # Equation 9.a

    ## Calculation of the EDR value
    mu_Y = Y_est
    #mu_Y = data.iloc[:, 1]
    #sigma = data.iloc[:, 2]

    mu_D = obs - mu_Y            # Equation 3.a
    s_D = sigma                  # Equation 3.b

    # Selection of an appropriate dmax value (Equation 8)
    d1c = np.abs(obs - (mu_Y - x*sigma))
    d1cmax = np.max(d1c)
    d2c = np.abs(obs - (mu_Y + x*sigma))
    d2cmax = np.max(d2c)
    dcmax = np.ceil(np.max((d1cmax, d2cmax)))

    dmax = dcmax - dd/2        # selected maximum discrete distance
    nd = len(np.arange(dmin, dmax+dd/2, dd))

    MDE = 0
    for j in np.arange(1, nd+1):
        di = dmin+(j-1)*dd
        d = np.full_like(obs, di)
        d1i = di-dd/2
        d1 = np.full_like(obs, d1i)
        d2i = di+dd/2
        d2 = np.full_like(obs, d2i)
        
        # # Calculations given in Equation 5
        P1 = norm.cdf((d1-mu_D)/s_D)-norm.cdf((-d1-mu_D)/s_D)
        P2 = norm.cdf((d2-mu_D)/s_D)-norm.cdf((-d2-mu_D)/s_D)
        
        # Calculations given in Equation 6
        P = P2-P1
        MDEi = P*d
        MDE = MDE + MDEi

    # Modified Euclidean Distance normalized by N (The first component of EDR index in terms of MDE)
    MDE_norm = np.sqrt(1/N*np.sum(MDE**2))
    # Square root of Kappa (The second component of EDR index in terms of Kappa)
    Kappa_sq = np.sqrt(Kappa)
    # Resultant EDR value (Equation 11)
    EDR = np.sqrt(Kappa*1/N*np.sum(MDE**2))

    RSLTs = np.array([MDE_norm, Kappa_sq, EDR])
    ## Store the results
    # # # # # # # # # # # # # # # # 
    return RSLTs

if __name__=="__main__":
    fname = './KaleAkkarCode/InputData.txt'
    fsep = '\t'
    mydata = pd.read_csv(fname, sep=fsep, header=None)

    obs = mydata.iloc[:,0]
    Y_est = mydata.iloc[:,1]
    sigma = mydata.iloc[:,2]
    RSLTs = EDR_Method_Kale_Akkar(obs, Y_est, sigma)
    print("Kale and Akkar MDE_norm: {}, Kappa_sq: {}, EDR: {}".format(RSLTs[0],RSLTs[1],RSLTs[2]))
    # breakpoint()
    # print("end")