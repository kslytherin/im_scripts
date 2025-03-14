import sys

sys.path.insert(0, "/Users/kksmith/src/stochastic-area-metric/")
from EDR_Method_Kale_Akkar import EDR_Method_Kale_Akkar
from multivariateNorm import multiNorm_ll
from test_lib import *
import pandas as pd
import numpy as np
import time
from datetime import datetime
from libcomcat.search import search
import areametric as am

time_events = search(
    starttime=datetime(1994, 1, 17, 12, 30),
    endtime=datetime(1994, 4, 18, 12, 35),
    updatedafter=datetime(1994, 2, 18, 12, 35),
)
start = time.time()

fsep = "\t"
N = 70
get_subset = True
gmm_list = ["SandikkayaAkkar2017Rjb", "CampbellBozorgnia2019"]
# gmm_list = ["CampbellBozorgnia2019"]
ver_no = "v2"
for region in ["Globe", "Hawaii", "Alaska", "EUS", "WUS"]:
    test_vals = {}

    mak_set = []
    mde_set = []
    ksq_set = []
    edr_set = []
    am_set = []
    llh_set = []
    chisq_set = []

    for g_i, gmm in enumerate(gmm_list):

        fname = (
            "/Users/kksmith/im_scripts/data/GMMout/Full_GMMvsData_"
            + gmm
            + "_"
            + region
            + "_"
            + ver_no
            + ".txt"
        )

        f2name = (
            "/Users/kksmith/im_scripts/data/GMMout/NumRecs"
            + "_"
            + gmm
            + "_"
            + region
            + "_"
            + ver_no
            + ".txt"
        )

        # fname = './KaleAkkarCode/InputData.txt'

        mydata = pd.read_csv(fname, sep=fsep, header=None)
        mydata2 = pd.read_csv(f2name, sep=fsep, header=None)

        num_sub = np.cumsum(mydata2.iloc[:, 1])[N]

        if get_subset:
            arr = np.arange(np.shape(mydata)[0])
            np.random.shuffle(arr)
            subset_i = arr[0:num_sub]
            mydata = mydata.iloc[subset_i, :]

        evtNum = mydata2.iloc[0 : N + 1, 1]

        obs = mydata.iloc[:, 0].values
        pred = mydata.iloc[:, 1].values
        TotSig = mydata.iloc[:, 2].values

        print("{}".format(region))
        if np.shape(mydata)[1] == 6 and True:
            bwEvtSig = mydata.iloc[:, 3].values
            wiEvtSig = mydata.iloc[:, 4].values
            # evtNum = mydata.iloc[:,5]
            V = multiNorm_ll(
                pred, obs, bwEvtSig, wiEvtSig, evtNum, useScipy=False, sparse=True
            )
            print("Mak log-likelihood is {}".format(V))
        RSLTs = EDR_Method_Kale_Akkar(obs, pred, TotSig)

        llh = LLH(obs, pred, TotSig)
        chisq = CHIMF(obs, pred, TotSig)

        # z = np.random.normal(0, TotSig[0], 1000)
        Y_sigma = []
        for i in range(len(obs)):
            z = np.random.normal(0, TotSig[i], 1000)
            for j in range(len(z)):
                Y_sigma.append(obs[i] + z[j])
        #
        AM = am.areaMe(obs, Y_sigma)

        print(
            "Kale and Akkar MDE_norm: {}  And for Kappa_sq: {}, EDR: {}, AM: {}, LLH: {}, CHISQ-MF: {}\n".format(
                RSLTs[0], RSLTs[1], RSLTs[2], AM, llh, chisq
            )
        )
        mak_set.append(V)
        mde_set.append(RSLTs[0])
        ksq_set.append(RSLTs[1])
        edr_set.append(RSLTs[2])
        am_set.append(AM)
        llh_set.append(llh)
        chisq_set.append(chisq)

        # test_vals_df.insert(0, "GMM", gmm_list)

    test_vals["GMM"] = gmm_list
    test_vals["Mak"] = mak_set
    test_vals["MDE_norm"] = mde_set
    test_vals["kappa_sq"] = ksq_set
    test_vals["EDR"] = edr_set
    test_vals["AM"] = am_set
    test_vals["LLH"] = llh_set
    test_vals["CHISQ-MF"] = chisq_set

    test_vals_df = pd.DataFrame(test_vals)
    output_csv = (
        f"/Users/kksmith/im_scripts/data/GMMout/Test_Results_{region}_{ver_no}.csv"
    )
    test_vals_df.to_csv(output_csv, index=False)
end = time.time()

print("Time of execution of above program is :", end - start, "s")
breakpoint()
print("end")
