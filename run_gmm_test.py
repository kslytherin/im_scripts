import sys
import os.path

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
import copy

time_events = search(
    starttime=datetime(1994, 1, 17, 12, 30),
    endtime=datetime(1994, 4, 18, 12, 35),
    updatedafter=datetime(1994, 2, 18, 12, 35),
)
start = time.time()

fsep = "\t"
N = 70  # seems like an arbitraty number, not sure why I picked it, must re-examine
get_subset = True
# gmm_list = ["SandikkayaAkkar2017Rjb", "CampbellBozorgnia2019"]
#    "LiuMacedo2022SInter-ParkerEtAl2020SInter",
gmm_list = [
    "CampbellBozorgnia2019",
    "SandikkayaAkkar2017Rjb",
    "MacedoEtAl2021-BooreEtAl2014",
    "LiuMacedo2022SSlab-ParkerEtAl2020SSlab",
    "LiuMacedo2022M1SSlab",
    "LiuMacedo2022M2SSlab",
    "LiuMacedo2022SInter-ParkerEtAl2020SInter",
    "LiuMacedo2022M1SInter",
    "LiuMacedo2022M2SInter",
]
# gmm_list = ["CampbellBozorgnia2019"]
ver_no = "v2"

# selected_tect_reg = [
#         "Active",
#         "Stable",
#         "SubductionInterface",
#         "SubductionIntraslab",
#         "Volcanic",
#     ]

selected_tect_reg = [
    "Active",
    "SubductionInterface",
    "SubductionIntraslab",
]

# for geo_reg in ["Globe", "Hawaii", "Alaska", "EUS", "WUS"]:
for geo_reg in ["Alaska", "WUS"]:
    print("{}".format(geo_reg))
    for tect_reg in selected_tect_reg:
        print("{}".format(tect_reg))
        test_vals = {}

        mak_set = []
        mde_set = []
        ksq_set = []
        edr_set = []
        am_set = []
        llh_set = []
        chisq_set = []
        updated_gmm_list = copy.deepcopy(gmm_list)
        if tect_reg == "SubductionIntraslab":
            updated_gmm_list.remove("LiuMacedo2022M1SInter")
            updated_gmm_list.remove("LiuMacedo2022M2SInter")
            updated_gmm_list.remove("LiuMacedo2022SInter-ParkerEtAl2020SInter")
        elif tect_reg in ["SubductionInterface", "Volcanic", "Stable"]:
            updated_gmm_list.remove("LiuMacedo2022M1SSlab")
            updated_gmm_list.remove("LiuMacedo2022M2SSlab")
            updated_gmm_list.remove("LiuMacedo2022SSlab-ParkerEtAl2020SSlab")
        elif tect_reg == "Active":
            updated_gmm_list.remove("LiuMacedo2022M1SInter")
            updated_gmm_list.remove("LiuMacedo2022M2SInter")
            updated_gmm_list.remove("LiuMacedo2022SInter-ParkerEtAl2020SInter")
            updated_gmm_list.remove("LiuMacedo2022M1SSlab")
            updated_gmm_list.remove("LiuMacedo2022M2SSlab")
            updated_gmm_list.remove("LiuMacedo2022SSlab-ParkerEtAl2020SSlab")

        for g_i, gmm in enumerate(updated_gmm_list):

            fname = (
                "/Users/kksmith/data/GMMout/Full_GMMvsData_"
                + gmm
                + "_"
                + geo_reg
                + "_"
                + tect_reg
                + "_"
                + ver_no
                + ".txt"
            )

            f2name = (
                "/Users/kksmith/data/GMMout/NumRecs"
                + "_"
                + gmm
                + "_"
                + geo_reg
                + "_"
                + tect_reg
                + "_"
                + ver_no
                + ".txt"
            )

            # fname = './KaleAkkarCode/InputData.txt'
            if not os.path.exists(fname):
                print("File " + fname + " does not exist.")
                breakpoint()
                continue
            mydata = pd.read_csv(fname, sep=fsep, header=None)
            mydata2 = pd.read_csv(f2name, sep=fsep, header=None)
            N = np.shape(mydata2)[1]
            try:
                num_sub = np.cumsum(mydata2.iloc[:, 1])[N]
            except:
                breakpoint()

            if get_subset:
                arr = np.arange(np.shape(mydata)[0])
                np.random.shuffle(arr)
                subset_i = arr[0:num_sub]
                mydata = mydata.iloc[subset_i, :]

            evtNum = mydata2.iloc[0 : N + 1, 1]

            obs = mydata.iloc[:, 0].values
            pred = mydata.iloc[:, 1].values
            TotSig = mydata.iloc[:, 2].values

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

            # test_vals_df.insert(0, "GMM", updated_gmm_list)

        test_vals["GMM"] = updated_gmm_list
        test_vals["Mak"] = mak_set
        test_vals["MDE_norm"] = mde_set
        test_vals["kappa_sq"] = ksq_set
        test_vals["EDR"] = edr_set
        test_vals["AM"] = am_set
        test_vals["LLH"] = llh_set
        test_vals["CHISQ-MF"] = chisq_set

        try:
            test_vals_df = pd.DataFrame(test_vals)
        except:
            print("No files for " + geo_reg + "-" + tect_reg)
            breakpoint()
            continue
        output_csv = (
            "/Users/kksmith/data/GMMout/Test_Results_"
            + geo_reg
            + "_"
            + tect_reg
            + "_"
            + ver_no
            + ".csv"
        )
        test_vals_df.to_csv(output_csv, index=False, float_format="%.4f")

        # Create a heatmap of the test_vals_df
end = time.time()

print("Time of execution of above program is :", end - start, "s")
breakpoint()
print("end")
