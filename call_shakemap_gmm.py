import numpy as np
import pandas as pd
from shakemap_gmm.shake_gmm import ShakeGMM
from openquake.hazardlib import imt
from tqdm import tqdm

test_path = "~/Documents/GMM_data/for_kyle_conus/"
fname1 = "CONUS_HW_AK_h1_records.csv"
fullfname1 = test_path + fname1
df = pd.read_csv(fullfname1)

uids = np.unique(df["EarthquakeId"])
nevents = len(uids)
nrows = df.shape[0]

cols = list(df.keys()[35:58])
ncol = len(cols)


def str_to_imt(x):
    if x == "PGA":
        return imt.PGA()
    elif x == "PGV":
        return imt.PGV()
    else:
        period = float(x.split("=")[1].split(",")[0])
        return imt.SA(period=period)


# breakpoint()

print("start")
for x in cols:
    print(x)
    # print(str_to_imt(x))

col_imts = [str_to_imt(x) for x in cols]

pred_mean = np.full((nrows, ncol), -999, dtype=np.float64)
pred_std = np.full((nrows, ncol), -999, dtype=np.float64)
pred_phi = np.full((nrows, ncol), -999, dtype=np.float64)
pred_tau = np.full((nrows, ncol), -999, dtype=np.float64)

failed_events = []

for ind in tqdm(range(nevents)):
    sel_id = uids[ind]
    row_ind = np.where(df["EarthquakeId"] == sel_id)[0]
    try:
        shakegmm = ShakeGMM.from_eqid(sel_id)
    except:
        failed_events.append([ind, uids[ind]])
        continue

    rjb = np.array(df["JoynerBooreDistance"][row_ind])
    rrup = np.array(df["RuptureDistance"][row_ind])
    vs30 = np.array(df["Preferred_Vs30"][row_ind])
    for i, IMT in enumerate(col_imts):
        try:
            mean, std, std_types, distances, component = shakegmm.compute(
                vs30, rrup, rjb, str(IMT)
            )
            pred_mean[row_ind, i] = mean
            for istd, istd_type in zip(std, std_types):
                if istd_type == "Total Uncertainty":
                    pred_std[row_ind, i] = istd
                elif istd_type == "Intra event Uncertainty":
                    pred_phi[row_ind, i] = istd
                elif istd_type == "Inter event Uncertainty":
                    pred_tau[row_ind, i] = istd
        except KeyError:
            pred_mean[row_ind, i] = np.nan
            pred_std[row_ind, i] = np.nan
            pred_phi[row_ind, i] = np.nan
            pred_tau[row_ind, i] = np.nan
    if (ind % 900) == 0:
        np.save("pred_mean.npy", pred_mean)
        np.save("pred_std.npy", pred_std)
        np.save("pred_phi.npy", pred_phi)
        np.save("pred_tau.npy", pred_tau)
        with open("last_save.txt", "w") as outfile:
            outfile.write(f"Last save index: {ind}")
        with open("failed_events.pkl", "wb") as fp:
            pickle.dump(failed_events, fp)
