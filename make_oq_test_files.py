"""
Replicating Figure 8 in Macedo et al 2021
Started on Mar 12, 2024
"""

import numpy as np
import matplotlib.pyplot as plt
import get_cgmm_im
import get_gmm_im
import pandas as pd
import copy
from mytools import *

from get_eq_stn_info import get_eq_stn_dist_vals, make_stn_array, get_event_from_rup
from gmprocess.utils.constants import UNITS, STATION_METRIC_UNITS
from openquake.hazardlib.imt import CAV, PGA, SA
import kaklamanos_relations

from esi_utils_time.ancient_time import HistoricTime

# Test over different magnitudes, rrup distances, vs30 values, HW effect(effects dips,rjb), conditional models?
# 5 Models used by Macedo et al 2021 for PGA
# should match up but need to check to make sure

# Papers of models, OQ code
# Campbell and Bozorgnia (2014), campbell_bozorgnia_2014.py
# Abrahamson et al (2014), abrahamson_2014.py
# Boore et al (2014), boore_2014.py
# Chiou and Youngs (2014), chiou_youngs_2014.py
# Idriss (2014), idriss_2014.py

# Make main file (macedo_2021_mean.csv) with the following:
# rup_mag,dist_rjb,rup_dip,dist_rrup,site_vs30,result_type,cav

################## Input ##################
makefiles = True
eid = "921"  # for Macedo replication

# 921 file from
# Rupture of downloads from:
# https://earthquake.usgs.gov/earthquakes/eventpage/usp0009eq0/shakemap/intensity
# which will lead here:
# https://earthquake.usgs.gov/product/shakemap/usp0009eq0/atlas/1594169383237/download/rupture.json


# main_gmm = "MacedoEtAl2021"
# main_gmm = "LiuMacedo2022M1SInter"
main_gmm = "LiuMacedo2022SInter"
# main_gmm = "CampbellBozorgnia2019"
# main_gmm = "SandikkayaAkkar2017Rjb"
# main_gmm = "BullockEtAl2021V1"

# con_model_set = None
# con_model_set = ["AbrahamsonEtAl2014"]  # for macedo replication
# con_model_set = ["AbrahamsonEtAl2014"]
con_model_set = ["ParkerEtAl2020SInter"]
##############################################


path_to_rupture_files = "/Users/kksmith/data/rup/"
rup_file = path_to_rupture_files + eid + "_rupture.json"

# def gmm_est(rup_file, main_gmm, con_model_set = None)

origin, event, dummy = get_event_from_rup(rup_file)

split_gmm_str = split_num_letters_caps(main_gmm)

# this will not catch all cases b/c sometimes the region is not the last word, but it catches most cases
tec_region = split_cap_words(main_gmm)[-1]
tec_str = ""
pre_tec_str = ""
if tec_region in ["Slab", "Inter", "Asc", "Stable", "Volc"]:
    if tec_region in ["Slab", "Inter"]:
        pre_tec_str = split_cap_words(main_gmm)[-2]
    tec_str = pre_tec_str.lower() + tec_region.lower() + "_"

compare = [i for i, k in enumerate(split_gmm_str) if "Et" in split_gmm_str[i]]
if len(compare) > 0:
    etal = split_gmm_str[compare[0]] + split_gmm_str[compare[0] + 1]
    gmm_str = split_gmm_str[0] + "_" + etal + "_" + split_gmm_str[-1]
    gmm_str = "_".join([split_gmm_str[0], etal, split_gmm_str[-1]]).lower()
else:
    gmm_str = "_".join(split_gmm_str).lower()
# exception so that it is consistent with other directory name
if split_gmm_str[0].lower() == "macedo":
    gmm_str = split_gmm_str[0].lower() + "_" + split_gmm_str[-1]

# con_model_set = []

# path_to_test_file_dir = "~/gmm_data/"
path_to_test_file_dir = (
    "~/oq-engine/openquake/hazardlib/tests/gsim/data/" + gmm_str + "/"
)

if main_gmm == "MacedoEtAl2021":
    # con_model_set = ["AbrahamsonEtAl2014"]
    path_to_test_file_dir = (
        "~/oq-engine/openquake/hazardlib/tests/gsim/data/" + gmm_str + "/"
    )

elif main_gmm == "MacedoEtAl2019SInter":
    con_model_set = ["AbrahamsonEtAl2018SInter"]

myimt = CAV()
conimt_set = [PGA(), SA(1.0)]
# conimt = PGA()
# conimt2 = SA(1.0)

acc1 = "%.9f"
acc2 = "%.9f"

# otime = HistoricTime(1999, 9, 20, 17, 47, 18)
# event = {"id": "usp0009eq0",
#          "lat": 23.772,
#          "lon": 120.982,
#          "depth": 33,
#          "mag": 7.7,
#          "time": otime,
#          "mech": "ALL",
#          "reference": "from usgs website",
#          "netid": "",
#          "network": "",
#          "productcode": "",
#          "locstring": ""
#          }

# show_plt = True

# breakpoint()
# make own data set
# arraydense = {"nlon": 5, "nlat": 3}
arraydense = {"nlon": 3, "nlat": 2}
lonlat_center = {"lon": event["lon"], "lat": event["lat"]}
lonlatbox = [-3, 3, -3, 3]
stn_locs = make_stn_array(lonlat_center, lonlatbox, arraydense)
# breakpoint()
rjbset, rhypset, rrupset, repiset, o_dis, eq_width, eq_ztor, eq_dip = (
    get_eq_stn_dist_vals(origin, rup_file, stn_locs)
)

eqmagset = np.array([copy.deepcopy(event["mag"])])
dipset = np.array([copy.deepcopy(eq_dip)])
vs30set = np.arange(100, 1000, 3)
# dipset = np.arange(0, 90, 10)
# eqmagset = np.arange(3, 8, 0.5)

# read-in data set
# get data
lnmean_set = []
mean_set = []
sigma_set = []
tau_set = []
phi_set = []

lnmean_set_con = []
mean_set_con = []
sigma_set_con = []
tau_set_con = []
phi_set_con = []
if len(conimt_set) > 1:
    lnmean_set_con_alt = []
    mean_set_con_alt = []
    sigma_set_con_alt = []
    tau_set_con_alt = []
    phi_set_con_alt = []

# Will need to make this better in the future
rxset = o_dis["rx"]
ry0set = o_dis["ry0"]
# num_dist = 5
# rjbset = np.logspace(-1, np.log10(400), num_dist)
# rrupset = np.logspace(-1, np.log10(400), num_dist)

# the following are the available GMM models in OQ that we want for
# Macedo et al (2021)
# con_leg_labels = ["BSSA_2014_nga"]
# con_col = ["b"]

case_i = 0
# breakpoint()
for mag_i, mag_num in enumerate(eqmagset):
    for v_i, vs30_val in enumerate(vs30set):
        # Specifying earthquake parameters
        myeqparams = {
            "mag": mag_num,
            "rrup": rrupset,
            "rhypo": rhypset,
            "hypo_depth": event["depth"],
            "vs30": vs30_val,  # 424,
            "ztor": eq_ztor,  # 20,
            "backarc": False,
            "rake": event["rake"],
            "rx": rrupset,
            "width": eq_width,
            "dip": eq_dip,
            "ry0": ry0set,
            "rjb": rjbset,
            "z2pt5": 0.2,
            "z1pt0": 0.1,
            "vs30measured": True,
        }
        # breakpoint()
        if con_model_set is not None:
            # if len(con_model_set) > 0:
            for m_i, con_model in enumerate(con_model_set):
                # myeqparams, main_gmm, con_model, myimt
                lnmean_pre, sigma, tau, phi, ctx = get_cgmm_im.get_cgmm_im(
                    myeqparams, main_gmm, con_model, myimt
                )
                # breakpoint()

                lnmean_con, sigma_con, tau_con, phi_con, ctx_con = (
                    get_gmm_im.get_gmm_im(myeqparams, con_model, conimt_set[0])
                )

                if len(conimt_set) > 1:
                    (
                        lnmean_con_alt,
                        sigma_con_alt,
                        tau_con_alt,
                        phi_con_alt,
                        ctx_con_alt,
                    ) = get_gmm_im.get_gmm_im(myeqparams, con_model, conimt_set[1])

                newctx = pd.DataFrame.from_records(ctx)
                newctx_con = pd.DataFrame.from_records(ctx_con)
                if case_i == 0:
                    master_ctx = copy.deepcopy(newctx)
                    master_ctx_con = copy.deepcopy(newctx_con)
                else:
                    master_ctx = pd.concat([master_ctx, newctx], ignore_index=True)
                    master_ctx_con = pd.concat(
                        [master_ctx_con, newctx_con], ignore_index=True
                    )

                lnmean = lnmean_pre
                mean = np.exp(lnmean)
                mean_con = np.exp(lnmean_con)
                mean_con_alt = np.exp(lnmean_con_alt)

                lnmean_set = np.concatenate((lnmean_set, lnmean))
                mean_set = np.concatenate((mean_set, mean))
                sigma_set = np.concatenate((sigma_set, sigma))
                tau_set = np.concatenate((tau_set, tau))
                phi_set = np.concatenate((phi_set, phi))

                lnmean_set_con = np.concatenate((lnmean_set_con, lnmean_con))
                mean_set_con = np.concatenate((mean_set_con, mean_con))
                sigma_set_con = np.concatenate((sigma_set_con, sigma_con))
                tau_set_con = np.concatenate((tau_set_con, tau_con))
                phi_set_con = np.concatenate((phi_set_con, phi_con))

                lnmean_set_con_alt = np.concatenate(
                    (lnmean_set_con_alt, lnmean_con_alt)
                )
                mean_set_con_alt = np.concatenate((mean_set_con_alt, mean_con_alt))
                sigma_set_con_alt = np.concatenate((sigma_set_con_alt, sigma_con_alt))
                tau_set_con_alt = np.concatenate((tau_set_con_alt, tau_con_alt))
                phi_set_con_alt = np.concatenate((phi_set_con_alt, phi_con_alt))
        else:
            lnmean_pre, sigma, tau, phi, ctx = get_gmm_im.get_gmm_im(
                myeqparams, main_gmm, myimt
            )

            newctx = pd.DataFrame.from_records(ctx)
            if case_i == 0:
                master_ctx = copy.deepcopy(newctx)
            else:
                master_ctx = pd.concat([master_ctx, newctx], ignore_index=True)

            lnmean = lnmean_pre
            mean = np.exp(lnmean)
            # breakpoint()
            lnmean_set = np.concatenate((lnmean_set, lnmean))
            mean_set = np.concatenate((mean_set, mean))
            sigma_set = np.concatenate((sigma_set, sigma))
            tau_set = np.concatenate((tau_set, tau))
            phi_set = np.concatenate((phi_set, phi))
        case_i += 1


# change names to appropriate naming convention
master_ctx = master_ctx.rename(
    columns={
        "mag": "rup_mag",
        "hypo_depth": "rup_hypo_depth",
        "rake": "rup_rake",
        "dip": "rup_dip",
        "width": "rup_width",
        "ztor": "rup_ztor",
        "rrup": "dist_rrup",
        "repi": "dist_repi",
        "rjb": "dist_rjb",
        "rhypo": "dist_rhypo",
        "rx": "dist_rx",
        "ry0": "dist_ry0",
        "z2pt5": "site_z2pt5",
        "z1pt0": "site_z1pt0",
        "vs30": "site_vs30",
        "backarc": "site_backarc",
        "vs30measured": "site_vs30measured",
    }
)

ResultTypes = ["MEAN", "TOTAL_STDDEV", "INTER_EVENT_STDDEV", "INTRA_EVENT_STDDEV"]
ftags = [x.lower() for x in ResultTypes]
# breakpoint()

for st_i, stat in enumerate([mean_set, sigma_set, tau_set, phi_set]):
    master_ctx["result_type"] = ResultTypes[st_i]
    # breakpoint()
    master_ctx["cav"] = stat
    if makefiles:
        master_ctx.to_csv(
            path_to_test_file_dir + gmm_str + "_" + tec_str + ftags[st_i] + ".csv",
            index=False,
            float_format=acc1,
        )
    master_ctx = master_ctx.drop(["cav", "result_type"], axis=1)

master_ctx = master_ctx.drop(
    [
        "rup_hypo_depth",
        "rup_rake",
        "rup_width",
        "rup_ztor",
        "dist_rhypo",
        "dist_ry0",
        "site_z1pt0",
        "site_z2pt5",
        "site_vs30measured",
        "site_backarc",
    ],
    axis=1,
)

master_ctx = master_ctx.rename(
    columns={
        "rup_mag": "mag",
        "dist_rrup": "rrup",
        "site_vs30": "vs30",
        "rup_dip": "dip",
        "dist_rjb": "rjb",
        "dist_rx": "rx",
    }
)


# New method puts this first rather than PGA stuff
master_ctx["Global_MEAN"] = mean_set
master_ctx["Global_SIG"] = sigma_set
master_ctx["Global_TAU"] = tau_set
master_ctx["Global_PHI"] = phi_set

if con_model_set is not None:
    master_ctx[str(conimt_set[0]) + "_MEAN"] = mean_set_con
    master_ctx[str(conimt_set[0]) + "_TOTAL_STDDEV"] = sigma_set_con
    master_ctx[str(conimt_set[0]) + "_INTER_EVENT_STDDEV"] = tau_set_con
    master_ctx[str(conimt_set[0]) + "_INTRA_EVENT_STDDEV"] = phi_set_con

    if len(conimt_set) > 1:
        master_ctx[str(conimt_set[1]) + "_MEAN"] = mean_set_con_alt
        master_ctx[str(conimt_set[1]) + "_TOTAL_STDDEV"] = sigma_set_con_alt
        master_ctx[str(conimt_set[1]) + "_INTER_EVENT_STDDEV"] = tau_set_con_alt
        master_ctx[str(conimt_set[1]) + "_INTRA_EVENT_STDDEV"] = phi_set_con_alt

if makefiles:

    master_ctx.to_csv(
        path_to_test_file_dir + gmm_str + "_" + tec_str + "conditioning_gmvs.csv",
        index=False,
        float_format=acc2,
    )

# print("sigma={}".format(sigma))

# if mag_i == 1:
#    breakpoint()

# # plotting means
# fig = plt.figure()
# ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# dist_meas = "rrup"
# for m_i, con_model in enumerate(con_model_set):
#     ax.loglog(
#         ctx[dist_meas],
#         mean_set[con_model],
#         color=con_col[m_i],
#         label=con_leg_labels[m_i],
#     )
# ax.set_xlabel("$R_{rup}$, " + STATION_METRIC_UNITS[dist_meas + "_mean"])
# #ax.set_ylabel(str(myimt) + ", " + UNITS[str(myimt).lower()])
# if str(myimt) == "CAV":
#     ax.set_ylabel(str(myimt) + ", m/s")
# ax.yaxis.set_ticks_position("both")
# ax.xaxis.set_ticks_position("both")

# ax.set_title(
#     main_gmm + ", M=" + str(myeqparams["mag"]) + ", Vs30="
#     + str(myeqparams["vs30"])
# )
# ax.set(xlim=(1, 100), ylim=(0.01, 100))
# ax.set_box_aspect(1)
# ax.legend()

# # plotting standard devations
# fig = plt.figure()
# ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# dist_meas = "rrup"
# for m_i, con_model in enumerate(con_model_set):
#     ax.semilogx(
#         ctx[dist_meas],
#         sigma_set[con_model],
#         color=con_col[m_i],
#         label=con_leg_labels[m_i],
#     )
# ax.set_xlabel("$R_{rup}$, " + STATION_METRIC_UNITS[dist_meas + "_mean"])
# # ax.set_ylabel(str(myimt) + " $ \sigma$, " + UNITS[str(myimt).lower()])
# ax.set_ylabel(str(myimt) + " $ \sigma$")
# ax.yaxis.set_ticks_position("both")
# ax.xaxis.set_ticks_position("both")

# ax.set_title(
#     main_gmm + ", M=" + str(myeqparams["mag"]) + ", Vs30="
#     + str(myeqparams["vs30"])
# )
# ax.set(xlim=(1, 1000), ylim=(0.35, 0.85))
# ax.set_box_aspect(1)
# ax.legend()
# if show_plt:
#    plt.show()

breakpoint()
print("end")
