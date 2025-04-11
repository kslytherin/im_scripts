import numpy as np
import os
import sys
import mytools
import matplotlib.pyplot as plt
import pandas as pd
import get_gmm_im
import get_cgmm_im
import kaklamanos_relations
from scipy.stats import norm
from gmprocess.utils.constants import UNITS, STATION_METRIC_UNITS
from openquake.hazardlib.imt import PGA, CAV, IA, PGV
import statsmodels.formula.api as smf

from convert_GMM_units import *
from get_eq_stn_info import get_eq_stn_dist_vals, get_event_from_rup
import time
import copy

# fix kaklamanos operations. Use gmprocess? will only really work if there is a rupture file for every earthquake

# record start time
start = time.time()
# will need to do imt labels later

############## Inputs ##############
makefiles = True
myfig_ops = {"makefigs": True, "showfigs": True, "printfigs": True}
myimt = CAV()
dta_ver = "v2"

# gmm_model = "CampbellBozorgnia2019"
# gmm_model = "SandikkayaAkkar2017Rjb"
# gmm_model = ("MacedoEtAl2021","BooreEtAl2014")
# gmm_model = ("LiuMacedo2022SInter", "ParkerEtAl2020SInter")
# gmm_model = ("LiuMacedo2022SInter", "ParkerEtAl2020SInter")
# gmm_model = "LiuMacedo2022M1SInter"
# gmm_model = "LiuMacedo2022M2SInter"
# gmm_model = "LiuMacedo2022M1SSlab"
# gmm_model = "LiuMacedo2022M2SSlab"
# gmm_model = ("MacedoEtAl2021","AbrahamsonEtAl2014")
# gmm_model = ("MacedoEtAl2021","ChiouYoungs2014")
# gmm_model = ("MacedoEtAl2021","CampbellBozorgnia2014")

#####################################

path_to_rupture_files = "/Users/kksmith/data/rup/"
out_path_0 = "~/data/GMMout/"

if dta_ver == "v2":
    test_path = "~/Documents/GMM_data/for_kyle_conus/"
    # out_path = "~/output_conus/"
    out_path = "~/output_conus_test/"
    USname = "CONUS"
    fname1 = "CONUS_HW_AK_h1_records.csv"
    fname2 = "dsgmd_h2.csv"
    magstr = "PreferredMagnitude"
    vs30str = "Preferred_Vs30"
    dta_ver = "v2"
elif dta_ver == "v1":
    test_path = "~/Documents/GMM_data/for_kyle/"
    out_path = "~/output/"
    USname = "West_US"
    fname1 = "wus_default_metrics_channels(component=h1)_filtered_vs30.csv"
    fname2 = "wus_default_metrics_channels(component=h2)_filtered.csv"
    magstr = "EarthquakeMagnitude"
    vs30str = "vs30"
    dta_ver = "v1"

dpi_res = 250
if gmm_model == "SandikkayaAkkar2017Rjb":
    pre_dist_meas = "Rjb"
else:
    pre_dist_meas = "Rrup"
dist_meas = pre_dist_meas.lower()
# dist_meas_label = pre_dist_meas
dist_meas_label = "Distance to Earthquake"

# ind_geo_sec = ["WUS", "EUS", "Alaska", "Hawaii"]
ind_geo_sec = ["WUS", "Alaska"]
# if False, it uses our custom boundaries for plotting, else it does it "on the fly":w

# main_gmm_set = [("MacedoEtAl2021","CampbellBozorgnia2014"), "CampbellBozorgnia2019", "SandikkayaAkkar2017Rjb"]
# gmm_type = ["CGMM", "TGMM", "TGMM"]
# main_leg_labels = ["M21-CB14", "CB19", "S17"]
# main_col = ["r", "g", "b"]

# if False, will use our own input limits

if isinstance(gmm_model, tuple):
    gmm_model_str = gmm_model[0] + "-" + gmm_model[1]
    if gmm_model[0] == "MacedoEtAl2021":
        if gmm_model[1] == "BooreEtAl2014":
            nice_gmm_model_str = "Macedo et al. (2021) [Boore et al. (2014)]"
        elif gmm_model[1] == "AbrahamsonEtAl2014":
            nice_gmm_model_str = "Macedo et al. (2021) [Abrahamson et al. (2014)]"
        elif gmm_model[1] == "ChiouYoungs2014":
            nice_gmm_model_str = "Macedo et al. (2021) [Chiou and Youngs (2014)]"
        elif gmm_model[1] == "CampbellBozorgnia2014":
            nice_gmm_model_str = "Macedo et al. (2021) [Campbell and Bozorgnia (2014)]"
        else:
            nice_gmm_model_str = "Macedo et al. (2021) [{}]".format(gmm_model[1])
    if gmm_model[0] == "LiuMacedo2022SInter":
        if gmm_model[1] == "ParkerEtAl2020SInter":
            nice_gmm_model_str = "Liu and Macedo (2021) [Parker et al. (2020)]"
    if gmm_model[0] == "LiuMacedo2022SSlab":
        if gmm_model[1] == "ParkerEtAl2020SSlab":
            nice_gmm_model_str = "Liu and Macedo (2021) [Parker et al. (2020)]"
else:
    gmm_model_str = gmm_model
    if gmm_model == "SandikkayaAkkar2017Rjb":
        nice_gmm_model_str = "Sandikkaya and Akkar (2017)"
    elif gmm_model == "CampbellBozorgnia2019":
        nice_gmm_model_str = "Campbell and Bozorgnia (2019)"
    elif gmm_model == "LiuMacedo2022M1SSlab":
        nice_gmm_model_str = "Liu and Macedo (2021) (M1, Slab)"
    elif gmm_model == "LiuMacedo2022M1SInter":
        nice_gmm_model_str = "Liu and Macedo (2021) (M1, Inter)"
    elif gmm_model == "LiuMacedo2022M2SSlab":
        nice_gmm_model_str = "Liu and Macedo (2021) (M2, Slab)"
    elif gmm_model == "LiuMacedo2022M2SInter":
        nice_gmm_model_str = "Liu and Macedo (2021) (M2, Inter)"

print("nice_gmm_model_str is " + nice_gmm_model_str)
# Get Data
# How is the geometric mean computed
fullfname1 = test_path + fname1
fullfname2 = test_path + fname2

# Get tectonic region information
tect_region = "CONUS_HW_AK_unique_eqs_tect_region_final.csv"
# tect_region = "CONUS_HW_AK_unique_eqs_tect_region.csv"
fulltec = test_path + tect_region
EQ_w_tec_info = pd.read_csv(fulltec)
# only works for the dataset: ~/Documents/GMM_data/for_kyle/wus_default_metrics_channels(component=h*)_filtered.csv

# fig_dir = "/Users/kksmith/test_figures/"
fig_dir = "/Users/kksmith/figures/GMMvsData/"
ftype = ".png"
tick_fontsize = 10
ax_fontsize = 12
leg_fontsize = 10
title_fontsize = 14

mydata1a = pd.read_csv(fullfname1)
mydata2a = pd.read_csv(fullfname2)
breakpoint()
# Add tectonic region information to mydata1a
mydata1a = mydata1a.merge(
    EQ_w_tec_info[["EarthquakeId", "Tect_region"]], on="EarthquakeId", how="left"
)

mydata2a = mydata2a.merge(
    EQ_w_tec_info[["EarthquakeId", "Tect_region"]], on="EarthquakeId", how="left"
)

mydata1a.sort_values(by=["EventLatitude"], ascending=True, inplace=True)
mydata2a = mydata2a.iloc[mydata1a.index]

bad_stn = "4O.PL01.00.HH"
bad_net = "YB"
mydata1a = mydata1a[mydata1a["StationID"] != bad_stn]
mydata1a = mydata1a[mydata1a["Network"] != bad_net]

# mydata1a["some_ind"] = np.arange(len(mydata1a))
# mydata1a.sort_values(by=['EarthquakeId'])

mydata2a = mydata2a[mydata2a["StationID"] != bad_stn]
mydata2a = mydata2a[mydata2a["Network"] != bad_net]

mydata1a.reset_index(inplace=True, drop=True)
mydata2a.reset_index(inplace=True, drop=True)


use_Globe_lims = True
geo_reg_cols = {}

geo_reg_cols[USname] = "k"
geo_reg_cols["Alaska"] = "g"
geo_reg_cols["Hawaii"] = "m"
geo_reg_cols["WUS"] = "r"
geo_reg_cols["EUS"] = "b"
geo_reg_cols["Global"] = "c"

if dta_ver == "v1":
    depth_plot_lims = [0, 160]
    mag_plot_lims = [4.5, 7]
    plot_dist_range = [0, 450]
    plot_logdist_range = [5 * 10**-1, 5 * 10**2]
    imlims = {"pga": [10**-3, 5 * 10**2], "cav": [10**-4, 10**1]}
elif dta_ver == "v2":
    depth_plot_lims = [-5, 210]
    mag_plot_lims = [3.5, 7.5]
    plot_dist_range = [0, 601]
    plot_logdist_range = [5 * 10**-1, 6 * 10**2]
    imlims = {"pga": [10**-4, 3 * 10**5], "cav": [10**-5, 10**3]}

bias_set = {}
num_stn_set = {}
num_eqs_set = {}
num_recs_set = {}

mydata1a["EarthquakeRake"] = 0
mydata1a["EarthquakeDip"] = 90
mydata1a["EarthquakeZtor"] = mydata1a["EventDepth"].values
mydata1a["Fault_Type"] = "SS"
mydata1a["alpha"] = 90
mydata1a["EarthquakeWidth"] = kaklamanos_relations.calc_width(
    mydata1a[magstr], mydata1a["Fault_Type"]
)

nonSS_tol = 45

# stn_locs = {'lat':mydata1a['StationLatitude'], 'lon': mydata1a['StationLongitude']}

o_start = time.time()
rup_mag_thresh = 7
# fill in information from rupture file
for eqid in np.unique(mydata1a["EarthquakeId"]):

    mycond = mydata1a["EarthquakeId"] == eqid
    mycond_i = [ii for ii, xx in enumerate(mycond) if xx]
    stn_locs = {
        "lat": mydata1a[mycond]["StationLatitude"].values,
        "lon": mydata1a[mycond]["StationLongitude"].values,
    }
    rup_fname = eqid + "_rupture.json"
    full_rup_fname = path_to_rupture_files + rup_fname
    # if (not os.path.isfile(full_rup_fname)) or (
    #     mydata1a.loc[mycond_i[0], "EarthquakeMagnitude"] < rup_mag_thresh
    # ):
    # breakpoint()
    # If file doesn't exist or mag <= rup_mag_thresh
    if (not os.path.isfile(full_rup_fname)) or (
        mydata1a["EarthquakeMagnitude"][mycond_i[0]] <= rup_mag_thresh
    ):
        # assume it is a point source or should I calculate width?
        # print(rup_fname + " does not exist")
        new_rx = np.full_like(stn_locs["lat"], 0)
        new_ry = np.full_like(stn_locs["lon"], 0)
    else:
        print("\n " + rup_fname + " exists!!")
        origin, event, ruptype = get_event_from_rup(full_rup_fname)

        # if ruptype == 'Point':
        # elif ruptype == 'MultiPolygon':

        rjbset, rhypset, rrupset, repiset, o_dis, W_use, ztor_know, dip_know = (
            get_eq_stn_dist_vals(origin, full_rup_fname, stn_locs)
        )
        new_rx = o_dis["rx"]
        new_ry = o_dis["ry"]
        mydata1a.loc[mycond, "EarthquakeRake"] = event["rake"]
        mydata1a.loc[mycond, "EarthquakeDip"] = int(dip_know)
        mydata1a.loc[mycond, "EarthquakeZtor"] = ztor_know
        mydata1a.loc[mycond, "EarthquakeWidth"] = W_use
        if np.abs(event["rake"] - 90) <= nonSS_tol:
            mydata1a.loc[mycond, "Fault_Type"] = "R"
        elif np.abs(event["rake"] + 90) <= nonSS_tol:
            mydata1a.loc[mycond, "Fault_Type"] = "N"
        else:
            mydata1a.loc[mycond, "Fault_Type"] = "SS"

    # fill-in data if it is not provided
    for o_i, mci in enumerate(mycond_i):
        if np.isnan(mydata1a.loc[mci, "GC2_rx"]):
            mydata1a.loc[mci, "GC2_rx"] = new_rx[o_i]
            mydata1a.loc[mci, "GC2_ry"] = new_ry[o_i]
    # W_use = kaklamanos_relations.calc_width(mydata1a[mycond][magstr].values, fault_type)

rx_know = mydata1a["GC2_rx"].values
ry_know = mydata1a["GC2_ry"].values
alpha_assume = kaklamanos_relations.calc_alpha(rx_know, ry_know)
fix_alpha_nan = mytools.find_nparray_nan(alpha_assume)

pre_end = time.time()
print("The time of execution of above program is :", pre_end - o_start, "s")
mydata1a.to_csv(
    test_path
    + "CONUS_HW_AK_h1_records_stn_filtered_wEQparams_M{}_thresh_h1.csv".format(
        rup_mag_thresh
    )
)

mydata2a.to_csv(test_path + "CONUS_HW_AK_h1_records_stn_filtered_h2.csv")

# breakpoint()
alpha_assume_other = 90
alpha_assume[fix_alpha_nan] = alpha_assume_other

# mydata1a["alpha"] = alpha_assume

# for geo_reg in ["Globe"]:
# for geo_reg in ["Globe", "Hawaii"]:
All_tect_sets = ["all_tect"] + list(EQ_w_tec_info["Tect_region"].unique())
All_tect_sets = [tect for tect in All_tect_sets if tect not in ["Volcanic", "Stable"]]
# for tect_reg in EQ_w_tec_info["Tect_region"].unique():
for tect_reg in All_tect_sets:
    print(tect_reg + " \n")
    if tect_reg == "all_tect":
        mydata1b = copy.deepcopy(mydata1a)
        mydata2b = copy.deepcopy(mydata2a)
    else:
        mydata1b = mydata1a[mydata1a["Tect_region"] == tect_reg]
        mydata2b = mydata2a[mydata2a["Tect_region"] == tect_reg]

    regional_cond = {}
    regional_cond["Globe"] = mydata1b["EventLatitude"] > -91
    regional_cond["Alaska"] = mydata1b["EventLatitude"] > 50
    regional_cond["Hawaii"] = mydata1b["EventLatitude"] < 25
    regional_cond["EUS"] = mydata1b["EventLongitude"] >= -104
    regional_cond["WUS"] = np.logical_and(
        mydata1b["EventLongitude"] > -125, mydata1b["EventLongitude"] < -104
    )
    # regional_cond['West_US'] = np.logical_and(25 < mydata1b["EventLatitude"], mydata1b["EventLatitude"] < 50)
    regional_cond[USname] = np.logical_and(
        25 < mydata1b["EventLatitude"], mydata1b["EventLatitude"] < 50
    )

    # for geo_reg in ["Globe", "Alaska", "Hawaii", USname, "WUS", "EUS"]:
    for geo_reg in ["Alaska", "WUS"]:
        print(geo_reg + " \n")
        # data subset by region
        if np.sum(regional_cond[geo_reg]) == 0:
            continue
        mydata1c = mydata1b[regional_cond[geo_reg]]
        mydata2c = mydata2b[regional_cond[geo_reg]]

        if geo_reg == "Globe" or dta_ver == "v2":
            calc_bnd = False
        else:
            calc_bnd = True

        if dta_ver == "v2":
            calc_bnd = True

        z2pt5_assume = 0.6
        z1pt0_assume = 0.1

        fix_nan = mytools.find_dict_nan(mydata1c)
        eqid_know = mydata1c["EarthquakeId"].values
        stnid_know = mydata1c["StationID"].values
        mag_know = mydata1c[magstr].values
        rrup_know = mydata1c["RuptureDistance"].values
        rjb_know = mydata1c["JoynerBooreDistance"].values
        eq_lat = (mydata1c["EventLatitude"],)
        eq_lon = mydata1c["EventLongitude"]
        vs30_know = mydata1c[vs30str].values

        hypo_depth_assume = mydata1c["EventDepth"].values

        # new_rx = kaklamanos_relations.calc_rx2(
        #     rrup_know, rjb_know, mydata1c['EarthquakeDip'], mydata1c['EarthquakeWidth'], mydata1c['EarthquakeZtor'], alpha_assume
        # )
        # breakpoint()

        # investigate why it is nan. Many entries have ztor > rrup, how is it possible?
        # rx_nan_i = mytools.find_nparray_nan(new_rx)
        # for nan_i in rx_nan_i:
        #     print("{6} has rx:{5} for rrup: {0} km, rjb: {1} km, dip: {2},W : {3}, Ztor: {4}".format(rrup_know[nan_i],rjb_know[nan_i],dip_use,W_use[nan_i],ztor_assume[nan_i],new_rx[nan_i],eqid_know[nan_i]))

        # new_rx = kaklamanos_relations.calc_rx(
        #    rrup_know, rjb_know, dip_use, W_use, ztor_assume, alpha_assume
        # )

        # rx_know[fix_nan["GC2_rx"]] = new_rx[fix_nan["GC2_rx"]]
        # rx_know[fix_nan["GC2_rx"]] = 0

        # new_ry = kaklamanos_relations.calc_ry(rx_know, rjb_know, alpha_assume)
        # ry_know[fix_nan["GC2_ry"]] = new_ry[fix_nan["GC2_ry"]]
        # ry_know[fix_nan["GC2_ry"]] = 0

        # alpha_assume[fix_nan["GC2_rx"]] = kaklamanos_relations.calc_alpha(
        #     rx_know[fix_nan["GC2_rx"]], ry_know[fix_nan["GC2_rx"]]
        # )

        print("\n start again")
        mytools.find_dict_nan(mydata1c)

        pga_vals = mydata1c["PGA"].values
        pga_vals2 = mydata2c["PGA"].values
        pga_gm = np.sqrt(pga_vals * pga_vals2)
        mydata_comp1 = {
            "eqid": eqid_know,
            "stnid": stnid_know,
            "mag": mag_know,
            "rjb": rjb_know,
            "repi": mydata1c["EpicentralDistance"].values,
            "rhypo": mydata1c["HypocentralDistance"].values,
            # "rx": rx_know,
            # "ry": ry_know,
            "rx": mydata1c["GC2_rx"].values,
            "ry": mydata1c["GC2_ry"].values,
            "ry0": mydata1c["GC2_ry0"].values,
            # assumptions, use Kaklamanos (2017) to derive?
            "ztor": mydata1c["EarthquakeZtor"],
            "rake": mydata1c["EarthquakeRake"],
            "width": mydata1c["EarthquakeWidth"],
            "rrup": mydata1c["RuptureDistance"].values,
            "cav": mydata1c["CAV"].values,
            "cav2": mydata2c["CAV"].values,
            "dip": mydata1c["EarthquakeDip"],
            "vs30": vs30_know,
            "alpha": mydata1c["alpha"],
            "hypo_depth": hypo_depth_assume,
            "z2pt5": z2pt5_assume,
            "z1pt0": z1pt0_assume,
            "vs30measured": True,
        }

        cav_gm = np.sqrt(mydata_comp1["cav"] * mydata_comp1["cav2"])
        mytools.get_dict_key_len(mydata_comp1)
        # for key_name in mydata_comp1.keys():
        #     print(key_name)
        #     print(len(mydata_comp1[key_name]))
        if isinstance(gmm_model, tuple):
            (
                lnmean_data_cmp,
                sigma_data_cmp,
                tau_data_cmp,
                phi_data_cmp,
                ctx_data_cmp,
            ) = get_cgmm_im.get_cgmm_im(mydata_comp1, gmm_model[0], gmm_model[1], myimt)
        else:
            (
                lnmean_data_cmp,
                sigma_data_cmp,
                tau_data_cmp,
                phi_data_cmp,
                ctx_data_cmp,
            ) = get_gmm_im.get_gmm_im(mydata_comp1, gmm_model, myimt)
        lnmean_data_cmp = convert_GMM_units(gmm_model, lnmean_data_cmp)
        print("tau min = {:0.2f}".format(np.min(tau_data_cmp)))
        print("tau max = {:0.2f}".format(np.max(tau_data_cmp)))
        print("phi min = {:0.2f}".format(np.min(phi_data_cmp)))
        print("phi max= {:0.2f}".format(np.max(phi_data_cmp)))
        # breakpoint()
        mydata1c["tau"] = tau_data_cmp
        mydata1c["phi"] = phi_data_cmp

        # breakpoint()
        # Convert from ln(mean)
        mean_data_cmp = np.exp(lnmean_data_cmp)

        # plt.show()
        # breakpoint()

        ### Investigative plots
        # IM vs Distance
        fig, axes = plt.subplots(2, 2, layout="constrained", figsize=(7.2, 7.2))
        for gm_i, (imstr, IM) in enumerate(zip(["PGA", "CAV"], [pga_gm, cav_gm])):
            # fig, axes = plt.subplots()
            temp_ax = axes[gm_i, 1]
            if geo_reg == "Globe":
                for or_i, geo_reg_alt in enumerate(ind_geo_sec):
                    leg_label = "{:s} (N={:0.0f})".format(
                        geo_reg_alt, regional_cond[geo_reg_alt].sum()
                    )
                    temp_ax.scatter(
                        ctx_data_cmp[dist_meas][regional_cond[geo_reg_alt]],
                        IM[regional_cond[geo_reg_alt]],
                        1,
                        geo_reg_cols[geo_reg_alt],
                        label=leg_label,
                    )
            else:
                temp_ax.scatter(ctx_data_cmp[dist_meas], IM, 1, "k")
            temp_ax.set_xlabel(dist_meas_label + ", km", fontsize=ax_fontsize)
            temp_ax.set_ylabel(
                imstr + ", " + UNITS[imstr.lower()], fontsize=ax_fontsize
            )
            temp_ax.set_yscale("log")
            temp_ax.set_xscale(
                "log"
            )  # good to have for a lot of data that is well spread out
            # axes.axhline(0, ls='--', c='k')
            # temp_ax.set_xlim(plot_logdist_range[0],plot_logdist_range[1])
            # temp_ax.set_ylim(imlims[imstr.lower()][0],imlims[imstr.lower()][1])
            temp_ax.set(
                xlim=(plot_logdist_range[0], plot_logdist_range[1]),
                ylim=(imlims[imstr.lower()][0], imlims[imstr.lower()][1]),
            )
            if gm_i == 1:
                temp_ax.legend(loc="lower left", fontsize=leg_fontsize, fancybox=True)
            # temp_ax.set_title(geo_reg + ": N = " + str(len(cav_gm)), fontsize=title_fontsize)
            # fname = imstr + "_vs_" + dist_meas + "_distance_from_data_" + geo_reg + ftype
            # fullname = fig_dir + fname
            # print(fullname)
            # if myfig_ops["printfigs"]:
            #     plt.savefig(fullname)

        # Mag vs Distance
        # fig, axes = plt.subplots()
        temp_ax = axes[0, 0]

        if geo_reg == "Globe":
            for or_i, geo_reg_alt in enumerate(ind_geo_sec):
                temp_ax.scatter(
                    ctx_data_cmp[dist_meas][regional_cond[geo_reg_alt]],
                    ctx_data_cmp["mag"][regional_cond[geo_reg_alt]],
                    1,
                    geo_reg_cols[geo_reg_alt],
                    label=geo_reg_alt,
                )
        else:
            temp_ax.scatter(ctx_data_cmp[dist_meas], ctx_data_cmp["mag"], 1, "k")
        temp_ax.set_xlabel(dist_meas_label + ", km", fontsize=ax_fontsize)
        temp_ax.set_ylabel("M", fontsize=ax_fontsize)
        temp_ax.set_xscale("log")
        temp_ax.set_ylim(mag_plot_lims)
        # temp_ax.set_xlim(plot_dist_range)
        temp_ax.set_xlim(plot_logdist_range)
        # temp_ax.set_title(geo_reg + ": N = " + str(len(ctx_data_cmp["mag"])), fontsize=title_fontsize)
        # fname = "M_vs_" + dist_meas + "_distance_from_data_" + geo_reg + ftype
        # fullname = fig_dir + fname
        # print(fullname)
        # if myfig_ops["printfigs"]:
        #     plt.savefig(fullname)

        # CAV vs PGA
        # fig, axes = plt.subplots()
        temp_ax = axes[1, 0]
        if geo_reg == "Globe":
            for or_i, geo_reg_alt in enumerate(ind_geo_sec):
                temp_ax.scatter(
                    pga_gm[regional_cond[geo_reg_alt]],
                    cav_gm[regional_cond[geo_reg_alt]],
                    1,
                    geo_reg_cols[geo_reg_alt],
                    label=geo_reg_alt,
                )
        else:
            temp_ax.scatter(pga_gm, cav_gm, 1, "k")
        temp_ax.set_xlabel("PGA, " + UNITS["pga"], fontsize=ax_fontsize)
        temp_ax.set_ylabel("CAV, " + UNITS["cav"], fontsize=ax_fontsize)
        temp_ax.set_yscale("log")
        temp_ax.set_xscale(
            "log"
        )  # good to have for a lot of data that is well spread out
        # temp_ax.axhline(0, ls='--', c='k')
        temp_ax.set(
            xlim=(imlims["pga"][0], imlims["pga"][1]),
            ylim=(imlims["cav"][0], imlims["cav"][1]),
        )
        # temp_ax.set_title(geo_reg + ": N = " + str(len(cav_gm)), fontsize=title_fontsize)
        # fname = "CAV_vs_PGA_from_data_" + geo_reg +  ftype
        # fullname = fig_dir + fname
        # print(fullname)
        # if myfig_ops["printfigs"]:
        #     plt.savefig(fullname)
        num_recs = len(cav_gm)
        fig.suptitle(
            geo_reg + "-" + tect_reg + ": Number of Records = " + str(num_recs),
            fontsize=title_fontsize,
        )
        num_recs_set[geo_reg] = num_recs

        fname = "Investigative_plots_" + geo_reg + "_" + tect_reg + ftype
        fullname = fig_dir + fname
        print("printing: " + fullname)
        if myfig_ops["printfigs"]:
            plt.savefig(fullname, dpi=dpi_res)
        # if tect_reg == "all_tect":
        #     break
        eq_numstns_set = []
        unique_eqs_filt = []
        eq_I_set = []

        # Plot GMM line with data
        n = 100
        dist_pts = np.logspace(plot_dist_range[0], np.log10(plot_dist_range[1]), n)
        # num_eqs = len(unique_eqs)
        # print("Number of EQ IDs: {} ".format(num_eqs))
        # unique_eqs = np.unique(mydata1d["EarthquakeId"])
        unique_eqs = np.unique(mydata1c["EarthquakeId"])
        # num_eqs_set[geo_reg] = num_eqs
        for eqid in unique_eqs:
            # get EQ parameters and do distance scaling of IM
            myeq = mydata_comp1["eqid"] == eqid
            # eq_numstns_set.append(sum(myeq))
            # eq_I_set.append(np.where(myeq)[0][0])
            mydata_comp2 = {}
            for key_name, mydata_val in mydata_comp1.items():
                # print(type(mydata_val))
                # print(key_name)
                # if key_name == "vs30" or key_name == "rx":
                #     breakpoint()
                # breakpoint()
                try:
                    mydata_comp2[key_name] = mydata_val[myeq][0]
                except TypeError:
                    try:
                        mydata_comp2[key_name] = mydata_val[myeq]
                    except TypeError:
                        mydata_comp2[key_name] = mydata_val
                except KeyError:
                    mydata_comp2[key_name] = mydata_val.values[myeq][0]

            mydata_comp2["rjb"] = dist_pts
            mydata_comp2["rx"] = (
                mydata_comp2["width"] * np.cos(np.radians(mydata_comp2["dip"]))
                + mydata_comp2["rjb"]
            )
            mydata_comp2["ry"] = kaklamanos_relations.calc_ry(
                mydata_comp2["rx"], mydata_comp2["rjb"], mydata_comp2["rjb"]
            )
            mydata_comp2["rrup"] = kaklamanos_relations.calc_rrup(
                mydata_comp2["rx"],
                mydata_comp2["ry"],
                mydata_comp2["rjb"],
                mydata_comp2["dip"],
                mydata_comp2["width"],
                mydata_comp2["ztor"],
            )

            if isinstance(gmm_model, tuple):
                lnmean_est, sigma_est, tau_est, phi_est, ctx_est = (
                    get_cgmm_im.get_cgmm_im(
                        mydata_comp2, gmm_model[0], gmm_model[1], myimt
                    )
                )
            else:
                lnmean_est, sigma_est, tau_est, phi_est, ctx_est = (
                    get_gmm_im.get_gmm_im(mydata_comp2, gmm_model, myimt)
                )
            mean_est = np.exp(lnmean_est)
            # print(eqid)

            if str(myimt) == "IA":
                unit_label = UNITS["ariasintensity"]
            else:
                unit_label = UNITS[str(myimt).lower()]
            add_unit_label = ", " + unit_label

            if myfig_ops["makefigs"] and False:
                fig, ax = plt.subplots()
                # ax.loglog(
                #     ctx_data_cmp[dist_meas][myeq], mean_data_cmp[myeq], "ro", label="estpts"
                # )  # should be this
                ax.loglog(mydata_comp2[dist_meas], mean_est, "r-", label="est")
                ax.loglog(
                    ctx_data_cmp[dist_meas][myeq], cav_gm[myeq], "ko", label="data"
                )

                title_str = (
                    geo_reg
                    + ": "
                    + eqid
                    + ", M="
                    + str(ctx_data_cmp["mag"][myeq][0])
                    + ", Vs30="
                    + str(ctx_data_cmp["vs30"][myeq][0])
                )

                ax.set_title(
                    title_str,
                    fontsize=title_fontsize,
                )
                plt.xlabel(dist_meas_label + " km", fontsize=ax_fontsize)
                plt.ylabel(str(myimt) + add_unit_label, fontsize=ax_fontsize)
                # plt.show()
                # breakpoint()
                print("how is it?")

        # Mixed Effects regression on residuals
        # from GM process tutorial
        # plot withing-event terms and between-terms
        res = cav_gm - mean_data_cmp
        # res_pct_diff = res / cav_gm * 100

        # residual in log space
        # print(geo_reg)
        resid_col_name = str(myimt) + "residuals"
        # resid_col_name = "residuals"
        # breakpoint()
        mydata1c[resid_col_name] = np.log(cav_gm) - np.log(mean_data_cmp)

        mydata1c["zscore"] = (np.log(cav_gm) - np.log(mean_data_cmp)) / sigma_data_cmp
        # mydata1c[resid_col_name] = res
        mydata1c.replace([np.inf, -np.inf], np.nan, inplace=True)
        non_nan_i = ~mydata1c[resid_col_name].isna()
        mydata1d = mydata1c[non_nan_i]
        # breakpoint()

        # plot this laters

        num_eqs = len(unique_eqs)
        print("Number of EQ IDs: {} ".format(num_eqs))
        num_eqs_set[geo_reg] = num_eqs

        num_recs = []
        for UE in unique_eqs:
            num_recs.append(np.sum(mydata1d["EarthquakeId"] == UE))
        NumRecs = pd.DataFrame(
            np.transpose([unique_eqs, num_recs]), columns=["EarthquakeId", "num_recs"]
        )

        infile_eid = "EarthquakeId"
        mydata1e = mydata1d.join(NumRecs.set_index(infile_eid), on=infile_eid)
        # [unique_eqs,num_recs]

        # final_vals = pd.DataFrame(np.transpose([cav_gm[non_nan_i],mean_data_cmp[non_nan_i],sigma_data_cmp[non_nan_i]]))
        # data stored in columns with columns corresponding to (1) observations (2) predictions (3) sigma (4) tau (5) phi
        if makefiles:
            full_final_vals = pd.DataFrame(
                np.transpose(
                    [
                        np.log(cav_gm[non_nan_i]),
                        np.log(mean_data_cmp[non_nan_i]),
                        sigma_data_cmp[non_nan_i],
                        tau_data_cmp[non_nan_i],
                        phi_data_cmp[non_nan_i],
                        mydata1e["num_recs"],
                    ]
                )
            )
            # breakpoint()
            # For statistical tests
            full_final_vals.to_csv(
                out_path_0
                + "Full_GMMvsData_{}_{}_{}_{}.txt".format(
                    gmm_model_str, geo_reg, tect_reg, dta_ver
                ),
                sep="\t",
                index=False,
                header=False,
            )
            NumRecs.to_csv(
                out_path_0
                + "NumRecs_{}_{}_{}_{}.txt".format(
                    gmm_model_str, geo_reg, tect_reg, dta_ver
                ),
                sep="\t",
                index=False,
                header=False,
            )
            mydata1d.to_csv(
                test_path
                + "CONUS_HW_AK_h1_records_wRes_{}_{}_{}_{}.csv".format(
                    gmm_model_str, geo_reg, tect_reg, dta_ver
                ),
            )

        # breakpoint()

        for eqid in unique_eqs:
            # print(eqid)
            # get filtered EQ parameter
            myeq = mydata1d["EarthquakeId"] == eqid
            eq_numstns = sum(myeq)
            if eq_numstns == 0:
                continue
            unique_eqs_filt.append(eqid)
            eq_numstns_set.append(eq_numstns)
            eq_I_set.append(np.where(myeq)[0][0])
        # mydata1c.dropna(subset=['residuals'],inplace=True)
        #
        mdf = smf.mixedlm(
            "%s ~ 1" % resid_col_name, mydata1d, groups=mydata1d["EarthquakeId"]
        ).fit()
        print(mdf.summary())
        bias = mdf.fe_params.Intercept
        print("bias = {:0.4f}".format(bias))
        # breakpoint()
        bias_set[geo_reg] = bias
        num_wth = len(mdf.resid)
        # breakpoint()

        if calc_bnd is False:
            if geo_reg == "Globe":
                if isinstance(gmm_model, tuple):
                    if gmm_model[0] == "MacedoEtAl2021":
                        if gmm_model[1] == "BooreEtAl2014":
                            if dta_ver == "v1":
                                btw_bnd = 1
                                wth_bnd = 3.5
                                zbx_bnd = 3
                                zwx_bnd = 7
                        if gmm_model[1] == "CampbellBozorgnia2014":
                            btw_bnd = 1.25
                            wth_bnd = 3.5
                            zbx_bnd = 4.5
                            zwx_bnd = 8
                else:
                    if gmm_model == "CampbellBozorgnia2019":
                        btw_bnd = 1.25
                        wth_bnd = 3.5
                        zbx_bnd = 4
                        zwx_bnd = 8
                    elif gmm_model == "SandikkayaAkkar2017Rjb":
                        if dta_ver == "v1":
                            btw_bnd = 1.25
                            wth_bnd = 4
                            zbx_bnd = 3.5
                            zwx_bnd = 6
                        # elif region == "CONUS":
                        #     btw_bnd = 1.25
                        #     wth_bnd = 4
                        #     zbx_bnd = 3.5
                        #     zwx_bnd = 6
                try:
                    zbx = np.linspace(-zbx_bnd, zbx_bnd, n)
                    zwx = np.linspace(-zwx_bnd, zwx_bnd, n)
                except NameError:
                    print(
                        "NameError because no custom bounds exist for this model(s) yet, try calc_bnd = True, then you can make set good boundaries"
                    )
                    sys.exit(1)
            else:
                print(
                    "NameError because no custom bounds exist for this region/model(s) yet, try calc_bnd = True, then you can make set good boundaries"
                )
                sys.exit(1)
        # re_array = [float(re.iloc[0]) for group, re in mdf.random_effects.items()]
        # fig = plt.figure()
        # plt.plot(re_array, "ko")
        # plt.xlabel("Index",fontsize=ax_fontsize)
        # plt.ylabel("Random effect",fontsize=ax_fontsize)
        # plt.ylim(btw_plot_val_lims)
        # plt.title(gmm_model_str + ", N = " + str(len(re_array)), fontsize=title_fontsize)
        bnd_pass_prop = 1.10
        if calc_bnd:
            wth_bnd = (
                np.max([abs(np.min(mdf.resid)), abs(np.max(mdf.resid))]) * bnd_pass_prop
            )
        # Within-event terms
        fig, axes = plt.subplots(4, 1, figsize=(5.5, 8), layout="constrained")
        for x_i, stn_met in enumerate([rjb_know[non_nan_i], vs30_know[non_nan_i]]):
            temp_ax = axes[x_i]
            # breakpoint()
            temp_ax.plot(stn_met, mdf.resid, "go", markersize=1)
            if x_i == 0:
                temp_ax.set_xlabel(dist_meas_label + ", km", fontsize=ax_fontsize)
                temp_ax.set_xlim(plot_logdist_range)
                temp_ax.set_xscale("log")
                temp_ax.set_title(
                    geo_reg + " - " + tect_reg + ": " + nice_gmm_model_str,
                    fontsize=title_fontsize,
                )
            elif x_i == 1:
                temp_ax.set_xlim([0, 2500])
                temp_ax.set_xlabel("Vs30, m/s", fontsize=ax_fontsize)
            temp_ax.set_ylabel("Within-Event Terms", fontsize=ax_fontsize)
            temp_ax.set_ylim([-wth_bnd, wth_bnd])
            temp_ax.axhline(0, ls="--", c="grey")
        # temp_ax.text(,,"N = " + str(num_wth))
        # fname = "Within-Event_Terms_vs_" + dist_meas + "_distance_" + gmm_model_str + "_" + geo_reg + "_" + tect_reg + ftype
        # fullname = fig_dir + fname
        # print(fullname)
        # breakpoint()
        # if myfig_ops["printfigs"]:
        #     plt.savefig(fullname)

        # Between-event terms?
        btw_event_terms = pd.DataFrame(mdf.random_effects).T
        mydata1d = mydata1d.merge(
            btw_event_terms, left_on="EarthquakeId", right_index=True
        )
        df_events = mydata1d.drop_duplicates(subset=["EarthquakeId"])

        num_btw = len(df_events["Group"])
        if calc_bnd:
            btw_bnd = (
                np.max(
                    [abs(np.min(df_events["Group"])), abs(np.max(df_events["Group"]))]
                )
                * bnd_pass_prop
            )

        proper_labels = ["Hypocentral Depth, km", "M"]
        temp_xlims = [depth_plot_lims, mag_plot_lims]

        for eqpar_i, data_labels in enumerate(["EventDepth", magstr]):
            # fig, axes = plt.subplots()
            temp_ax = axes[eqpar_i + 2]
            temp_ax.scatter(df_events[data_labels], df_events["Group"], color="g")
            temp_ax.set_xlabel(proper_labels[eqpar_i], fontsize=ax_fontsize)
            temp_ax.set_xlim(temp_xlims[eqpar_i])
            temp_ax.set_ylabel("Between-Event Terms", fontsize=ax_fontsize)
            temp_ax.set_ylim([-btw_bnd, btw_bnd])
            temp_ax.axhline(0, ls="--", c="grey")
            # temp_ax.set_title(geo_reg + ": " + nice_gmm_model_str + ", N = " + str(num_btw), fontsize=title_fontsize)
            # plt.tight_layout()
            # fname = "Between-Event_terms_vs_" + data_labels + "_" + gmm_model_str + "_" + geo_reg + "_" + tect_reg + ftype
            # fullname = fig_dir + fname
            # print(fullname)
            # if myfig_ops["printfigs"]:
            #     plt.savefig(fullname)

        fname = (
            "Exploration_of_Wtn-Event_and_Btw-Event_terms_"
            + gmm_model_str
            + "_"
            + geo_reg
            + "_"
            + tect_reg
            + ftype
        )
        fullname = fig_dir + fname
        print("printing: " + fullname)
        if myfig_ops["printfigs"]:
            plt.savefig(fullname, dpi=dpi_res)
        # Histograms of residuals
        num_bins = [20, 20]
        res_plt_lims = [[-btw_bnd, btw_bnd], [-wth_bnd, wth_bnd]]
        res_labels = ["Between-Event Terms", "Within-Event Terms"]

        # plt.show()
        # breakpoint()

        res_file_labels = ["Between-Event", "Within-Event"]
        for res_i, residuals in enumerate([df_events["Group"], mdf.resid]):
            fig, axs = plt.subplots(1, 1, tight_layout=True)
            stat_label = "mean={:0.2f}, st. dev={:0.2f}".format(
                np.mean(residuals), np.std(residuals)
            )
            axs.hist(
                residuals,
                bins=num_bins[res_i],
                density=True,
                color="k",
                label=stat_label,
            )
            data_res_std = np.std(residuals)
            print("{} std dev={:0.4f}".format(res_labels[res_i], data_res_std))
            if res_i == 0:
                data_res_mean = np.mean(residuals)
            else:
                data_res_mean = 0

            # xmin, xmax = axs.get_xlim()
            # absmax = np.max([np.abs(xmin), np.abs(xmax)])
            # axs.set_xlim(-absmax,absmax)
            # x = np.linspace(-absmax,absmax,100)
            axs.set_xlim(res_plt_lims[res_i])
            x = np.linspace(res_plt_lims[res_i][0], res_plt_lims[res_i][1], n)
            # axs.plot(x, norm.pdf(x,data_res_mean,data_res_std), 'r', linewidth = 2.5)
            axs.set_title(
                "{:s} - {:s}: {:s}\n {:s}\n{:s}, N = {:0.0f}".format(
                    geo_reg,
                    tect_reg,
                    nice_gmm_model_str,
                    res_labels[res_i],
                    stat_label,
                    len(residuals),
                ),
                fontsize=title_fontsize,
            )

            # plt.legend(loc="upper left")
            # data_res_mean,data_res_std,
            # mean = {:0.2f}, std = {:0.2f},
            fname = (
                res_file_labels[res_i]
                + "_terms_hist"
                + "_"
                + gmm_model_str
                + "_"
                + geo_reg
                + "_"
                + tect_reg
                + ftype
            )
            fullname = fig_dir + fname
            print("printing: " + fullname)
            if myfig_ops["printfigs"]:
                plt.savefig(fullname, dpi=dpi_res)

        # Plot Standardized Version
        x = np.linspace(-3, 3, n)
        Z_B = df_events["Group"] / df_events["tau"]
        Z_W = mdf.resid / mydata1d["phi"]

        if calc_bnd:
            zbx_bnd = np.max([np.abs(np.min(Z_B)), np.abs(np.max(Z_B))])
            zwx_bnd = np.max([np.abs(np.min(Z_W)), np.abs(np.max(Z_W))])
            zbx = np.linspace(-zbx_bnd, zbx_bnd, n)
            zwx = np.linspace(-zwx_bnd, zwx_bnd, n)

        fig, axs = plt.subplots(1, 1, tight_layout=True)
        leg_label = "mean={:0.2f}\nst. dev={:0.2f}".format(np.mean(Z_W), np.std(Z_W))
        axs.hist(Z_W, bins=20, density=True, color="g", label=leg_label)
        axs.plot(zwx, norm.pdf(zwx), "b", linewidth=2.5)
        axs.set_xlim(-zwx_bnd, zwx_bnd)
        axs.set_title(
            "{:s} - {:s}: {:s}\n Standardized Within-Event Terms, N = {:0.0f}".format(
                geo_reg, tect_reg, nice_gmm_model_str, num_wth
            )
        )
        plt.legend(loc="upper left")
        fname = (
            "Std_Within-Event_terms_hist"
            + "_"
            + gmm_model_str
            + "_"
            + geo_reg
            + "_"
            + tect_reg
            + ftype
        )
        fullname = fig_dir + fname
        print("printing: " + fullname)
        if myfig_ops["printfigs"]:
            plt.savefig(fullname, dpi=dpi_res)

        fig, axs = plt.subplots(1, 1, tight_layout=True)
        leg_label = "mean={:0.2f}\nst. dev={:0.2f}".format(np.mean(Z_B), np.std(Z_B))
        axs.hist(Z_B, bins=20, density=True, color="g", label=leg_label)
        axs.plot(zbx, norm.pdf(zbx), "b", linewidth=2.5)
        axs.set_xlim(-zbx_bnd, zbx_bnd)
        axs.set_title(
            "{:s} - {:s}: {:s}\n Standardized Between-Event Terms, N = {:0.0f}".format(
                geo_reg, tect_reg, nice_gmm_model_str, num_btw
            )
        )
        # axs.plot([], [], ' ', label=leg_label)
        plt.legend(loc="upper left")
        #
        fname = (
            "Std_Between-Event_terms_hist"
            + "_"
            + gmm_model_str
            + "_"
            + geo_reg
            + "_"
            + tect_reg
            + ftype
        )
        fullname = fig_dir + fname
        print("printing: " + fullname)
        if myfig_ops["printfigs"]:
            plt.savefig(fullname, dpi=dpi_res)

        # fig, axs = plt.subplots(1, 1, tight_layout=True)
        # axs.hist(mydata1d["zscore"], bins=20,density=True,color='b')
        # zscore_avg = np.mean(mydata1d["zscore"])
        # zscore_std = np.std(mydata1d["zscore"])
        # axs.set_title("{:s}: {:s}\n Zscore (data_mean-$\mu$)/$\sigma$ mean = {:0.2f}, std = {:0.2f}, N = {:0.0f}".format(geo_reg, nice_gmm_model_str, zscore_avg,zscore_std,len(mydata1d["zscore"])), fontsize=title_fontsize)
        # axs.plot(x, norm.pdf(x), 'r', linewidth = 2.5)
        # axs.set_title("{:s}: {:s}, mean = {:0.2f}, std = {:0.2f}, N = {:0.0f}".format(geo_reg, res_labels[res_i],data_res_mean,data_res_std,len(residuals)), fontsize=title_fontsize)

        # Make EQ CSV file, later add in average between event values
        deq = {
            # "eqid": unique_eqs,
            "eqid": unique_eqs_filt,
            "eqlat": mydata1d["EventLatitude"].values[eq_I_set],
            "eqlon": mydata1d["EventLongitude"].values[eq_I_set],
            "numstns": eq_numstns_set,
            "mag": mydata1d[magstr].values[eq_I_set],
            "depth": mydata1d["EventDepth"].values[eq_I_set],
        }

        deq_copy = copy.deepcopy(deq)

        dfeq = pd.DataFrame(data=deq_copy)

        if makefiles:
            dfeq.to_csv(
                out_path
                + "EQ_"
                + gmm_model_str
                + "_"
                + geo_reg
                + "_"
                + tect_reg
                + "_"
                + dta_ver
                + ".csv",
                index=False,
            )

        # Print Within-event terms to station CSV file
        unique_stations = np.unique(mydata1d["StationID"])
        num_stn_ids = len(unique_stations)
        print("Number of station IDs: {}".format(num_stn_ids))
        num_stn_set[geo_reg] = num_stn_ids
        stn_numeqs = []
        stn_avg = []
        stn_id_set = []
        for stn in unique_stations:
            stn_id_TF_set = stn == mydata1d["StationID"]
            # stn_id_TF_set = stn == mydata_comp1["stnid"]
            stn_numeqs.append(sum(stn_id_TF_set))
            stn_avg.append(np.mean(mdf.resid[stn_id_TF_set]))
            stn_id_set.append(np.where(stn_id_TF_set)[0][0])

        dstn = {
            "stnid": unique_stations,
            "stnlat": mydata1d["StationLatitude"].values[stn_id_set],
            "stnlon": mydata1d["StationLongitude"].values[stn_id_set],
            "numeqs": stn_numeqs,
            "dS2S": stn_avg,
        }
        dfstn = pd.DataFrame(data=dstn)
        if makefiles:
            dfstn.to_csv(
                out_path
                + "STN_"
                + gmm_model_str
                + "_"
                + geo_reg
                + "_"
                + tect_reg
                + "_"
                + dta_ver
                + ".csv",
                index=False,
            )

        print("end of region: " + geo_reg)
        # if myfig_ops["showfigs"]:
        #     plt.show()
        # breakpoint()

print("bias_set:")
print(bias_set)

print("num_eq_set:")
print(num_eqs_set)

print("num_stn_set:")
print(num_stn_set)

print("num_recs_set:")
print(num_recs_set)
if myfig_ops["showfigs"]:
    plt.show()

    # record end time
end = time.time()

# print the difference between start
# and end time
print("The time of execution of above program is :", end - start, "s")
print("end")
