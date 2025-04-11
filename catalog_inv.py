import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
import time
from gmprocess.utils.constants import UNITS, STATION_METRIC_UNITS
import matplotlib as mpl
import copy

# make files of unique eqs and unique stns of weird data
make_txt_files = False
bprint = False

start = time.time()

test_path = "/Users/kksmith/Documents/GMM_data/for_kyle_conus/"
# fname1 = "CONUS_HW_AK_h1_records.csv"
fname1 = "CONUS_HW_AK_h1_records_wRes_SandikkayaAkkar2017Rjb_Globe_v2.csv"
fname2 = "dsgmd_h2.csv"
fullfname1 = test_path + fname1
fullfname2 = test_path + fname2

# only works for the dataset: /Users/kksmith/Documents/GMM_data/for_kyle/wus_default_metrics_channels(component=h*)_filtered.csv

fig_dir = "/Users/kksmith/figures/debug/"

ftype = ".png"
dpi_res = 250
tick_fontsize = 10
ax_fontsize = 12
leg_fontsize = 10
title_fontsize = 14


# 10**-1 < mydata1a[mydata1a['CAV'] < 10 ** 3]
# mydata1a[mydata1a['CAV'] > 10 ** -1]

# filtering
datasets = {}
# filtering out by station
bad_stn = "4O.PL01.00.HH"
bad_net = "YB"
# filtering by subset
# subsets = ["original", "almost_original", "dur_pga_small_blob"]
# subsets = ["original", "almost_original", "pga_cav_blob", "dur_pga_small_blob", "pga_cav_small"]
# subsets = ["original", "almost_original", "dur_pga_small_blob"]
subsets = ["original"]
for datasubset in subsets:
    if datasubset == "original":
        mydata1a = pd.read_csv(fullfname1)
        mydata2a = pd.read_csv(fullfname2)
        datasets[datasubset] = copy.copy(mydata1a)
        continue
    elif datasubset == "almost_original":
        mydata1a = pd.read_csv(fullfname1)
        mydata2a = pd.read_csv(fullfname2)

        mydata1ba = mydata1a[mydata1a["StationID"] != bad_stn]
        mydata1b = mydata1ba[mydata1ba["Network"] != bad_net]
        # mydata1b = copy.copy(mydata1ba)
        datasets[datasubset] = copy.copy(mydata1b)
        continue
    elif datasubset == "pga_cav_blob":
        mydata1a = pd.read_csv(fullfname1)
        mydata2a = pd.read_csv(fullfname2)
        # PGA vs CAV
        mydata1ba = mydata1a[mydata1a["StationID"] != bad_stn]
        mydata1b = mydata1ba[mydata1ba["Network"] != bad_net]
        weird_data = mydata1b[
            (10**-1 < mydata1b["CAV"])
            & (mydata1b["CAV"] < 10**3)
            & (mydata1b["PGA"] < 10**-1)
        ]
        datasets[datasubset] = weird_data
        # datasets[datasubset] = copy.copy(mydata1a)
        continue
    elif datasubset == "dur_pga_small_blob":
        # PGA vs Dur
        mydata1ba = mydata1a[mydata1a["StationID"] != bad_stn]
        mydata1b = mydata1ba[mydata1ba["Network"] != bad_net]
        weird_data = mydata1b[
            (mydata1b["PGA"] < 10**-5) & (mydata1b["Duration(5-95)"] < 10**-5)
        ]
        datasets[datasubset] = weird_data
    elif datasubset == "pga_cav_small":
        # PGA vs CAV
        weird_data = mydata1a[(mydata1a["CAV"] < 10**-16) & (mydata1a["PGA"] < 10**-16)]
        datasets[datasubset] = weird_data

    # making files of weird records
    if make_txt_files:
        # records for EQ and station
        weird_data.to_csv(
            "~/output_conus/weird_records_{}.txt".format(datasubset),
            index=False,
            sep="\t",
        )

        unique_eqs = np.unique(weird_data["EarthquakeId"])
        unique_stns = np.unique(weird_data["StationID"])

        orig_num_recs_per_eq = []
        for ue in unique_eqs:
            orig_num_recs_per_eq.append(sum(mydata1a["EarthquakeId"] == ue))

        orig_num_recs_per_stn = []
        for us in unique_stns:
            orig_num_recs_per_stn.append(sum(mydata1a["StationID"] == us))

        end = time.time()
        print("Time of execution of above program is :", end - start, "s")
        # weird records by EQ
        num_recs_per_eq = []
        for ue in unique_eqs:
            num_recs_per_eq.append(sum(weird_data["EarthquakeId"] == ue))
        UE = pd.DataFrame(
            np.transpose(
                [
                    unique_eqs,
                    num_recs_per_eq,
                    orig_num_recs_per_eq,
                    np.array(num_recs_per_eq) / np.array(orig_num_recs_per_eq) * 100,
                ]
            ),
            columns=[
                "EarthquakeId",
                "Number_of_weird_records",
                "Total_number_of_records",
                "%_Weird_records",
            ],
        )
        # breakpoint()
        UE.to_csv(
            "~/output_conus/EQs_of_weird_records_{}.txt".format(datasubset),
            index=False,
            sep="\t",
        )

        # weird records by stn
        num_recs_per_stn = []
        for us in unique_stns:
            num_recs_per_stn.append(sum(weird_data["StationID"] == us))
        US = pd.DataFrame(
            np.transpose(
                [
                    unique_stns,
                    num_recs_per_stn,
                    orig_num_recs_per_stn,
                    np.array(num_recs_per_stn) / np.array(orig_num_recs_per_stn) * 100,
                ]
            ),
            columns=[
                "StationID",
                "Number_of_weird_records",
                "Total_number_of_records",
                "%_Weird_records",
            ],
        )
        US.to_csv(
            "~/output_conus/STNs_of_weird_records_{}.txt".format(datasubset),
            index=False,
            sep="\t",
        )

        # pga filt
        # mydata1a['PGA'] < 10 ** -1

#############################################################
# initial investigate

#
# low_pga_dur = np.logical_and(mydata1a["PGA"]<1, mydata1a["Duration(5-95)"] < 10)
# #low_pga_dur_hi_cav = np.logical_and(low_pga_dur,np.log(mydata1a["CAV"])>5)

# weird plot
# plt.scatter(mydata1a['PGA'], mydata1a['CAV'], c=mag_cols, cmap='plasma')
# plt.scatter(mydata1a['PGA'], mydata1a['CAV'], c=mydata1a['PreferredMagnitude'], cmap='plasma',label="Mag")
# mydata1a.dropna(subset=['PreferredMagnitude'])
# Inv_Val = "CAV"
Inv_Val = "residuals"

if Inv_Val == "CAV":
    plt_scale = "log"
elif Inv_Val == "residuals":
    plt_scale = "linear"

for datasubset in subsets:
    eqval = Inv_Val
    print(np.min(np.log10(datasets[datasubset][eqval])))


# Investigative Plots
for datasubset in subsets:
    eqval = Inv_Val
    if Inv_Val == "residuals":
        cav_dtaset = datasets[datasubset][eqval]
        xlab = eqval
    else:
        cav_dtaset = np.log10(datasets[datasubset][eqval])
        xlab = "log10({})".format(eqval)
    fig, axs = plt.subplots(1, 1, tight_layout=True)
    stat_label = "mean={:0.2f}, st. dev={:0.2f}".format(
        np.mean(cav_dtaset), np.std(cav_dtaset)
    )
    axs.hist(cav_dtaset, bins=20, density=False, color="k", label=stat_label)
    axs.set_title(datasubset)
    axs.set_xlabel(xlab)

scatterplotsets = [
    ["PGA", "Duration(5-95)", Inv_Val],
    ["PGA", Inv_Val, "RuptureDistance"],
    ["PGA", Inv_Val, "EarthquakeMagnitude"],
    ["RuptureDistance", Inv_Val, "EarthquakeMagnitude"],
    ["RuptureDistance", "Duration(5-95)", "EarthquakeMagnitude"],
    ["Duration(5-95)", Inv_Val, "RuptureDistance"],
    ["EarthquakeMagnitude", "Duration(5-95)", Inv_Val],
]

ctr = -1
set_std_xlim = {}
for datasubset in subsets:
    ctr += 1

    if datasubset == "pga_cav_blob" and ctr == 0:
        breakpoint()

    somedata = datasets[datasubset]
    for xydata in scatterplotsets:
        fig, axes = plt.subplots()
        if "EarthquakeMagnitude" == xydata[2]:
            mycmap = "tab20b"
        elif "RuptureDistance" == xydata[2]:
            mycmap = "plasma"
        elif Inv_Val == xydata[2]:
            mycmap = "viridis"
        else:
            mycmap = "turbo"

        if xydata[2] == "residuals":
            line = axes.scatter(
                somedata[xydata[0]],
                somedata[xydata[1]],
                s=1,
                c=somedata[xydata[2]],
                cmap=mycmap,
            )
            fig.colorbar(line, ax=axes, label=xydata[2])
        # if xydata[2] == Inv_Val and Inv_Val == "CAV":
        # if xydata[2] == "CAV":

        else:
            line = axes.scatter(
                somedata[xydata[0]],
                somedata[xydata[1]],
                s=1,
                c=np.log10(somedata[xydata[2]]),
                cmap=mycmap,
            )
            fig.colorbar(line, ax=axes, label="log10({})".format(xydata[2]))
        if xydata[1] == "residuals":
            myxmin, myxmax = axes.get_xlim()
            max_y = np.max(np.abs(axes.get_ylim()))
            plot_max_y = np.ceil(max_y * 1.05)
            axes.set_ylim(-plot_max_y, plot_max_y)
            axes.hlines(0, myxmin, myxmax, colors="gray", linestyles="dashed")
        # else:

        magloc = np.where("EarthquakeMagnitude" in xydata)
        if len(magloc) > 0 or magloc == 2:
            axes.set_xscale("log")
            axes.set_yscale(plt_scale)
        elif magloc == 0:
            axes.set_yscale(plt_scale)
        elif magloc == 1:
            axes.set_xscale("log")
        # breakpoint()
        # plt.scatter(mydata1a['Duration(5-95)'], mydata1a['CAV'], c=mydata1a['PreferredMagnitude'], cmap='viridis')
        # plt.scatter(mydata1a['PGA'], mydata1a['CAV'], c=mydata1a['RuptureDistance'], cmap='viridis')
        # plt.scatter(mydata1_cav_out['PGA'], mydata1_cav_out['Duration(5-95)'], s=1, c=np.log(mydata1_cav_out['CAV']), cmap='viridis',vmin=-10, vmax=10)
        # plt.scatter(mydata1a['PGA'], mydata1a['Duration(5-95)'], s=1, c=np.log10(mydata1a['CAV']), cmap='viridis')
        # plt.scatter(mydata1a['PGA'], mydata1a['CAV'], c=mydata1a['Duration(5-95)'], cmap='turbo')
        # plt.scatter(mydata1a['PGA'], mydata1a['CAV'], s=1, c=mydata1a['Duration(5-95)'], cmap='viridis')
        # plt.scatter(mydata1a['PGA'], mydata1a['CAV'], c=mydata1a['EarthquakeMagnitude'], cmap='tab20b')
        # plt.scatter(mydata1a['RuptureDistance'], mydata1a['CAV'], s=1, c=mydata1a['EarthquakeMagnitude'], cmap='tab20b')
        axes.set_xlabel(xydata[0], fontsize=ax_fontsize)
        axes.set_ylabel(xydata[1], fontsize=ax_fontsize)
        axes.set_title(datasubset, fontsize=ax_fontsize + 2)
        axes.tick_params(axis="both", top=True, right=True)
        use_almost_original_axes = False
        plotcombo = "_".join(xydata)
        if use_almost_original_axes:
            if datasubset != "original" or datasubset != "pga_cav_small":
                if datasubset == "almost_original":
                    set_std_xlim[plotcombo] = axes.get_xlim()
                else:
                    axes.set_xlim(set_std_xlim[plotcombo])

        print(datasubset)
        print(xydata)
        plt.tight_layout()
        fullname = fig_dir + datasubset + "_" + plotcombo + ftype
        print(fullname)
        # if datasubset == "pga_cav_small":
        #     plt.show()
        #     breakpoint()
        plt.savefig(fullname, dpi=dpi_res)
        # imlims = {"pga": [10**-4, 3*10**5], "cav": [10**-5, 10**3]}
        # plt.axis((10**-4, 1, 10**-4, 10))
        # plt.show()


#############################################################

plt.show()
breakpoint()
print("end")
