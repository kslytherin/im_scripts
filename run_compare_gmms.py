'''
Distance scaling plots GMM comparisons, consider doing magnitude scaling
plotting mean, sigma, phi, and tua for multiple GMMs
Started on Apr 2, 2024
'''

import numpy as np
import matplotlib.pyplot as plt
import get_cgmm_im
import get_gmm_im
import compare_gmms

from gmprocess.utils.constants import UNITS, STATION_METRIC_UNITS
from openquake.hazardlib.imt import CAV, IA

fig_dir = "/Users/kksmith/figures/"
ftype = ".png"
tick_fontsize = 10
ax_fontsize = 12
leg_fontsize = 10
title_fontsize = 14

pre_dist_meas = "Rjb"
dist_meas = pre_dist_meas.lower()
#dist_meas_label = pre_dist_meas
dist_meas_label = "Distance to Earthquake"

myfig_ops = {"makefigs": True, "showfigs": True, "printfigs": False}

myimt = CAV()
# IA comparison
if str(myimt) == "IA":
    main_gmm_set = [
        ("MacedoEtAl2019SInter","BooreEtAl2014"),
        "CampbellBozorgnia2019",
        "SandikkayaAkkar2017Rjb",
    ]
    # con_model_set = ["AtkinsonBoore2003SInter"]
    # main_leg_labels = ["M19-AB03", "CB19", "S17"]
    main_leg_labels = ["M19-B14", "CB19", "SA17"]
    nice_gmm_model_str = main_leg_labels
# CAV comparison
elif str(myimt) == "CAV":
    #main_gmm_set = [("MacedoEtAl2021","CampbellBozorgnia2014"), "CampbellBozorgnia2019", "SandikkayaAkkar2017Rjb"]
    main_gmm_set = [("MacedoEtAl2021","BooreEtAl2014"), "CampbellBozorgnia2019", "SandikkayaAkkar2017Rjb"]
    main_leg_labels = ["M21-B14", "CB19", "SA17"]
    nice_gmm_model_str = ["Macedo et al. (2021) \n[Boore et al. (2014)]","Campbell and Bozorgnia (2019)","Sandikkaya and Akkar (2017)"]

main_col = ["r","g","b"]
main_ls = ["-","--","-."]
n = 100
eqmagset = np.array([4, 5, 6, 7, 8])
#eqmagset = np.array([4])
dist_min = 0.1 
dist_max = 300 
dist_pts = np.logspace(np.log10(dist_min), np.log10(dist_max), n)

for e_i, mag_num in enumerate(eqmagset):
    # Specifying earthquake parameters
    myeqparams = {
        "mag": mag_num,
        "rrup": dist_pts,
        "rhypo": dist_pts,
        "hypo_depth": 15,
        "vs30": 424,
        "ztor": 15,
        "backarc": False,
        "rake": 0,
        "rx": dist_pts,
        "width": 100,
        "dip": 90,
        "ry0": dist_pts,
        "rjb": dist_pts,
        "z2pt5": 0.6,
        "z1pt0": 0.1,
        "vs30measured": True,
    }
    #breakpoint()
    mean_set, lnmean_set, sigma_set, phi_set, tau_set, ctx, myplots = compare_gmms.compare_gmms(myimt,main_gmm_set,myeqparams,dist_meas,fig_ops={"makefigs": False, "showfigs": False, "printfigs": False})
    #myplots["mean"]

    # Plotting data
    stat_sets = [mean_set, sigma_set, phi_set, tau_set]
    # should all have the same units?
    simple_stat_type_str = ["mean", "sigma", "phi", "tau"]
    stat_type_str = ["mean", "$\sigma$", "$\phi$", r"$\tau$"]
    myplots = {}
    if myfig_ops["makefigs"]:
#,layout="compressed"
        fig, ax = plt.subplots(2,2,sharex='col',squeeze=True,layout="constrained",figsize=(7.2,7.2))
        #breakpoint()
        for st_i, stat_type in enumerate(stat_sets):
            ax1 = int(np.floor(st_i / 2))
            ax2 = int(st_i % 2)
            # plot mean and uncertainties
            #fig = plt.figure()
            for m_i, main_model in enumerate(main_gmm_set):
                if simple_stat_type_str[st_i] == "mean":
                    ax[ax1,ax2].loglog(
                        ctx[dist_meas],
                        stat_sets[st_i][m_i],
                        linestyle=main_ls[m_i],
                        #color=main_col[m_i],
                        color='g',
                        linewidth=2,
                        label=nice_gmm_model_str[m_i],
                    )
                else:
                    ax[ax1,ax2].semilogx(
                        ctx[dist_meas],
                        stat_sets[st_i][m_i],
                        linestyle=main_ls[m_i],
                        #color=main_col[m_i],
                        color='g',
                        linewidth=2,
                        label=nice_gmm_model_str[m_i]                    )

            if stat_type_str[st_i] == "mean":
                if str(myimt) == "IA":
                    unit_label = UNITS["ariasintensity"]
                else:
                    unit_label = UNITS[str(myimt).lower()]
                add_unit_label = ", " + unit_label
            else:
                if str(myimt) == "IA":
                    unit_label = "ln(" + UNITS["ariasintensity"] + ")"
                else:
                    unit_label = "ln(" + UNITS[str(myimt).lower()] + ")"
                # override, since std dev in log space is more of a pct difference
                unit_label = ""
                add_unit_label = unit_label
            ax[ax1,ax2].set_xlim(dist_min,dist_max)
            ax[ax1,ax2].set_ylabel(
                str(myimt) + " " + stat_type_str[st_i] + add_unit_label,
                fontsize=ax_fontsize,
            )
            if ax1 == 1: 
                ax[ax1,ax2].set_xlabel(
                    dist_meas_label + ", " + STATION_METRIC_UNITS[dist_meas + "_mean"],
                    fontsize=ax_fontsize,
                )
            ax[ax1,ax2].yaxis.set_ticks_position("both")
            ax[ax1,ax2].xaxis.set_ticks_position("both")


            if simple_stat_type_str[st_i] != "mean":
                #plt.axhline(0, color="lightgrey", linestyle="--")
                ax[ax1,ax2].hlines(0, dist_min, dist_max, color="lightgrey", linestyle="--")
            #plt.axhline(0, color="lightgrey", linestyle="--")

            ax[ax1,ax2].set_box_aspect(1)
            if st_i == 1:
                ax[ax1,ax2].legend(loc="lower left",fontsize=leg_fontsize,fancybox=True)
            ax[ax1,ax2].tick_params(axis="both", which="major", labelsize=tick_fontsize)
        fig.suptitle(
            "M=" + str(myeqparams["mag"]) + ", Vs30=" + str(myeqparams["vs30"]),
            fontsize=title_fontsize
        )
            # subfname = (
            #     str(myimt) + "_" + simple_stat_type_str[st_i] + "_EQ_" + str(e_i) + ftype
            # )
        fname = (
            str(myimt) + "_EQ_" + str(e_i) + ftype
        )
        fullfname = fig_dir + fname
        print(fullfname)
        if myfig_ops["printfigs"]:
            plt.savefig(fullfname, dpi=300)

            #myplots[simple_stat_type_str[st_i]] = ax 


if myfig_ops["showfigs"]:
    plt.show()

breakpoint()
print("end")
