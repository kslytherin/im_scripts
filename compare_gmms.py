'''
Distance scaling plots GMM comparisons, consider doing magnitude scaling
plotting mean, sigma, phi, and tua for multiple GMMs
Started on Apr 2, 2024
'''

import numpy as np
import matplotlib.pyplot as plt
import get_cgmm_im
import get_gmm_im

from gmprocess.utils.constants import UNITS, STATION_METRIC_UNITS
from openquake.hazardlib.imt import CAV, IA

fig_dir = "/Users/kksmith/figures/"
ftype = ".png"
tick_fontsize = 10
ax_fontsize = 12
title_fontsize = 14

def compare_gmms(myimt,main_gmm_set,myeqparams,dist_meas_plot,fig_ops={"makefigs": False, "showfigs": False, "printfigs": False}):
    # getting data!
    mean_set = []
    lnmean_set = []
    phi_set = []
    tau_set = []
    sigma_set = []
    for m_i, main_model in enumerate(main_gmm_set):
        print(main_model)
        if isinstance(main_model,tuple):
            lnmean, sigma, tau, phi, ctx = get_cgmm_im.get_cgmm_im(
                myeqparams, main_model[0], main_gmm_set[1], myimt
            )
        else:
            lnmean, sigma, tau, phi, ctx = get_gmm_im.get_gmm_im(
                myeqparams, main_model, myimt
            )
        #breakpoint()
        # Convert from ln(mean)
        mean = np.exp(lnmean)
        lnmean_set.append(lnmean)
        mean_set.append(mean)
        phi_set.append(phi)
        tau_set.append(tau)
        sigma_set.append(sigma)

    # Plotting data
    stat_sets = [mean_set, sigma_set, phi_set, tau_set]
    # should all have the same units?
    simple_stat_type_str = ["mean", "sigma", "phi", "tau"]
    stat_type_str = ["mean", "$\sigma$", "$\phi$", r"$\tau$"]
    myplots = {}
    if fig_ops["makefigs"]:
        for st_i, stat_type in enumerate(stat_sets):

            # plot mean and uncertainties
            fig, ax = plt.subplots()
            for m_i, main_model in enumerate(main_gmm_set):
                if simple_stat_type_str[st_i] == "mean":
                    ax.loglog(
                        ctx[dist_meas_plot],
                        stat_sets[st_i][m_i],
                        linewidth=2,
                        label=main_gmm_set[m_i],
                    )
                else:
                    ax.semilogx(
                        ctx[dist_meas_plot],
                        stat_sets[st_i][m_i],
                        linewidth=2,
                        label=main_gmm_set[m_i],
                    )

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
                # override, since std dev, in log space, is like a % difference
                unit_label = ""
                add_unit_label = unit_label

            ax.set_ylabel(
                str(myimt) + " " + stat_type_str[st_i] + add_unit_label,
                fontsize=ax_fontsize,
            )
            ax.set_xlabel(
                dist_meas_plot + ", " + STATION_METRIC_UNITS[dist_meas_plot + "_mean"],
                fontsize=ax_fontsize,
            )
            ax.yaxis.set_ticks_position("both")
            ax.xaxis.set_ticks_position("both")

            ax.set_title(
                "M=" + str(myeqparams["mag"]) + ", Vs30=" + str(myeqparams["vs30"]),
                fontsize=title_fontsize,
            )

            if simple_stat_type_str[st_i] != "mean":
                plt.axhline(0, color="lightgrey", linestyle="--")

            ax.set_box_aspect(1)
            ax.legend(loc="lower left")
            ax.tick_params(axis="both", which="major", labelsize=tick_fontsize)
            fname = (
                str(myimt) + "_" + simple_stat_type_str[st_i] + ftype
            )
            fullfname = fig_dir + fname
            if fig_ops["printfigs"]:
                plt.savefig(fullfname, dpi=300)

            myplots[simple_stat_type_str[st_i]] = ax 


    if fig_ops["showfigs"]:
        plt.show()

    return mean_set, lnmean_set, sigma_set, phi_set, tau_set, ctx, myplots,
