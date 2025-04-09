"""
Replicating Figure 8 in Macedo et al 2021
Started on Mar 12, 2024
"""

import numpy as np
import matplotlib.pyplot as plt
import get_cgmm_im

from gmprocess.utils.constants import UNITS, STATION_METRIC_UNITS
from openquake.hazardlib.imt import CAV

# 5 Models used by Macedo et al 2021 for PGA
# should match up but need to check to make sure

# Papers of models, OQ code
# Campbell and Bozorgnia (2014), campbell_bozorgnia_2014.py
# Abrahamson et al (2014), abrahamson_2014.py
# Boore et al (2014), boore_2014.py
# Chiou and Youngs (2014), chiou_youngs_2014.py
# Idriss (2014), idriss_2014.py


main_gmm = "MacedoEtAl2021"
myimt = CAV()

show_plt = True
n = 100

eqmagset = np.array([4.75, 5.75, 7, 7.6])
# eqmagset = np.array([3, 4, 5, 6, 7])
# eqmagset = np.array([4.75])

# the following are the available GMM models in OQ that we want for
# Macedo et al (2021)
# con_model_set = [
#     "BooreEtAl2014",
# ]
# con_leg_labels = ["BSSA_2014_nga"]
# con_col = ["b"]
con_model_set = [
    "AbrahamsonEtAl2014",
    "BooreEtAl2014",
    "CampbellBozorgnia2014",
    "ChiouYoungs2014",
    "Idriss2014",
]
con_leg_labels = [
    "ASK14",
    "BSSA_2014_nga",
    "CB_2014_nga",
    "CY_2014_nga",
    "I_2014_gna",
]
con_col = ["r", "g", "b", "m", "y"]
dist_pts = np.logspace(-1, np.log10(300), n)

for mag_i, mag_num in enumerate(eqmagset):
    # Specifying earthquake parameters
    myeqparams = {
        "mag": mag_num,
        "rrup": dist_pts,
        "rhypo": dist_pts,
        "hypo_depth": 15,
        "vs30": 424,  # 424 for fig 8, 760 for fig 9
        "ztor": 17,  # 20,
        "backarc": False,
        "rake": 0,
        "rx": dist_pts,
        "width": 100,
        "dip": 90,
        "ry0": dist_pts,
        "rjb": dist_pts,
        "z2pt5": 0.2,
        "z1pt0": 0.1,
        "vs30measured": True,
    }

    # get data
    lnmean_set = {}
    mean_set = {}
    sigma_set = {}
    tau_set = {}
    phi_set = {}
    for m_i, con_model in enumerate(con_model_set):
        lnmean_pre, sigma, tau, phi, ctx = get_cgmm_im.get_cgmm_im(
            myeqparams, main_gmm, con_model, myimt
        )

        gfact = 9.81
        conv_fact = 0
        # convert g-s to m/s
        if str(myimt) == "CAV":
            conv_fact = -np.log(1 / gfact)

        lnmean = lnmean_pre + conv_fact

        # Convert from ln(mean)
        mean = np.exp(lnmean)

        lnmean_set[con_model] = lnmean
        mean_set[con_model] = mean
        sigma_set[con_model] = sigma
        tau_set[con_model] = tau
        phi_set[con_model] = phi
        # print("sigma={}".format(sigma))

    # if mag_i == 1:
    #    breakpoint()

    # plotting means
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    dist_meas = "rrup"
    for m_i, con_model in enumerate(con_model_set):
        ax.loglog(
            ctx[dist_meas],
            mean_set[con_model],
            linestyle="--",
            color=con_col[m_i],
            label=con_leg_labels[m_i],
        )
    ax.set_xlabel("$R_{rup}$, " + STATION_METRIC_UNITS[dist_meas + "_mean"])
    # ax.set_ylabel(str(myimt) + ", " + UNITS[str(myimt).lower()])
    if str(myimt) == "CAV":
        ax.set_ylabel(str(myimt) + ", m/s")
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")

    ax.set_title(
        main_gmm + ", M=" + str(myeqparams["mag"]) + ", Vs30=" + str(myeqparams["vs30"])
    )
    ax.set(xlim=(1, 100), ylim=(0.01, 100))
    ax.set_box_aspect(1)
    ax.legend()

    # plotting standard devations
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    dist_meas = "rrup"
    for m_i, con_model in enumerate(con_model_set):
        ax.semilogx(
            ctx[dist_meas],
            sigma_set[con_model],
            color=con_col[m_i],
            linestyle="--",
            label=con_leg_labels[m_i],
        )
    ax.set_xlabel("$R_{rup}$, " + STATION_METRIC_UNITS[dist_meas + "_mean"])
    # ax.set_ylabel(str(myimt) + " $ \sigma$, " + UNITS[str(myimt).lower()])
    ax.set_ylabel(str(myimt) + " $ \sigma$")
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")

    ax.set_title(
        main_gmm + ", M=" + str(myeqparams["mag"]) + ", Vs30=" + str(myeqparams["vs30"])
    )
    ax.set(xlim=(1, 1000), ylim=(0.35, 0.85))
    ax.set_box_aspect(1)
    ax.legend()
if show_plt:
    plt.show()

breakpoint()
print("end")
