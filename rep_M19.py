"""
Replicating Figure 5 in Macedo et al 2019
Started on Mar 12, 2024
"""

import numpy as np
import matplotlib.pyplot as plt
import get_cgmm_im

from gmprocess.utils.constants import UNITS, STATION_METRIC_UNITS
from openquake.hazardlib.imt import IA

# 5 Models used by Macedo for PGA

# Papers of models, OQ code
# Zhao et al (2006), zhao_2006.py
# Montalva et al (2017), montalva_2017.py
# Abrahamson et al (2018), abrahamson_2018.py (BC hydro)
# Atkinson and Boore (2003,2008), atkinson_boore_2003.py
# BChydro 2016, bchydro_2016_epistemic.py (Not sure!!, maybe no)

# the following are the available GMM models in OQ that we want for Macedo etal 2019,
# con_model_set = ["AbrahamsonEtAl2018SInter","ZhaoEtAl2006SInter","AtkinsonBoore2003SInter","MontalvaEtAl2017SInter"]


main_gmm = "MacedoEtAl2019SInter"
myimt = IA()

show_plt = True
num_dist_pts = 100
# eqmagset = np.array([4.75, 5.75, 7, 7.6])
eqmagset = np.array([7, 7.5, 8.5, 9])
eqvs30set = np.array([300, 300, 760, 760])
# eqvs30set = np.array([300, 300, 300, 300])

# Montalva is not shown in Fig 5, so I left it out
con_model_set = [
    "AtkinsonBoore2003SInter",
    "ZhaoEtAl2006SInter",
    "AbrahamsonEtAl2018SInter",
]
con_leg_labels = ["G_AB2003,2008", "G_Z2006", "G_BCH2018"]
con_line_type = ["-.", ":", "-"]
dist_pts = np.logspace(-1, np.log10(300), num_dist_pts)

for mag_num, vs30_num in zip(eqmagset, eqvs30set):
    # Specifying earthquake parameters
    myeqparams = {
        "mag": mag_num,
        "rrup": dist_pts,
        "rhypo": dist_pts,
        "hypo_depth": 20,
        "vs30": vs30_num,
        "backarc": False,
        "rake": 0,
        "rx": dist_pts,
        "width": 100,
        "dip": 90,
        "rjb": dist_pts,
        "z2pt5": 0.2,
        "z1pt0": 0.1,
    }

    dist_meas = "rrup"

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    for m_i, con_model in enumerate(con_model_set):
        lnmean, sigma, tau, phi, ctx = get_cgmm_im.get_cgmm_im(
            myeqparams, main_gmm, con_model, myimt
        )
        mean = np.exp(lnmean)

        # IM vs Distance
        ax.loglog(
            ctx[dist_meas],
            mean,
            linestyle=con_line_type[m_i],
            label=con_leg_labels[m_i],
            color="k",
        )

    ax.set_xlabel("$R_{rup}$, " + STATION_METRIC_UNITS[dist_meas + "_mean"])
    if str(myimt) == "IA":
        ax.set_ylabel("IA, " + UNITS["ariasintensity"])
    else:
        ax.set_ylabel(str(myimt) + UNITS[str(myimt)])

    ax.set_title(
        main_gmm + ", M=" + str(myeqparams["mag"]) + ", Vs30=" + str(myeqparams["vs30"])
    )
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")
    ax.set(xlim=(10, 1000), ylim=(0.001, 100))
    ax.set_box_aspect(1)
    ax.legend()
if show_plt:
    plt.show()
breakpoint()
print("end")