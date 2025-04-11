import get_gmm_im
import numpy as np
from openquake.hazardlib.imt import PGA, CAV
from gmprocess.utils.constants import UNITS, STATION_METRIC_UNITS
import matplotlib.pyplot as plt

# User input
# gmm_model = "CampbellBozorgnia2010"
gmm_model = "CampbellBozorgnia2008"
# copy Fig 4 of CB10
myimt = CAV()

# Specifying earthquake parameters
n = 300
minmax_dist = [1, 100]
dist_pts = np.logspace(np.log10(minmax_dist[0]), np.log10(minmax_dist[1]), n)
eqmagset = np.array([5, 6, 7, 8])


fig = plt.figure()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
dist_meas = "rrup"

for mag_i, mag_num in enumerate(eqmagset):
    myeqparams = {
        "mag": mag_num,
        "rrup": dist_pts,
        "vs30": 760,
        "ztor": 0,
        "rake": 0,
        "hypo_depth": 5,
        "rx": -1,
        "width": 10 ** (-0.76 + 0.27 * mag_num),
        "dip": 90,
        "rjb": dist_pts,
        "z2pt5": 2,
    }

    # if mag_i == 3:
    #    print(myeqparams["width"])
    # breakpoint()
    lnmean, sigma, tau, phi, ctx = get_gmm_im.get_gmm_im(myeqparams, gmm_model, myimt)

    # Convert from ln(g)
    mean = np.exp(lnmean)
    # breakpoint()
    ax.loglog(ctx[dist_meas], mean)
    ax.text(ctx[dist_meas][0], mean[0], "M=" + str(myeqparams["mag"]))

ax.yaxis.set_ticks_position("both")
ax.xaxis.set_ticks_position("both")
ax.set_xlabel("$R_{rup}$, " + STATION_METRIC_UNITS[dist_meas + "_mean"])
ax.set_ylabel(str(myimt) + ", " + UNITS[str(myimt).lower()])

ax.set_title(gmm_model + ", Vs30=" + str(myeqparams["vs30"]))
ax.set(xlim=(minmax_dist[0], minmax_dist[1]), ylim=(0.01, 10))
ax.set_box_aspect(1)
ax.legend()
plt.show()
