import get_gmm_im
import numpy as np
from openquake.hazardlib.imt import CAV, IA
from gmprocess.utils.constants import UNITS
import matplotlib.pyplot as plt

# try reproducing fig 5 mag 7.5 event

myimt = CAV()
dist_meas_set = [
    "rjb",
    "repi",
    "rhypo",
]
dist_meas_col = ["k", "r", "b"]

conv_fact = 0
conv_fact2 = 0
gfact = 9.81
ofact = 100

if str(myimt) == "CAV":
    myylim = [1, 10000]
    # convert from g-s to m/s
    conv_fact = -np.log(1 / gfact)
    # convert from m/s to cm/s
    conv_fact2 = -np.log(1 / ofact)
elif str(myimt) == "IA":
    myylim = [10**-6, 100]

# Specifying earthquake parameters
n = 100
minmax_dist = [1, 200]
dist_pts = np.logspace(np.log10(minmax_dist[0]), np.log10(minmax_dist[1]), n)

eqmagset = np.array([5, 7.5])
# eqmagset = np.array([4.5, 6, 7.5])
# eqmagset = np.array([7.5])
# eqmagset = np.array([5])

fig = plt.figure()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

main_name = "SandikkayaAkkar2017"
for dm_i, dist_meas in enumerate(dist_meas_set):
    if dist_meas == "rjb":
        dist_name = "Rjb"
    elif dist_meas == "repi":
        dist_name = "Repi"
    elif dist_meas == "rhypo":
        dist_name = "Rhyp"

    gmm_model = main_name + dist_name
    for m_i, mag_num in enumerate(eqmagset):
        myeqparams = {
            "mag": mag_num,
            "vs30": 400,  # 424, 400 is orig
            dist_meas: dist_pts,
            "rake": 0,
        }

        lnmean_pre, sigma, tau, phi, ctx = get_gmm_im.get_gmm_im(
            myeqparams, gmm_model, myimt
        )
        # convert units g-s to cm/s
        lnmean = lnmean_pre + conv_fact + conv_fact2
        # Convert from ln(mean)
        mean = np.exp(lnmean)
        mean_alt = np.exp(lnmean_pre) * gfact * ofact
        (line,) = ax.loglog(ctx[dist_meas], mean, dist_meas_col[dm_i])
        if m_i == 1:
            line.set_label(dist_meas)

ax.yaxis.set_ticks_position("both")
ax.xaxis.set_ticks_position("both")
ax.set_xlabel("R, km")
# breakpoint()
if str(myimt) == "IA":
    ax.set_ylabel(str(myimt) + ", " + UNITS["arias"])
    # ax.set_ylabel(str(myimt) + ", " + UNITS["arias"])
else:
    ax.set_ylabel(str(myimt) + ", cm/s")

ax.set_title(main_name + ", Vs30=" + str(myeqparams["vs30"]) + " m/s")
ax.set(xlim=(minmax_dist[0], minmax_dist[1]), ylim=(myylim[0], myylim[1]))
ax.set_box_aspect(1)
ax.legend()
# plt.savefig("")
# plt.show()
