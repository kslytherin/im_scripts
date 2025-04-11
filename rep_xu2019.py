import get_gmm_im
import numpy as np
from openquake.hazardlib.imt import CAV
from gmprocess.utils.constants import UNITS
import matplotlib.pyplot as plt

# copy fig 13
shallow = False
if shallow:
    gmm_model = "Xu2019Shallow"
    hdepth = 15
    myylim = [0.01, 10]
else:
    gmm_model = "Xu2019Deep"
    hdepth = 75
    myylim = [0.01, 1]

myimt = CAV()
dist_meas = "repi"

# Specifying earthquake parameters
n = 100
minmax_dist = [1, 200]
dist_pts = np.logspace(np.log10(minmax_dist[0]), np.log10(minmax_dist[1]), n)

eqmagset = np.array([5, 6, 7, 8])
# eqmagset = np.array([5])

# fig_ops = {"makefigs": True, "showfigs": True, "printfigs": False}
fig_ops = {"makefigs": False, "showfigs": False, "printfigs": False}

fig = plt.figure()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

for mag_num in eqmagset:
    myeqparams = {
        "mag": mag_num,
        "vs30": 233,
        "hypo_depth": hdepth,
        "repi": dist_pts,
    }

    # breakpoint()
    lnmean, sigma, tau, phi, ctx = get_gmm_im.get_gmm_im(
        myeqparams, gmm_model, myimt, fig_ops
    )
    # breakpoint()
    # Convert from ln(g)
    mean = np.exp(lnmean)
    ax.loglog(ctx[dist_meas], mean)
    ax.text(ctx[dist_meas][0], mean[0], "M=" + str(myeqparams["mag"]))

ax.yaxis.set_ticks_position("both")
ax.xaxis.set_ticks_position("both")
ax.set_xlabel("$R_{" + dist_meas + "}$, km")
ax.set_ylabel(str(myimt) + ", " + UNITS[str(myimt).lower()])

ax.set_title(
    gmm_model + ", Vs30=" + str(myeqparams["vs30"]) + ", Depth=" + str(hdepth) + " km"
)
ax.set(xlim=(minmax_dist[0], minmax_dist[1]), ylim=(myylim[0], myylim[1]))
ax.set_box_aspect(1)
ax.legend()
plt.show()
