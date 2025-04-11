import get_gmm_im
import numpy as np
from openquake.hazardlib.imt import PGA, CAV, IA, PGV

# from gmprocess.utils.constants import UNITS, STATION_METRIC_UNITS
import matplotlib.pyplot as plt

# User input
# gmm_model = "CampbellBozorgnia2014"

gmm_model = "CampbellBozorgnia2019"
# replicate Fig 3(c) of CB19
myimt = CAV()
# Specifying earthquake parameters
# n = 4
n = 100
# Figure 3c has these limits
minmax_dist = [1, 1000]

# minmax_dist = [5, 1000]

myylim = [0.001, 100]
make_plot = False

# specified in CB14 table C1 of supp and is for vert SS EQ
eqmagset = np.array([3.5, 4.5, 5.5, 6.5, 7.5, 8.5])
zhypset = np.array([7.5626, 8.2623, 7.2779, 8.8705, 10.2267, 10.2267])
ztorset = np.array([7.1449, 7.1449, 4.2887, 0.8741, 0, 0])
zbotset = np.ones_like(ztorset) * 15
Wset = np.array([0.5119, 1.6572, 5.3653, 14.1259, 15, 15])

temp_vs30 = 300
# Figure 3c has these specs
if str(myimt) == "CAV":
    temp_vs30 = 760
temp_z2p5 = 0.6068

# # specified in CB14 table C2 of supp of reverse SS EQ
# eqmagset = np.array([3.5, 4.5, 5.5, 6.5, 7.5, 8.5])
# zhypset = np.array([7.6374, 8.3663, 10.3008, 11.6288, 10.6889])
# ztorset = np.array([7.3116, 7.3116, 7.3116, 3.6324, 0.4622])
# zbotset = np.ones_like(ztorset) * 15
# Wset = np.array([0.5119, 1.6572, 5.3653, 16.0763, 20.5595])

# ztorset = np.array([0])

# Wset = np.array([15])

# tests from CB sheet
# eqmagset = np.array([7.1])
# zhypset = np.array([10.2909])
# ztorset = np.array([0.0642])
# Wset = np.array([18.2168])
# temp_vs30 = 1563
# temp_z2p5 = 0.2660
if str(myimt) == "CAV":
    myylim = [0.001, 100]
elif str(myimt) == "IA":
    myylim = [10**-8, 100]

dist_pts = np.logspace(np.log10(minmax_dist[0]), np.log10(minmax_dist[1]), n)

# "rx": np.flip(np.arange(-299, 1, 1)),
# "rjb": np.arange(0, 300, 1),
# "rrup": np.sqrt(np.arange(0, 300, 1) ** 2 + ztor_num**2),

fig = plt.figure()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
dist_meas = "rx"
for mag_num, ztor_num, zhyp_num, Wset_num in zip(eqmagset, ztorset, zhypset, Wset):
    myeqparams = {
        "mag": mag_num,
        "ztor": ztor_num,
        "hypo_depth": zhyp_num,
        "rx": dist_pts,
        "rjb": dist_pts,
        "rrup": dist_pts,
        "width": Wset_num,
        "rake": 0,
        "dip": 90,
        "vs30": temp_vs30,
        "z2pt5": temp_z2p5,  # is km/s?
    }

    lnmean_pre, sigma, tau, phi, ctx = get_gmm_im.get_gmm_im(
        myeqparams, gmm_model, myimt
    )
    gfact = 9.81
    conv_fact = 0
    # convert g-s to m/s
    # if str(myimt) == "CAV":
    #     conv_fact = -np.log(1 / gfact)

    lnmean = lnmean_pre + conv_fact
    # onvert from ln(mean)
    mean = np.exp(lnmean)
    ax.loglog(ctx[dist_meas], mean)
    ax.text(ctx[dist_meas][0], mean[0], "M=" + str(myeqparams["mag"]))

ax.yaxis.set_ticks_position("both")
ax.xaxis.set_ticks_position("both")
ax.set_xlabel("$R_{x}$, km")

if str(myimt) == "IA":
    # ax.set_ylabel(str(myimt) + ", " + UNITS["arias"])
    ax.set_ylabel(str(myimt) + ", arias units")
else:
    ax.set_ylabel(str(myimt) + ", m/s")

ax.set_title(gmm_model + ", Vs30=" + str(myeqparams["vs30"]))
ax.set(xlim=(minmax_dist[0], minmax_dist[1]), ylim=myylim)
ax.set_box_aspect(1)
ax.legend()
plt.show()
breakpoint()
print("End")
