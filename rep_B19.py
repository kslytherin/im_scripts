import get_gmm_im
import numpy as np
from openquake.hazardlib.imt import CAV
from gmprocess.utils.constants import UNITS
import matplotlib.pyplot as plt
from mytools import get_dict_key_len

# try reproducing Fig 12, Fig 15

# Paper does not specify:
# Tectonic region
# Rtvz


myimt = CAV()

gfact = 9.81
ofact = 100

mag_cols = {4:'k', 6:'k', 8:'k'}
mag_sym = {4:'-.', 6:'--', 8:'-'}
eqmagset = np.array([4, 6, 8])

# mag_cols = {4:'b'}
# mag_sym = {4:'--'}
# eqmagset = np.array([4])

# Specifying earthquake parameters
n = 100

minmax_dist = [1, 200]
myylim = [0.1, 10 ** 4]

# convert from g-s to m/s
# conv_fact = -np.log(1 / gfact)
conv_fact = 0
# convert from m/s to cm/s
conv_fact2 = 0

sigma_set = {}
phi_set = {}
tau_set = {}
mean_set = {}
ctx_set = {}

main_name = "BullockEtAl2019AscRotD100"
dist_meas = "rjb"
dist_name = "Rjb"

gmm_model = main_name

for m_i, mag_num in enumerate(eqmagset):
    dist_pts = np.logspace(np.log10(minmax_dist[0]), np.log10(minmax_dist[1]), n)
    dist_filter = dist_pts <= (77.5 * mag_num - 220)
    dist_pts = dist_pts[dist_filter]

    myeqparams = {
        "mag": mag_num,
        "vs30": 1100,  # from text 
        "ztor": 8,
        "rake": 0,
        "z1pt0": 0, 
        dist_meas: dist_pts
    }
    lnmean_pre, sigma, tau, phi, ctx = get_gmm_im.get_gmm_im(
        myeqparams, gmm_model, myimt
    )
    ctx_set[mag_num] = ctx
    sigma_set[mag_num] = sigma
    phi_set[mag_num] = phi
    tau_set[mag_num] = tau

    # convert units g-s to cm/s
    lnmean = lnmean_pre + conv_fact + conv_fact2
    mean = np.exp(lnmean)
    mean_set[mag_num] = mean
    #mean_alt = np.exp(lnmean_pre) * gfact * ofact
    mean_alt = np.exp(lnmean_pre)
    #breakpoint()
    print(m_i)

stats = {'mean': mean_set, 'phi': phi_set, 'tau': tau_set, 'sigma': sigma_set}
yes = get_dict_key_len(sigma_set)
yes2 = get_dict_key_len(phi_set)
yes3 = get_dict_key_len(stats)

#stat_names = ['mean', 'phi']
#my_y_lims = {'mean': myylim, 'phi': [0, 2]}
stat_names = ['mean']
my_y_lims = {'mean': myylim}

stat_str = {'mean': str(myimt) + ', cm/s', 'phi': '$\phi$'}
for st_i, SN in enumerate(stat_names):
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    #ax = fig.add_subplot()
    for m_i, mag_num in enumerate(eqmagset):

        new_x_lims = minmax_dist
        kwargs= {'color': mag_cols[mag_num], 'linestyle': mag_sym[mag_num], 'label': "M{}".format(mag_num)}

        if SN == 'mean':
            ax.loglog(ctx_set[mag_num][dist_meas], stats[SN][mag_num], **kwargs)
        else:
            ax.plot(ctx_set[mag_num][dist_meas], stats[SN][mag_num], **kwargs)

    ax.set_xlabel("Rrjb, km")
    ax.set_ylabel(stat_str[SN])
    ax.set(xlim=(new_x_lims[0], new_x_lims[1]), ylim=(my_y_lims[SN][0], my_y_lims[SN][1]))
    ax.set_title(main_name + ", Vs30=" + str(myeqparams["vs30"]) + " m/s")
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")
    #ax.set_box_aspect(1)
    plt.legend()
    #plt.tight_layout()
plt.show()
breakpoint()
print("end")
