"""
plotting a simple GMM
started on Mar 12, 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from openquake.hazardlib import valid
from openquake.hazardlib.contexts import simple_cmaker

# from gmprocess.utils.constants import UNITS, STATION_METRIC_UNITS

eq_param_names = [
    ("backarc", bool),
    ("dip", float),
    ("hypo_depth", float),
    ("mag", float),
    ("rake", float),
    ("repi", float),
    ("rhypo", float),
    ("rjb", float),
    ("rrup", float),
    ("rx", float),
    ("ry0", float),
    ("vs30", float),
    ("vs30measured", bool),
    ("width", float),
    ("ztor", float),
    ("z1pt0", float),
    ("z2pt5", float),
]


def get_gmm_im(
    EQ_params,
    main_gmm,
    imt,
    fig_ops={"makefigs": False, "showfigs": False, "printfigs": False},
    reg=None,
):

    imt_str = str(imt)
    if reg == None:
        gsim = valid.gsim("[" + main_gmm + "]")
    else:
        gsim = valid.gsim("[" + main_gmm + "]\nregion='" + reg + "'")

    # Build ctx file automatically given the EQ_parm.keys() from all available types given in eq_param_names
    pre_ctx_ids = []
    pci_l = []
    for w_i in range(np.shape(eq_param_names)[0]):
        if eq_param_names[w_i][0] in EQ_params.keys():
            pre_ctx_ids.append((eq_param_names[w_i][0], eq_param_names[w_i][1]))
            pci_l.append(eq_param_names[w_i][0])

    left_over = list(EQ_params.keys() - pci_l)
    for lo_i in left_over:
        if lo_i in ["eqid", "stnid"]:
            pre_ctx_ids.append((lo_i, str))
        else:
            pre_ctx_ids.append((lo_i, float))

    # pre_ctx_ids = []
    # for p_i, pname in enumerate(EQ_params.keys()):
    #     for w_i in range(np.shape(eq_param_names)[0]):
    #         if pname == eq_param_names[w_i][0]:
    #             pre_ctx_ids.append((pname, eq_param_names[w_i][1]))

    n = 0
    for mykey in EQ_params.keys():
        try:
            if n < len(EQ_params[mykey]):
                n = len(EQ_params[mykey])
        except:
            continue

    ctx = np.recarray(
        n,
        dtype=np.dtype(pre_ctx_ids),
    )

    for pname in EQ_params.keys():
        try:
            ctx[pname] = EQ_params[pname]
        except ValueError as inst:
            print(inst)
    mags = np.ones((n,)) * EQ_params["mag"]
    cmaker = simple_cmaker([gsim], [str(imt)], mags=["%2.f" % m for m in mags])

    # Evaluate the GMM.
    result = cmaker.get_mean_stds([ctx])
    # an array of shape (4, G, M, N) with mean and stddevs
    # M = len(self.imts)
    # G = len(self.gsims)
    # N = no. of pts
    # 4 = [mean,phi,tau,std]
    lnmean = result[0][0][0]
    lnstd = result[3][0][0]
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    tau = result[2][0][0]
    # phi = result[1][0][0]
    # sigma = lnstd

    # things are not matching up sigma and phi are switched
    sigma = result[1][0][0]
    phi = lnstd
    mean = np.exp(lnmean)
    if fig_ops["makefigs"]:
        # dist_meas = "repi"
        # dist_meas = "rjb"
        dist_meas = "rrup"
        ### Investigative plots
        # IM vs Distance
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.loglog(ctx[dist_meas], mean, "k")
        ax.loglog(ctx[dist_meas], np.exp(lnmean + lnstd), "0.8", linestyle="--")
        ax.loglog(ctx[dist_meas], np.exp(lnmean - lnstd), "0.8", linestyle="--")
        # try:
        #     ax.set_xlabel(
        #         "$R_{" + dist_meas + "}$, " + STATION_METRIC_UNITS[dist_meas + "_mean"]
        #     )
        # except:
        #     ax.set_xlabel("$R_{" + dist_meas + "}$")
        ax.set_xlabel(dist_meas)

        # ax.set_ylabel(imt_str + ", " + UNITS[imt_str.lower()])
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_ticks_position("both")
        # ax.set(xlim=(1, 100), ylim=(0.01, 10))
        ax.set_title(
            main_gmm
            + ", M="
            + str(EQ_params["mag"][0])
            + ", D="
            + str(EQ_params["hypo_depth"][0])
            + " km?"
        )
        if fig_ops["showfigs"]:
            plt.show()
    return lnmean, sigma, tau, phi, ctx
