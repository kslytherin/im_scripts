"""
Getting estimated im from a gmm but for only one EQ at a time, one model at a time
Started on Mar 12, 2024
"""

import numpy as np
import matplotlib.pyplot as plt

from gmprocess.utils.constants import UNITS, STATION_METRIC_UNITS


# Specifying as many earthquake parameters as possible (Master list)
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


def get_cgmm_im(
    EQ_params,
    main_gmm,
    cgmm,
    imt,
    fig_ops={"makefigs": False, "showfigs": False, "printfigs": False},
    # con_data = None,
):
    # Build ctx file automatically given the EQ_parm.keys() from all available types given in eq_param_names
    pre_ctx_ids = []
    for p_i, pname in enumerate(EQ_params.keys()):
        for w_i in range(np.shape(eq_param_names)[0]):
            if pname == eq_param_names[w_i][0]:
                pre_ctx_ids.append((pname, eq_param_names[w_i][1]))

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

    # Evaluate the GMM.
    lnmean = np.zeros([1, n])
    sigma = np.zeros_like(lnmean)
    tau = np.zeros_like(lnmean)
    phi = np.zeros_like(lnmean)

    if main_gmm == "MacedoEtAl2019SInter":
        from openquake.hazardlib.gsim.macedo_2019 import MacedoEtAl2019SInter

        gsim = MacedoEtAl2019SInter(gmpe={cgmm: {}})
    elif main_gmm == "MacedoEtAl2021":
        from openquake.hazardlib.gsim.macedo_2021 import MacedoEtAl2021

        gsim = MacedoEtAl2021(gmpe={cgmm: {}})
    elif main_gmm == "LiuMacedo2022SSlab":
        from openquake.hazardlib.gsim.liu_macedo_2022 import LiuMacedo2022SSlab

        gsim = LiuMacedo2022SSlab(gmpe={cgmm: {}})
    elif main_gmm == "LiuMacedo2022SInter":
        from openquake.hazardlib.gsim.liu_macedo_2022 import LiuMacedo2022SInter

        gsim = LiuMacedo2022SInter(gmpe={cgmm: {}})
    gsim.compute(ctx, [imt], lnmean, sigma, tau, phi)
    # if str(imt) == 'CAV':
    #    breakpoint()

    # Convert from ln(mean)
    mean = np.exp(lnmean)

    if fig_ops["makefigs"]:
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        dist_meas = "rrup"

        ### Investigative plots
        # IM vs Distance
        ax.loglog(ctx[dist_meas], mean[0], color="k")
        ax.loglog(ctx[dist_meas], np.exp(lnmean[0] + sigma[0]), "0.8", linestyle="--")
        ax.loglog(ctx[dist_meas], np.exp(lnmean[0] - sigma[0]), "0.8", linestyle="--")

        ax.set_xlabel("$R_{rup}$, " + STATION_METRIC_UNITS[dist_meas + "_mean"])
        if str(imt) == "IA":
            ax.set_ylabel(str(imt) + ", " + UNITS["ariasintensity"])
        else:
            ax.set_ylabel(str(imt) + ", " + UNITS[str(imt).lower()])

        ax.set_title(
            main_gmm
            + ", M="
            + str(EQ_params["mag"])
            + ", Vs30="
            + str(EQ_params["vs30"])
        )
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_ticks_position("both")
        ax.legend()
        if fig_ops["showfigs"]:
            plt.show()

    return lnmean[0], sigma[0], tau[0], phi[0], ctx
