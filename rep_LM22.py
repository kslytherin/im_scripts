import sys

import scipy.constants
import get_gmm_im
import get_cgmm_im
import numpy as np
from openquake.hazardlib.imt import CAV, PGA, SA
from kbcg_cav import *
from bchdro_cav_new import *

# from gmprocess.utils.constants import UNITS
import matplotlib.pyplot as plt
from mytools import *
import scipy
import pandas as pd

# try reproducing Fig 8, 10
# maybe try Fig 14, 16

# Have an issue with replicating paper example since it does not specify:
# cav1100

final_imt = CAV()

# tgmm_set = ["LiuMacedo2022M2SSlab", "LiuMacedo2022M2SInter", "LiuMacedo2022M2SSlab", "LiuMacedo2022M2SInter"]

n = 100
gfact = scipy.constants.g
ofact = 100
# repfig = "likeF10b"
repfig = "likeF8a"

minmax_dist = [10, 200]
dist_name = "Rrup"
dist_meas = dist_name.lower()

myX = dist_meas
# myX = "vs30"

dist_pts = np.logspace(np.log10(minmax_dist[0]), np.log10(minmax_dist[1]), n)

myylim = [10**-3, 10**3]  # for replicating
# myylim = [0.1, 30]  # for diagnosis

fig_ops = {"makefigs": False, "showfigs": False, "printfigs": False}
#    vs30set = np.array([470, 400, 470, 400])
#    ztorset = np.array([60, 30, 60, 30])
#    eqmagset = np.array([6, 7, 7.5, 8.5])

if repfig == "likeF8a":
    if myX == "rrup":
        dist_pts = np.logspace(np.log10(minmax_dist[0]), np.log10(minmax_dist[1]), n)
        vs30set = 400
        eqmagset = 7
    elif myX == "vs30":
        dist_pts = 100
        vs30set = np.linspace(200, 1000, n, endpoint=True)
        eqmagset = 7
    elif myX == "mag":
        dist_pts = 100
        vs30set = 500
        eqmagset = np.linspace(4, 8, 20, endpoint=True)
    ztorset = 30
    # tgmm_set = ["LiuMacedo2022M1SInter"]
    # tgmm_set = ["LiuMacedo2022M1SSlab"]
    # tgmm_set = ["LiuMacedo2022M2SInter"]
    tgmm_set = ["LiuMacedo2022M1SInter", "LiuMacedo2022M2SInter"]
    main_cgmm_name_set = ["LiuMacedo2022SInter"]
    con_model_set = ["ParkerEtAl2020SInter", "KuehnEtAl2020SInter"]
elif repfig == "likeF8b":
    vs30set = 400
    ztorset = 30
    eqmagset = 8.5
    # tgmm_set = ["LiuMacedo2022M1SInter"]
    tgmm_set = ["LiuMacedo2022M1SInter", "LiuMacedo2022M2SInter"]
    main_cgmm_name_set = ["LiuMacedo2022SInter"]
    con_model_set = ["ParkerEtAl2020SInter"]
elif repfig == "likeF10a":
    vs30set = 470
    ztorset = 60
    eqmagset = 6
    tgmm_set = ["LiuMacedo2022M1SSlab", "LiuMacedo2022M2SSlab"]
    main_cgmm_name_set = ["LiuMacedo2022SSlab"]
    con_model_set = ["ParkerEtAl2020SSlab"]
elif repfig == "likeF10b":
    vs30set = 470
    ztorset = 60
    eqmagset = 7.5
    tgmm_set = ["LiuMacedo2022M1SSlab", "LiuMacedo2022M2SSlab"]
    main_cgmm_name_set = ["LiuMacedo2022SSlab"]
    con_model_set = ["ParkerEtAl2020SSlab", "KuehnEtAl2020SSlab"]
else:
    print("not in set!")
    breakpoint()
# elif repfig = "likeF14":
#
# elif repfig = "likeF16":

# pre_computed_dta from matlab scripts, the issue must be how the units are outputted in the conditioned, model
gmm_pre_computed_dta_cols = ["g", "r"]
cgmm_cols = ["g", "r"]
# mod_colors = ["b", "c", "orange", "darkorange"]
mod_colors = ["b", "purple", "orange", "darkorange"]
mod_symbols = ["--", "--", ":", ":"]

tmod_cols = {}
mod_sym = {}
Omedian_col = {}
for tmod_i, tmodel in enumerate(tgmm_set):
    if tmodel in ["LiuMacedo2022M1SInter", "LiuMacedo2022M1SSlab"]:
        tmod_cols[tmodel] = "b"
        Omedian_col[tmodel] = "c"
    elif tmodel in ["LiuMacedo2022M2SInter", "LiuMacedo2022M2SSlab"]:
        tmod_cols[tmodel] = "purple"
        Omedian_col[tmodel] = "r"
    else:
        tmod_cols[tmodel] = mod_colors[tmod_i]
    mod_sym[tmodel] = mod_symbols[tmod_i]

# convert from g-s to m/s
# conv_fact = -np.log(1 / gfact)
conv_fact = 0
# convert from m/s to cm/s
conv_fact2 = 0

imlist = [PGA(), SA(1.0)]

median_pre_computed_dta_set = {}
median_cgmm_set = {}

sigma_set = {}
phi_set = {}
tau_set = {}
median_set = {}
ctx_set = {}
Omedian = {}

# region = "1_Alaska"
# region = "2_Cascadia"
# region = "3_CentralAmerica&Mexico"
# region = "4_Japan"
# region = "5_NewZealand"
# region = "6_SouthAmerica"
# region = "7_Taiwan"
region = "0_global"

for m_i, tgmm in enumerate(tgmm_set):
    # a1 = 0
    # a2 = 100
    # Rmax = a1+a2*tgmm
    # dist_filter = dist_pts <= Rmax
    # dist_pts = dist_pts[dist_filter]
    auth = split_cap_words_crop(con_model_set[m_i])[0]
    tect = split_cap_words_crop(con_model_set[m_i])[-1]

    # "cav1100": tgmm/dist_pts,
    # From MATLAB script
    datafile = "/Users/kksmith/Downloads/parker_and_nico/"
    full_dtafile_name = datafile + "pre_im_values_" + auth + "_" + tect + ".csv"
    data = pd.read_csv(full_dtafile_name)
    print("reading " + full_dtafile_name)
    myeqparams = {
        # "mag": eqmagset,
        # "vs30": vs30set,  # from text
        # "ztor": ztorset,
        # "hypo_depth": ztorset,
        # dist_meas: dist_pts,
        "mag": data["M"].values,
        "vs30": data["Vs30"].values,  # from text
        "ztor": data["hyperD"].values,
        "hypo_depth": data["hyperD"].values,
        "PGA_MEAN": data["PGAmean"].values - 1,
        "SA(1.0)_MEAN": data["SA1p0mean"].values - 1,
        dist_meas: data["Rrup"].values,
        "PGA_TOTAL_STDDEV": data["PGAsigma"].values,
        "SA(1.0)_TOTAL_STDDEV": data["SA1p0sigma"].values,
    }

    if 1 == 0:
        lnmedian_pre_C = {}
        sigma_C = {}
        tau_C = {}
        phi_C = {}
        ctx_C = {}
        imylim = {0: [0, 0.75], 1: [0, 0.3]}
        for im_i, im in enumerate(imlist):
            imstr = str(im)
            for cmodel in con_model_set:
                (
                    lnmedian_pre_C[imstr],
                    sigma_C[imstr],
                    tau_C[imstr],
                    phi_C[imstr],
                    ctx_C[imstr],
                ) = get_gmm_im.get_gmm_im(myeqparams, cmodel, im)
                fig1 = plt.figure()
                ax1 = fig1.add_subplot()
                ax1.semilogx(myeqparams["rrup"], np.exp(lnmedian_pre_C[imstr]))
                ax1.set_ylabel(imstr)
                ax1.set_xlabel("Rrup, km")
                ax1.set_xlim([10, 200])
                ax1.set_ylim(imylim[im_i])
                ax1.set_title(
                    cmodel
                    + ", M{}".format(eqmagset)
                    + ", Vs30={} ".format(myeqparams["vs30"])
                    + "m/s"
                )
        plt.show()

    # Using precomputed data
    (
        lnmedian_pre_computed_dta,
        sigma_pre_computed_dta,
        tau_pre_computed_dta,
        phi_pre_computed_dta,
        ctx_pre_computed_dta,
    ) = get_gmm_im.get_gmm_im(
        myeqparams,
        main_cgmm_name_set[0],
        final_imt,
        fig_ops,
    )

    lnmedian_pre_computed_dta = lnmedian_pre_computed_dta + conv_fact + conv_fact2
    median_pre_computed_dta = np.exp(lnmedian_pre_computed_dta)
    median_pre_computed_dta_set[con_model_set[m_i]] = median_pre_computed_dta

    # Using TGMMs
    lnmedian_pre, sigma, tau, phi, ctx = get_gmm_im.get_gmm_im(
        myeqparams, tgmm, final_imt
    )
    lnmedian = lnmedian_pre + conv_fact + conv_fact2
    median = np.exp(lnmedian)
    median_set[tgmm] = median
    # median_alt = np.exp(lnmedian_pre) * gfact * ofact
    # median_alt = np.exp(lnmedian_pre)

    # Using CGMMs like normal
    # only works because (len(con_model_set) == len(tgmm_set)) and (len(main_cgmm_name_set) == 1)
    lnmedian_pre_cgmm, sigma_cgmm, tau_cgmm, phi_cgmm, ctx_cgmm = (
        get_cgmm_im.get_cgmm_im(
            myeqparams, main_cgmm_name_set[0], con_model_set[m_i], final_imt
        )
    )
    lnmedian_cgmm = lnmedian_pre_cgmm + conv_fact + conv_fact2
    median_cgmm = np.exp(lnmedian_cgmm)
    # median_cgmm_set[con_model_set[m_i]] = median_cgmm * gfact
    median_cgmm_set[con_model_set[m_i]] = median_cgmm

    # breakpoint()

    # print("Final")
    # print(con_model_set[m_i])
    # print("tau")
    # print(tau_cgmm[0])
    # print("phi")
    # print(phi_cgmm[0])
    # print("sigma")
    # print(sigma_cgmm[0])
    if split_cap_words_crop(tgmm)[-1] == "Inter":
        mechanism = "interface"
    elif split_cap_words_crop(tgmm)[-1] == "Slab":
        mechanism = "intraslab"
    # else:
    #     breakpoint()

    if split_cap_words(tgmm)[3] == "M1":
        Omedian[tgmm] = np.array(
            kbcg_cav(
                myeqparams["mag"],
                myeqparams["rrup"],
                myeqparams["vs30"],
                myeqparams["ztor"],
                mechanism,
                region,
            )
        )
    elif split_cap_words(tgmm)[3] == "M2":
        Omedian[tgmm] = np.array(
            bchdro_cav_new(
                myeqparams["mag"],
                myeqparams["rrup"],
                myeqparams["vs30"],
                myeqparams["ztor"],
                mechanism,
                region,
            )
        )
    ctx_set[tgmm] = ctx
    sigma_set[tgmm] = sigma
    phi_set[tgmm] = phi
    tau_set[tgmm] = tau

    # plotting of with pre-specified pga and SA data
    # lnmedian = lnmedian_pre_computed_dta + conv_fact + conv_fact2

    print(m_i)

# breakpoint()

stats = {"median": median_set, "phi": phi_set, "tau": tau_set, "sigma": sigma_set}
yes = get_dict_key_len(sigma_set)
yes2 = get_dict_key_len(phi_set)
yes3 = get_dict_key_len(stats)
stat_names = ["median"]
my_y_lims = {"median": myylim, "phi": [0.4, 0.8]}

stat_str = {"median": str(final_imt) + ", m/s", "phi": "$\phi$"}
# stat_str = {"median": str(final_imt) + ", g-s", "phi": "$\phi$"}
for st_i, SN in enumerate(stat_names):
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    for m_i, tgmm in enumerate(tgmm_set):
        new_x_lims = minmax_dist
        kwargs_tgmm = {
            "color": tmod_cols[tgmm],
            # "linestyle": mod_sym[tgmm],
            "linestyle": "-",
            "label": tgmm,
        }
        kwargs_cgmm = {
            "color": cgmm_cols[m_i],
            # "linestyle": mod_sym[tgmm],
            "linestyle": "--",
            "label": main_cgmm_name_set[0] + "-" + con_model_set[m_i],
        }
        kwargs_precdtacgmm = {
            "color": cgmm_cols[m_i],
            "linestyle": ":",
            "label": "pre-computed data of "
            + main_cgmm_name_set[0]
            + "-"
            + con_model_set[m_i],
        }
        kwargs_2 = {
            "color": Omedian_col[tgmm],
            "linestyle": mod_sym[tgmm],
            # "label": tgmm,
        }
        if SN == "median":
            # plotting TGMM
            ax.semilogy(ctx_set[tgmm][myX], stats[SN][tgmm], **kwargs_tgmm)

            # plotting CGMMs as is
            ax.semilogy(
                ctx_set[tgmm][myX], median_cgmm_set[con_model_set[m_i]], **kwargs_cgmm
            )

            # plotting CGMM w/pre-computed data
            ax.semilogy(
                ctx_set[tgmm][myX],
                median_pre_computed_dta_set[con_model_set[m_i]],
                gmm_pre_computed_dta_cols[m_i],
                **kwargs_precdtacgmm
            )

            # ax.semilogy(ctx_set[tgmm][myX], Omedian[tgmm], **kwargs_2)
        elif SN == "phi":
            ax.plot(ctx_set[tgmm][myX], stats[SN][tgmm], **kwargs_tgmm)
    if myX == "rrup":
        ax.set_xlabel("Rrup, km")
        ax.set_xscale("log")
        ax.set_title(
            "M{}".format(eqmagset) + ", Vs30={} ".format(myeqparams["vs30"][0]) + "m/s"
        )
    elif myX == "vs30":
        ax.set_xlabel("Vs30, m/s")
        ax.set_title(
            "M{}".format(eqmagset) + ", Rrup={} ".format(myeqparams["rrup"]) + "km"
        )
    else:
        ax.set_xlabel(myX)

    ax.set_ylabel(stat_str[SN])
    ax.set(
        xlim=(new_x_lims[0], new_x_lims[1]), ylim=(my_y_lims[SN][0], my_y_lims[SN][1])
    )
    # ax.set_xlim(new_x_lims[0], new_x_lims[1])
    ax.yaxis.set_ticks_position("both")
    ax.xaxis.set_ticks_position("both")
    plt.tight_layout()
    plt.legend()
plt.show()
# breakpoint()
print("end")
