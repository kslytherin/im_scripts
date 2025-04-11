"""
Testing OQ input and excel sheet from CB19
Script to read in known EQ input values and 'true' values 
of mean and std devs from CB excel sheets 
starting with CB14_...Apr_22_2023.xlsx 
"""

import get_gmm_im
import get_cgmm_im
import numpy as np
from openquake.hazardlib.imt import PGA, CAV, IA, PGV
import pandas as pd

# read EQ values and CB results
gmm_model = "CampbellBozorgnia2019"

# test_path = "/Users/kksmith/Documents/CB_DATA/V3/"
test_path = "/Users/kksmith/oq-engine/openquake/hazardlib/tests/gsim/data/CB19/"

n = -3
perc_thresholds = {
    "mean": 10**n,
    "sigma": 10**n,
    "phi": 10**n,
    "tau": 10**n,
}

fail_count = 0
total_count = 0
stat_tag_list = ["MEAN", "STD_TOTAL", "STD_INTER", "STD_INTRA"]
for ctr, (stat_tag, stat_type) in enumerate(
    zip(stat_tag_list, ["mean", "sigma", "tau", "phi"])
):

    fname = test_path + "CB19_" + stat_tag + ".csv"
    eqparams = pd.read_csv(fname)
    if ctr == 0:
        myeqparams = {
            "mag": eqparams["rup_mag"],
            "ztor": eqparams["rup_ztor"],
            "hypo_depth": eqparams["rup_hypo_depth"],
            "rx": eqparams["dist_rx"],
            "rjb": eqparams["dist_rjb"],
            "rrup": eqparams["dist_rrup"],
            "width": eqparams["rup_width"],
            "rake": eqparams["rup_rake"],
            "dip": eqparams["rup_dip"],
            "vs30": eqparams["site_vs30"],
            "z2pt5": eqparams["site_z2pt5"],
        }

        my_calc_stats = {}
        for myimt in [CAV(), IA(), PGA()]:
            lnmean_pre, sigma, tau, phi, ctx = get_gmm_im.get_gmm_im(
                myeqparams, gmm_model, myimt
            )
            # convert g-s to m/s
            gfact = 9.81
            conv_fact = 0
            if str(myimt) == "CAV":
                conv_fact = -np.log(1 / gfact)
            lnmean = lnmean_pre + conv_fact
            mean = np.exp(lnmean)
            mean_alt = np.exp(lnmean_pre) * gfact
            my_calc_stats[str(myimt).lower() + "_lnmean"] = lnmean
            # Convert from ln(mean)
            my_calc_stats[str(myimt).lower() + "_mean"] = mean
            my_calc_stats[str(myimt).lower() + "_sigma"] = sigma
            my_calc_stats[str(myimt).lower() + "_tau"] = tau
            my_calc_stats[str(myimt).lower() + "_phi"] = phi

    # breakpoint()

    # Do check
    for myimt in [CAV(), IA(), PGA()]:
        CBval = eqparams[str(myimt).lower()]
        per_diff = np.abs(
            (CBval - my_calc_stats[str(myimt).lower() + "_" + stat_type])
            / CBval * 100
        )
        Tval = per_diff < perc_thresholds[stat_type]

        for tval_i, tval in enumerate(Tval):
            total_count += 1
            if not tval:
                fail_count += 1
                print(
                    str(myimt).upper()
                    + " "
                    + stat_type.upper()
                    + " is FAILURE for EQ on line "
                    + str(tval_i + 2)
                    + ": M="
                    + str(eqparams["rup_mag"][tval_i])
                    + ", Rjb="
                    + str(eqparams["dist_rjb"][tval_i])
                    + " km, Vs30 = "
                    + str(eqparams["site_vs30"][tval_i])
                    + " m/s "
                )
                print(
                    "CBval: "
                    + str(CBval[tval_i])
                    + ", myval: "
                    + str(my_calc_stats[str(myimt).lower() + "_"
                                        + stat_type][tval_i])
                    + ", percent diff is "
                    + "{:.1f}".format(per_diff[tval_i])
                    + "%"
                    + "\n"
                )
                #breakpoint()

    if fail_count != 0:
        print("\n")

if fail_count == 0:
    print("Total SUCCESS!!")
else:

    print(
        "Failure rate = {}/{}={:.1%}".format(
            fail_count, total_count, fail_count / total_count
        )
    )
