"""
debugging unit test for macedo
"""

import get_gmm_im
import get_cgmm_im
import numpy as np
from openquake.hazardlib.imt import PGA, CAV
import pandas as pd

# read EQ values and CB results
main_gmm = "MacedoEtAl2021"
con_model = "AbrahamsonEtAl2014"

test_path = "/Users/kksmith/oq-engine/openquake/hazardlib/tests/gsim/data/macedo_2021/"

n = -2
perc_thresholds = {
    "mean": 10**n,
    "sigma": 10**n,
    "phi": 10**n,
    "tau": 10**n,
}

fail_count = 0
total_count = 0
#stat_tag_list = ["MEAN", "STD_TOTAL", "STD_INTER", "STD_INTRA"]
stat_tag_list = ["mean", "total_stddev", "inter_event_stddev", "intra_event_stddev"]

for ctr, (stat_tag, stat_type) in enumerate(
    zip(stat_tag_list, ["mean", "sigma", "tau", "phi"])
):
    fname = test_path + "macedo_2021_" + stat_tag + ".csv"
    eqparams = pd.read_csv(fname)
    if ctr == 0:
        myeqparams = {
            "mag": eqparams["rup_mag"],
            "rrup": eqparams["dist_rrup"],
            "rhypo": eqparams["dist_rhypo"],
            "hypo_depth": eqparams["rup_hypo_depth"],
            "vs30": eqparams["site_vs30"],
            "ztor": eqparams["rup_ztor"],
            "rake": eqparams["rup_rake"],
            "rx": eqparams["dist_rx"],
            "width": eqparams["rup_width"],
            "dip": eqparams["rup_dip"],
            "ry0": eqparams["dist_ry0"],
            "rjb": eqparams["dist_rjb"],
            "z2pt5": eqparams["site_z2pt5"],
            "z1pt0": eqparams["site_z1pt0"],
            "vs30measured": eqparams["site_vs30measured"],
            "backarc": eqparams["site_backarc"]
        }

        my_calc_stats = {}
        for myimt in [CAV(), PGA()]:
            if str(myimt) == "PGA":
                lnmean_pre, sigma, tau, phi, ctx = get_gmm_im.get_gmm_im(
                    myeqparams, con_model, myimt
                )
            else:
                lnmean_pre, sigma, tau, phi, ctx = get_cgmm_im.get_cgmm_im(
                    myeqparams, main_gmm, con_model, myimt
                )
            lnmean = lnmean_pre
            mean = np.exp(lnmean)
            my_calc_stats[str(myimt).lower() + "_lnmean"] = lnmean
            my_calc_stats[str(myimt).lower() + "_mean"] = mean
            my_calc_stats[str(myimt).lower() + "_sigma"] = sigma
            my_calc_stats[str(myimt).lower() + "_tau"] = tau
            my_calc_stats[str(myimt).lower() + "_phi"] = phi

    # Do check
    #for myimt in [CAV(), PGA()]:
    for myimt in [CAV()]:
        true_val = eqparams[str(myimt).lower()]
        per_diff = np.abs(
            (true_val - my_calc_stats[str(myimt).lower() + "_" + stat_type])
            / true_val * 100
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
                    + " true value is "
                    + str(true_val[tval_i])
                    + " calculated value is "
                    + str(my_calc_stats[str(myimt).lower() + "_" + stat_type][tval_i])
                    + ": M="
                    + str(eqparams["rup_mag"][tval_i])
                    + ", Rjb="
                    + str(eqparams["dist_rjb"][tval_i])
                    + " km, Vs30 = "
                    + str(eqparams["site_vs30"][tval_i])
                    + " m/s "
                )
                breakpoint()
                print(
                    "True val: "
                    + str(true_val[tval_i])
                    + ", my calculated value: "
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