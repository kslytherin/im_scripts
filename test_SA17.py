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
import copy

# read EQ values and CB results

test_path = "/Users/kksmith/oq-engine/openquake/hazardlib/tests/gsim/data/SA17/"

n = -3
perc_thresholds = {
    "mean": 10**n,
}

fail_count = 0
total_count = 0
stat_tag = "MEAN"
stat_type = stat_tag.lower()

#breakpoint()
# Get Data
fname = test_path + "SA17_REPI_" + stat_tag + ".csv"
eqparams_repi = pd.read_csv(fname)
fname_rhypo = test_path + "SA17_RHYPO_" + stat_tag + ".csv"
eqparams_rhypo = pd.read_csv(fname_rhypo)
fname_rjb = test_path + "SA17_RJB_" + stat_tag + ".csv"
eqparams_rjb = pd.read_csv(fname_rjb)

myeqparams = {
    "mag": eqparams_repi["rup_mag"],
    "rjb": eqparams_rjb["dist_rjb"],
    "repi": eqparams_repi["dist_repi"],
    "rhypo": eqparams_rhypo["dist_rhypo"],
    "rake": eqparams_repi["rup_rake"],
    "vs30": eqparams_repi["site_vs30"],
}

# Calculate on my end
my_calc_stats = {}
for dm_i, dist_meas in enumerate(["rjb", "repi", "rhypo"]):
    if dist_meas == "rjb":
        dist_name = "Rjb"
    elif dist_meas == "repi":
        dist_name = "Repi"
    elif dist_meas == "rhypo":
        dist_name = "Rhyp"

    gmm_model = "SandikkayaAkkar2017" + dist_name
    for myimt in [CAV(), IA()]:
        lnmean_pre, sigma, tau, phi, ctx = get_gmm_im.get_gmm_im(
            myeqparams, gmm_model, myimt
        )
        conv_fact = 0
        conv_fact2 = 0

        gfact = 9.81
        ofact = 100

        if str(myimt) == 'CAV':
            # convert g-s to m/s
            conv_fact = -np.log(1 / gfact)
            # convert m/s to cm/s
            conv_fact2 = -np.log(1 / ofact)
        # lnmean = lnmean_pre + conv_fact + conv_fact2
        # mean_alt_units = np.exp(lnmean)
        # mean_alt_same = np.exp(lnmean_pre) * gfact * ofact
        lnmean = lnmean_pre
        mean = np.exp(lnmean_pre)
        # breakpoint()
        my_calc_stats[str(myimt).lower() + "_" + dist_meas + "_lnmean"] = lnmean
        my_calc_stats[str(myimt).lower() + "_" + dist_meas + "_mean"] = mean
        my_calc_stats[str(myimt).lower() + "_" + dist_meas + "_sigma"] = sigma
        my_calc_stats[str(myimt).lower() + "_" + dist_meas + "_tau"] = tau
        my_calc_stats[str(myimt).lower() + "_" + dist_meas + "_phi"] = phi

#breakpoint()
# Do check
for dm_i, dist_meas in enumerate(["rjb", "repi", "rhypo"]):
    for myimt in [CAV(), IA()]:
        if dist_meas == "repi":
            eqparams = copy.deepcopy(eqparams_repi)
        elif dist_meas == "rjb":
            eqparams = copy.deepcopy(eqparams_rjb)
        elif dist_meas == "rhypo":
            eqparams = copy.deepcopy(eqparams_rhypo)

        CBval = eqparams[str(myimt).lower()]
        per_diff = np.abs(
            (
                CBval
                - my_calc_stats[str(myimt).lower() + "_" + dist_meas + "_"
                                + stat_type]
            )
            / CBval
            * 100
        )
        Tval = per_diff < perc_thresholds[stat_type]
        #breakpoint()
        for tval_i, tval in enumerate(Tval):
            total_count += 1
            if not tval:
                fail_count += 1
                print(
                    str(myimt).upper()
                    + " "
                    + stat_type.upper()
                    + " "
                    + dist_meas
                    + " is FAILURE for EQ on line "
                    + str(tval_i + 2)
                    + ": M="
                    + str(eqparams["rup_mag"][tval_i])
                    + ", Rjb="
                    + str(eqparams["dist_" + dist_meas][tval_i])
                    + " km, Vs30 = "
                    + str(eqparams["site_vs30"][tval_i])
                    + " m/s "
                )
                key_name = myimt.string.lower() + "_" + dist_meas + "_" + stat_type
                            
                print(
                    "CBval: "
                    + str(CBval[tval_i])
                    + ", myval: "
                    + str(
                        my_calc_stats[key_name][tval_i]
                    )
                    + ", percent diff is "
                    + "{:.1f}".format(per_diff[tval_i])
                    + "%"
                    + "\n"
                )
                breakpoint()
                print("do investigate")
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
