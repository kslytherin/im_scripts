import pandas as pd
import numpy as np
import compare_gmms
import matplotlib.pyplot as plt
from openquake.hazardlib.imt import CAV, IA

rup_mag_thresh = 7
dist_name = "RuptureDistance"  # "RuptureDistance", "JoynerBooreDistance", "EpicentralDistance", "HypocentralDistance"
pre_dist_meas = "Rrup"
dist_meas = pre_dist_meas.lower()
myimt = CAV()

# take a look at earthquakes to choose from
u_eq_dir = "~/output_conus_test/"
generic_prefix = "EQ_MacedoEtAl2021-BooreEtAl2014_"
geo_reg = "Alaska"  # EUS, WUS, Alaska, CONUS, Globe
tect_reg = "SubductionIntraslab"  #  Active, Stable, Volcanic, SubductionIntraslab, SubductionIntraslab,
U_EQ_file = u_eq_dir + generic_prefix + geo_reg + "_" + tect_reg + "_v2.csv"
u_eq = pd.read_csv(U_EQ_file)
# u_eq = u_eq.sort_values(by="numstns", ascending=False)
u_eq = u_eq.sort_values(by="mag", ascending=False)
# choosing the following earthquakes
# eq_choose = "usp000jwyb"  # M 6.3
# eq_choose = "ci38457511"  # ridgecrest
# eq_choose = "ci38443183"  #
eq_choose = "nn00234425"
eq_choose = "ak0159nc9dk8"
eq_choose = "ak021gbh4rso"  # AK, sub-intraslab
# eq_choose = "ak018fe58vlz"  # AK, sub-intraslab
# eq_choose = "ak018fcnsk91"  # AK, sub-interface, anchorage
# eq_choose = "ak021bt4ffvw"  # AK, active, tanana valley
eq_choose = "nc73827571"  # WUS, interface
# eq_choose = "uw61965081"  # WUS, slab
# mydata1 = mydata1.sort_values(by="EarthquakeTime") # eathquake since 1992 to 2024
# breakpoint()
# eq_choose = "ak20419010"  # Anchorage Earthquake

u_eq_wtect_file = "~/Documents/GMM_data/for_kyle_conus/CONUS_HW_AK_unique_eqs_tect_region_final.csv"
u_eq_wtect = pd.read_csv(u_eq_wtect_file)
u_eq_wtect = u_eq_wtect[u_eq_wtect["EarthquakeId"] == eq_choose]
eq_tect = u_eq_wtect['Tect_region'].values[0]
# eq_geo = u_eq_wtect['Geo_region'].values[0] # wish list
# for active crust 
if eq_tect == "Active":
    main_gmm_set = [
        "CampbellBozorgnia2019",
        ("MacedoEtAl2021", "BooreEtAl2014"),
        "SandikkayaAkkar2017Rjb",
    ]
    nice_gmm_model_str = [
        "Campbell and Bozorgnia (2019)",
        "Macedo et al. (2021) \n[Boore et al. (2014)]",
        "Sandikkaya and Akkar (2017)",
    ]

# for slab
elif eq_tect == "SubductionIntraslab":
    main_gmm_set = [
        "CampbellBozorgnia2019",
        "LiuMacedo2022M1SSlab",
        "LiuMacedo2022M2SSlab",
        ("LiuMacedo2022SSlab", "ParkerEtAl2020SSlab"),
        ("MacedoEtAl2021", "BooreEtAl2014"),
        "SandikkayaAkkar2017Rjb",
    ]
    nice_gmm_model_str = [
        "Campbell and Bozorgnia (2019)",
        "Liu and Macedo (2022) M1 Slab",
        "Liu and Macedo (2022) M2 Slab",
        "Liu and Macedo (2022) Slab \n[Parker et al. (2020)]",
        "Macedo et al. (2021) \n[Boore et al. (2014)]",
        "Sandikkaya and Akkar (2017)",
    ]

# for interface
elif eq_tect == "SubductionInterface":
    main_gmm_set = [
        "CampbellBozorgnia2019",
        "LiuMacedo2022M1SInter",
        "LiuMacedo2022M2SInter",
        ("LiuMacedo2022SInter", "ParkerEtAl2020SInter"),
        ("MacedoEtAl2021", "BooreEtAl2014"),
        "SandikkayaAkkar2017Rjb",
    ]
    nice_gmm_model_str = [
        "Campbell and Bozorgnia (2019)",
        "Liu and Macedo (2022) M1 Inter",
        "Liu and Macedo (2022) M2 Inter",
        "Liu and Macedo (2022) Inter \n[Parker et al. (2020)]",
        "Macedo et al. (2021) \n[Boore et al. (2014)]",
        "Sandikkaya and Akkar (2017)",
    ]
main_col = ["m", "b", "c", "k", "g", "r"]
main_ls = ["-", "-", "-", "--", "--", "-"]


test_path = "~/Documents/GMM_data/for_kyle_conus/"
records_wEQinfo = (
    test_path
    + "CONUS_HW_AK_h1_records_stn_filtered_wEQparams_M{}_thresh.csv".format(
        rup_mag_thresh
    )
)

mydata1 = pd.read_csv(records_wEQinfo)
mydata2 = pd.read_csv(test_path + "CONUS_HW_AK_h1_records_stn_filtered_h2.csv")

# Has to match GMMvsData
z2pt5_assume = 0.6
z1pt0_assume = 0.1
hypo_depth_assume = mydata1["EarthquakeDepth"].values
magstr = "PreferredMagnitude"
vs30str = "Preferred_Vs30"

n = 100
dist_min = 1
dist_max = 500
dist_pts = np.logspace(np.log10(dist_min), np.log10(dist_max), n)

pga_vals = mydata1["PGA"].values
pga_vals2 = mydata2["PGA"].values
pga_gm_dta = np.sqrt(pga_vals * pga_vals2)
cav_gm_dta = np.sqrt(mydata1["CAV"].values * mydata2["CAV"].values)
mydata1["CAV_GM"] = cav_gm_dta

mydata1 = mydata1[mydata1["EarthquakeId"] == eq_choose]


# do the following by earthquake
myeqparams = {
    "mag": mydata1[magstr].values[0],
    "rjb": dist_pts,
    # "repi": mydata1["EpicentralDistance"].values,
    "rhypo": dist_pts,
    "rx": dist_pts,
    # "ry": mydata1["GC2_ry"].values,
    "ry0": dist_pts,
    "ztor": mydata1["EarthquakeZtor"].values[0],
    "rake": mydata1["EarthquakeRake"].values[0],
    "width": mydata1["EarthquakeWidth"].values[0],
    "rrup": dist_pts,
    "dip": mydata1["EarthquakeDip"].values[0],
    "vs30": np.median(mydata1[vs30str].values),
    "alpha": mydata1["alpha"].values[0],
    "hypo_depth": np.median(hypo_depth_assume),
    "z2pt5": z2pt5_assume,
    "z1pt0": z1pt0_assume,
    "vs30measured": True,
}


mean_set, lnmean_set, sigma_set, phi_set, tau_set, ctx, myplots = (
    compare_gmms.compare_gmms(
        myimt,
        main_gmm_set,
        myeqparams,
        dist_meas,
        fig_ops={"makefigs": False, "showfigs": False, "printfigs": False},
    )
)

stat_sets = [mean_set, sigma_set, phi_set, tau_set]

# Scatter plot of cav_gm_dta vs dist_pts with axes plotting
fig, ax = plt.subplots(figsize=(8, 6))
scatter = ax.scatter(
    mydata1["RuptureDistance"],
    mydata1["CAV_GM"],
    alpha=0.7,
    edgecolors="k",
    label="data",
)
for m_i, mean_line in enumerate(mean_set):
    ax.loglog(
        ctx[dist_meas],
        mean_set[m_i],
        linestyle=main_ls[m_i],
        color=main_col[m_i],
        linewidth=2,
        label=nice_gmm_model_str[m_i],
    )
ax.legend(loc="best", fontsize=10, frameon=True)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Distance to Earthquake, km")
ax.set_ylabel("CAV (g-s)")
ax.set_xlim(dist_min, dist_max)
ax.set_title(
    f"{geo_reg} - {eq_tect}: {eq_choose}, M{mydata1[magstr].values[0]}, Vs30={int(myeqparams['vs30'])} m/s"
)
plt.tight_layout()
plt.show()


breakpoint()
print("end")
