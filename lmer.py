"""
slightly modified version of James' script:
https://ghsc.code-pages.usgs.gov/esi/groundmotion-processing/contents/tutorials/lme.html
added in plt.show() commands and change data_path
copied on Mar 12, 2024
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

from openquake.hazardlib import valid
from openquake.hazardlib.imt import SA
from openquake.hazardlib.imt import PGA
from openquake.hazardlib.imt import PGV
from openquake.hazardlib.contexts import simple_cmaker
from gmprocess.utils.constants import DATA_DIR

# GMM_model = "AbrahamsonEtAl2014RegTWN"
GMM_model = "CampbellBozorgnia2014"
# imt = SA(1.0)
imt = PGV()
# breakpoint()
# will need to somehow use the previous info to build the following instead of doing it "by hand" each time
# IMT_type = "SA(1.000)"
IMT_type = "PGV"
resid_col_name = "SA_1_ASK14_residuals"
data_column_name = "SA_rotd50.0_2020.03.31.csv"
show_plt = True

# Path to example data
# source_path = "/Users/kksmith/src/groundmotion-processing/" # where example data is found
# data_path = source_path / DATA_DIR / "lme" / data_column_name
data_path = (
    "/Users/kksmith/gm_projects/testing4/data/New_Proj_default_metrics_rotd50.0.csv"
)

# read data and set up gmm
df = pd.read_csv(data_path)
gsim = valid.gsim("[" + GMM_model + "]")

predicted_data = []

cmaker = simple_cmaker(
    [gsim], [str(imt)], mags=["%2.f" % mag for mag in df.EarthquakeMagnitude]
)
n = df.shape[0]
ctx = cmaker.new_ctx(n)
ctx["mag"] = df.EarthquakeMagnitude
ctx["dip"] = 90.0
ctx["rake"] = 0.0
ctx["width"] = 10 ** (-0.76 + 0.27 * df.EarthquakeMagnitude)
ctx["ztor"] = df.EarthquakeDepth
ctx["vs30"] = 760.0
# ctx["vs30measured"] = False ## won't work with CB14
ctx["rrup"] = df.RuptureDistance
ctx["rjb"] = df.JoynerBooreDistance
# ctx["z1pt0"] = 48.0 ## won't work with CB14
ctx["rx"] = -1
# ctx["ry0"] = -1 ## won't work with CB14

# CB14 needs this
ctx["z2pt5"] = 10000
ctx["hypo_depth"] = df.EarthquakeDepth

# Evaluate the GMM.
mean = cmaker.get_mean_stds([ctx])[0][0][0]
# Convert from ln(g) to %g
predicted_data = 100.0 * np.exp(mean)

df[resid_col_name] = np.log(df[IMT_type]) - np.log(predicted_data)

mdf = smf.mixedlm("%s ~ 1" % resid_col_name, df, groups=df["EarthquakeId"]).fit()
print(mdf.summary())


### Investigative plots
# IM vs Distance
fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
axes.scatter(df["RuptureDistance"], df[IMT_type])
axes.set_xlabel("$R_{rup}$")
# axes.set_xscale('log')
# axes.set_xscale('log') good to have for a lot of data that is well spread out
axes.set_ylabel(IMT_type)
# axes.axhline(0, ls='--', c='k')
# axes.set_ylim(-2.5, 2.5);
axes.set_title("N = " + str(len(df_events["Group"])))
# if show_plt:
#    plt.show()

# Mag vs Distance
fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
axes.scatter(df["RuptureDistance"], df["EarthquakeMagnitude"])
axes.set_xlabel("$R_{rup}$")
axes.set_ylabel("Mag")
# axes.set_xscale('log') good to have for a lot of data that is well spread out
# if show_plt:
#    plt.show()

### Regression and residual analysis
# plotting random effects of between event terms?
re_array = [float(re.iloc[0]) for group, re in mdf.random_effects.items()]
fig = plt.figure()
plt.plot(re_array, "o")
plt.xlabel("Index")
plt.ylabel("Random effect")
plt.title("N = " + str(len(re_array)))
# if show_plt:
#    plt.show()

# plotting within-event residuals
fig = plt.figure()
plt.plot(df["EpicentralDistance"], mdf.resid, "o")
plt.xlabel("Epicentral Distance, km")
plt.ylabel("Within-event residual")
plt.title("N = " + str(len(mdf.resid)))
# if show_plt:
#    plt.show()


# plotting between-event terms?
btw_event_terms = pd.DataFrame(mdf.random_effects).T
df = df.merge(btw_event_terms, left_on="EarthquakeId", right_index=True)
df_events = df.drop_duplicates(subset=["EarthquakeId"])


fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
axes.scatter(df_events["EarthquakeDepth"], df_events["Group"])
axes.set_xlabel("Hypocentral Depth (km)")
axes.set_ylabel("Event term (T = 1 s)")
axes.axhline(0, ls="--", c="k")
axes.set_ylim(-2.5, 2.5)
axes.set_title("N = " + str(len(df_events["Group"])))
if show_plt:
    plt.show()


breakpoint()
