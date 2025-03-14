import pandas as pd
import numpy as np
'''
get best magnitudes for records and write out file of all EQs without duplicates
'''
out_dir = "/Users/kksmith/Documents/GMM_data/for_kyle_conus/"
# file with records
infile = "~/im_scripts/data/combine_files_here/dsgmd_h1_wVs30.csv"
stn_dta = pd.read_csv(infile)
infile_eid = 'EarthquakeId'

# EQ file with better mags
qmlfile = "~/im_scripts/data/quakeml/CONUS_AK_HW_EQs_latlon_box.csv"
# "EventID" or "PreferredOriginID"?
evid_col = "EventID"
quake_dta_all = pd.read_csv(qmlfile)
quake_dta = quake_dta_all.loc[:, ["EventID", "EventTime", "EventLatitude", "EventLongitude", \
                        "EventDepth", "EventType", "PreferredMagnitude", \
                        "PreferredMagnitudeSource", "PreferredMagnitudeMethod", "PreferredOriginID"]]
# final_dta = stn_dta.set_index('EarthquakeId').join(quake_dta.set_index(evid_col))

# Combined file of records with best mags
final_dta = stn_dta.join(quake_dta.set_index(evid_col), on=infile_eid)
final_dta.to_csv(out_dir + 'CONUS_HW_AK_h1_records.csv', index=False)

# make file of just earthquakes
eq, eq_i = np.unique(stn_dta[infile_eid], return_index=True)
stn_dta.iloc[eq_i].to_csv(out_dir + 'CONUS_HW_AK_unique_eqs.csv', index=False)

# make file of just stations
stn, stn_i = np.unique(stn_dta["StationID"], return_index=True)
stn_dta.iloc[stn_i].to_csv(out_dir + 'CONUS_HW_AK_unique_stns.csv', index=False)
breakpoint()
print("end")
