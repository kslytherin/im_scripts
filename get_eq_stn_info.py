import numpy as np

from esi_utils_rupture.origin import Origin
from esi_utils_rupture.factory import get_rupture
import json
from esi_utils_time.ancient_time import HistoricTime

"""
Input: QuadRupture, json rupture file from usgs
Output: distances
"""


def get_event_from_rup(rup_file):
    with open(rup_file, "r") as file:
        eq_json = json.load(file)
    json_event = eq_json["metadata"]
    rup_type = eq_json["features"][0]["geometry"]["type"]

    for js_key in ["locstring", "network", "netid", "productcode"]:
        json_event[js_key] = ""

    origin = Origin(json_event)

    return origin, json_event, rup_type


def get_eq_stn_dist_vals(origin, rup_file, stn_locs):
    rup = get_rupture(origin, rup_file)
    stn_lon = stn_locs["lon"]
    stn_lat = stn_locs["lat"]
    depth = np.ones_like(stn_lat)

    rjb = rup.computeRjb(stn_lon, stn_lat, depth)[0]
    rhyp = rup.computeRhyp(stn_lon, stn_lat, depth)
    rrup = rup.computeRrup(stn_lon, stn_lat, depth)[0]
    repi = rup.computeRepi(stn_lon, stn_lat, depth)
    o_dis = rup.computeGC2(stn_lon, stn_lat, depth)
    width = rup.getWidth()
    ztor = rup.getDepthToTop()
    dip = rup.getDip()
    return rjb, rhyp, rrup, repi, o_dis, width, ztor, dip


def make_stn_array(lonlat_center, lonlatbox, arraydense):
    # makes a rectangular grid of locations using the center of the grid (lonlat_center) and deviation from it w/lonlatbox w/density arraydense
    # lonlatbox = [lonmin, lonmax, latmin, latmax]
    stn_lon_pre = lonlat_center["lon"] + np.linspace(
        lonlatbox[0], lonlatbox[1], arraydense["nlon"]
    )
    stn_lat_pre = lonlat_center["lat"] + np.linspace(
        lonlatbox[2], lonlatbox[3], arraydense["nlat"]
    )

    stn_lon_set, stn_lat_set = np.meshgrid(stn_lon_pre, stn_lat_pre)

    stn_lon = np.concatenate(stn_lon_set)
    stn_lat = np.concatenate(stn_lat_set)
    stn_locs = {"lon": stn_lon, "lat": stn_lat}
    # breakpoint()
    # return stn_lon, stn_lat
    return stn_locs
