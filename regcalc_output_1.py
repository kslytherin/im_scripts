#!/usr/bin/env python

# mike says it will be fixed but I can call this script like this:
# $python regcalc_output.py ~/tmp/CONUS_HW_AK_unique_eqs_rename_short_simple.csv output.csv

import sys

import pandas as pd
from strec.calc import select_regions

if __name__ == "__main__":
    input_file = sys.argv[1]  # a CSV file
    output_file = sys.argv[2]  # desired output csv file
    regions_frame = select_regions(
        input_file,  # eqinfo
        None,  # eqinfo
        None,  # moment_info
        None,  # event_id
        ("EarthquakeLatitude","EarthquakeLongitude","EarthquakeDepth","EarthquakeMagnitude"),
        None,  # id_column
        None,  # verbose
    )
#        None,  # hypo_columns ()
    regions_frame.to_csv(output_file, index=False)
    print(f"{len(regions_frame)} records written to {output_file}.")
