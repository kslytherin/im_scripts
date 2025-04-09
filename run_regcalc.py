import strec.utils as su
from strec.calc import select_regions


eq_dir = "~/Documents/GMM_data/for_kyle_conus/"

input_file = eq_dir + "CONUS_HW_AK_unique_eqs_rename_short.csv"
# input_file = eq_dir + "CONUS_HW_AK_unique_eqs_rename.csv"
hypo_columns = ["latitude", "longitude", "depth", "magnitude"]
id_column = "EarthquakeId"
yay, idcol, latcol, loncol, depcol, magcol = su.read_input_file(
    input_file, hypo_columns, id_column
)

regions_frame = select_regions(
    args.input_file,
    args.eqinfo,
    args.moment_info,
    args.event_id,
    args.hypo_columns,
    args.id_column,
    args.verbose,
)
