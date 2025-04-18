import numpy as np
import mytools
from scipy.constants import g as gfact

"""
assuming already in log space
first assuming that models have m/s units but converting to g-s
figure out how to do this better, for now it might only work for CAV
"""


def convert_GMM_units(single_gmm_inst, lnmean):

    if isinstance(single_gmm_inst, tuple):
        main_model = single_gmm_inst[0]
    else:
        main_model = single_gmm_inst
    print("original name is " + main_model)
    model_split_name = mytools.split_num_letters(main_model)
    model_shorter_name = "".join(model_split_name)
    print("joined name is " + model_shorter_name)
    # if model_shorter_name in ["CampbellBozorgnia2019", "LiuMacedo2022"]:
    if model_shorter_name in ["CampbellBozorgnia2019"]:
        print("converting for " + main_model)
        # convert from m/s to g-s
        conv_fact = np.log(1 / gfact)
    else:
        conv_fact = 0
    conv_lnmean = lnmean + conv_fact
    return conv_lnmean
