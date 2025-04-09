import pandas as pd
import matplotlib.pyplot as plt
from openquake.hazardlib.imt import PGA, SA
from mytools import *
import get_gmm_im
import scipy.constants

con_model_set = ["ParkerEtAl2020SInter", "KuehnEtAl2020SInter"]
dist_name = "Rrup"
dist_meas = dist_name.lower()
for cmodel in con_model_set:
    auth = split_cap_words_crop(cmodel)[0]
    tect = split_cap_words_crop(cmodel)[-1]

    # from matlab files
    datafile = "/Users/kksmith/Downloads/parker_and_nico/"
    # full_dtafile_name = datafile + "pre_im_values_" + auth + "_" + tect + ".csv"
    full_dtafile_name = datafile + "ln_pre_im_values_" + auth + "_" + tect + ".csv"
    data = pd.read_csv(full_dtafile_name)

    myeqparams = {
        # "mag": eqmagset,
        # "vs30": vs30set,  # from text
        # "ztor": ztorset,
        # "hypo_depth": ztorset,
        # dist_meas: dist_pts,
        "mag": data["M"].values,
        "vs30": data["Vs30"].values,  # from text
        "ztor": data["hyperD"].values,
        "hypo_depth": data["hyperD"].values,
        "PGA_MEAN": data["PGAmean"].values,
        "SA(1.0)_MEAN": data["SA1p0mean"].values,
        dist_meas: data["Rrup"].values,
        "PGA_TOTAL_STDDEV": data["PGAsigma"].values,
        "SA(1.0)_TOTAL_STDDEV": data["SA1p0sigma"].values,
    }

    imlist = [PGA(), SA(1.0)]
    for im_i, im in enumerate(imlist):
        just_im = "".join(split_cap_words_crop(str(im)))
        if just_im == "SA":
            im_lab = just_im + "1p0mean"
        else:
            im_lab = just_im + "mean"
        (
            lnmedian_pre_C,
            sigma_C,
            tau_C,
            phi_C,
            ctx_C,
        ) = get_gmm_im.get_gmm_im(myeqparams, cmodel, im)

        gfact = scipy.constants.g

        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(
            data["Rrup"].values,
            # data[im_lab].values - 1,
            # "b--",
            # label="mlab",
            # np.log(data[im_lab].values - 1),
            # np.log(np.exp(data[im_lab].values) - 1),
            np.exp(data[im_lab].values) - 1,
            # data[im_lab].values,
            "b--",
            label="mlab",
        )
        ax.plot(data["Rrup"].values, np.exp(lnmedian_pre_C), "r--", label="oq")
        # ax.plot(data["Rrup"].values, lnmedian_pre_C, "r--", label="oq")
        ax.set_xlabel("Rrup, km")
        ax.set_title(cmodel + " " + im_lab)
        ax.legend()
        plt.show()
        # print(np.log(np.exp(data[im_lab].values) - 1) - lnmedian_pre_C)
        print(np.exp(data[im_lab].values) - np.exp(lnmedian_pre_C))
        print(np.mean(np.exp(data[im_lab].values) - np.exp(lnmedian_pre_C)))
        breakpoint()
