import sys
import os.path
import matplotlib.pyplot as plt
import seaborn as sns
import mytools

import pandas as pd
import numpy as np

printfigs = True
make_tick_labels = True
ver_no = "v2"

tick_tag = ""
if make_tick_labels:
    tick_tag = "_wticklabels"

selected_tect_reg = [
    "Active",
    "SubductionInterface",
    "SubductionIntraslab",
]

# for geo_reg in ["Globe", "Hawaii", "Alaska", "EUS", "WUS"]:
for geo_reg in ["Alaska", "WUS"]:
    print("{}".format(geo_reg))
    for tect_reg in selected_tect_reg:
        print("{}".format(tect_reg))
        output_csv = (
            "/Users/kksmith/data/GMMout/Test_Results_"
            + geo_reg
            + "_"
            + tect_reg
            + "_"
            + ver_no
            + ".csv"
        )
        if os.path.exists(output_csv):
            df = pd.read_csv(output_csv, index_col=0)
            print("reading " + output_csv)
            # print(df.head())
        else:
            print(f"File {output_csv} does not exist.")
            continue

        df.rename(columns={"Mak": "MLLH", "LLH": "ULLH"}, inplace=True)
        df = df[
            [
                "CHISQ-MF",
                "ULLH",
                "MLLH",
                "MDE_norm",
                "kappa_sq",
                "EDR",
                "AM",
            ]
        ]
        df = df.drop(columns=["MDE_norm", "kappa_sq"])
        # print(df.index)
        # NL_split = mytools.split_num_letters(df.index[0])
        # for nl in NL_split:
        #     print(mytools.split_cap_words_crop(nl))

        if not df.empty:
            plt.figure(figsize=(10, 8))
            for i, column in enumerate(df.columns):
                plt.subplot(1, len(df.columns), i + 1)
                sns.heatmap(
                    df[[column]], annot=True, fmt=".3f", cmap="Blues_r", cbar=False
                )
                plt.gca().xaxis.set_ticks_position("none")
                plt.gca().yaxis.set_ticks_position("none")
                plt.xticks(fontsize=14)
                plt.ylabel("")
                plt.yticks([])  # for NRC plots
                if not make_tick_labels:
                    plt.xticks([])
                    plt.yticks([])
                if i > 0:
                    plt.yticks([])
                if i == 3:
                    plt.title(f"{geo_reg} - {tect_reg}")
                plt.title("")  # for NRC plots
            plt.tight_layout()
            fig_file_name = f"/Users/kksmith/figures/heatmaps/heatmap_{geo_reg}_{tect_reg}_{ver_no}{tick_tag}.png"
            print(fig_file_name)
            if printfigs:
                plt.savefig(fig_file_name)
            # plt.show()
            # breakpoint()
            # plt.close()
