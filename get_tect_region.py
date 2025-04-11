import pandas as pd

# Get tectonic region information
test_path = "~/Documents/GMM_data/for_kyle_conus/"
tect_region = "CONUS_HW_AK_unique_eqs_tect_region.csv"
fulltec = test_path + tect_region
EQ_w_tec_info = pd.read_csv(fulltec)
tec_prob_list = [
    "ProbabilityActive",
    "ProbabilityStable",
    "ProbabilitySubduction",
    "ProbabilityVolcanic",
]

# Get maximum probability and associated name
EQ_w_tec_info["MaxProbability"] = EQ_w_tec_info[tec_prob_list].max(axis=1)
EQ_w_tec_info["MaxProbabilityName"] = EQ_w_tec_info[tec_prob_list].idxmax(axis=1)
# Check if MaxProbabilityName is "ProbabilitySubduction"
subduction_mask = EQ_w_tec_info["MaxProbabilityName"] == "ProbabilitySubduction"

# Find the max of the subduction probabilities
subduction_probs = EQ_w_tec_info.loc[
    subduction_mask,
    [
        "ProbabilitySubductionCrustal",
        "ProbabilitySubductionInterface",
        "ProbabilitySubductionIntraslab",
    ],
]

EQ_w_tec_info.loc[subduction_mask, "MaxSubductionProbability"] = subduction_probs.max(
    axis=1
)
EQ_w_tec_info.loc[subduction_mask, "MaxSubductionProbabilityName"] = (
    subduction_probs.idxmax(axis=1)
)

EQ_w_tec_info["FinalProbabilityName"] = EQ_w_tec_info[
    "MaxSubductionProbabilityName"
].fillna(EQ_w_tec_info["MaxProbabilityName"])

EQ_w_tec_info["FinalProbability"] = EQ_w_tec_info["MaxSubductionProbability"].fillna(
    EQ_w_tec_info["MaxProbability"]
)
# Remove "Probability" string and rename "SubductionIntraslab" to "Active"
EQ_w_tec_info["Tect_region"] = EQ_w_tec_info["FinalProbabilityName"].str.replace(
    "Probability", ""
)

EQ_w_tec_info["Tect_region"] = EQ_w_tec_info["Tect_region"].replace(
    "SubductionCrustal", "Active"
)
output_file = test_path + "CONUS_HW_AK_unique_eqs_tect_region_final.csv"
EQ_w_tec_info.to_csv(output_file, index=False)
